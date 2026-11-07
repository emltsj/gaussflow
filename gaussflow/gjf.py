"""Gaussian input (.gjf) generation with a controlled route whitelist.

The Agent never hand-writes a route line into a shell.  A route is split into
a *base* (method/basis + tolerated extra keywords) and the *stage* keyword
(opt / freq / none for sp) that the pipeline appends.  Unknown methods, basis
sets or keywords are flagged and can only be submitted with explicit
approval, and resource allocations come from the central config.
"""
from __future__ import annotations

import re

from .structure import formula_of_coords

STAGE_KEYWORDS = {"opt", "freq", "force", "stable", "scan", "irc", "sp"}

# extra route keywords that are considered safe to pass through
_EXTRA_OK = re.compile(
    r"^(int|grid|pop|density|nosymm|symmetry|scf|maxdisk|guess|geom|6d|10f|"
    r"temperature|pressure|td|cis|nstates|root|nmr|polar|volume|scrf|smd|"
    r"empiricaldispersion|iop|freq|opt|force|stable|scan|irc)\b", re.I)

_SEMIEMPIRICAL = {"PM6", "PM7", "PM3", "AM1", "AM1D", "CNDO", "INDO", "MNDO",
                  "PDDG", "ZINDO"}


def split_route(route: str) -> list[str]:
    """Split a route string into whitespace-separated tokens (no '#')."""
    route = (route or "").strip()
    route = re.sub(r"^#+\s*[pPnNtT]?\s*", "", route)
    return [t for t in re.split(r"\s+", route) if t]


def strip_stage_keywords(tokens: list[str]) -> tuple[list[str], list[str]]:
    base, removed = [], []
    for t in tokens:
        low = t.lower()
        if low in STAGE_KEYWORDS or low.startswith(("opt=", "freq=")):
            removed.append(t)
        else:
            base.append(t)
    return base, removed


# Gaussian 16 is case-insensitive for methods but QPErr-fails on the
# hyphenated spellings that many tools/LLMs use.
_ROUTE_ALIASES = {
    "DEF2-TZVP": "def2TZVP",
    "DEF2-TZVPP": "def2TZVPP",
    "DEF2-SVP": "def2SVP",
    "DEF2-SVPP": "def2SVPP",
    "DEF2-QZVP": "def2QZVP",
    "M06-2X": "M062X",
    "M06-2X-D3": "M062X-D3",
    "M06-L": "M06L",
}


def normalize_route(route: str) -> str:
    """Normalise known Gaussian spelling traps (def2-TZVP -> def2TZVP, …)."""
    out = []
    for token in split_route(route):
        if "/" in token:
            method, basis = token.split("/", 1)
            method = _ROUTE_ALIASES.get(method.upper(), method)
            basis = _ROUTE_ALIASES.get(basis.upper(), basis)
            out.append(f"{method}/{basis}")
        else:
            out.append(_ROUTE_ALIASES.get(token.upper(), token))
    return " ".join(out)


def _method_basis(token: str) -> tuple[str, str]:
    if "/" in token:
        method, basis = token.split("/", 1)
        return method, basis
    return token, ""


def validate_route(base_tokens: list[str], cfg) -> dict:
    """Return {'ok': bool, 'problems': [...], 'warnings': [...], 'needs_approval': bool}."""
    problems: list[str] = []
    warnings: list[str] = []
    approval_cfg = cfg.approval

    if not base_tokens:
        return {"ok": False, "problems": ["路由为空"], "warnings": [],
                "needs_approval": False}

    method, basis = _method_basis(base_tokens[0])
    allowed_methods = {m.upper().replace("-", "") for m in approval_cfg["allowed_methods"]}
    if method.upper().replace("-", "") not in allowed_methods:
        problems.append(f"方法 '{method}' 不在允许列表（approval.allowed_methods）")

    if not basis:
        if method.upper() not in _SEMIEMPIRICAL:
            problems.append(f"路由 '{base_tokens[0]}' 缺少基组（例如 B3LYP/6-31G(d)）")
    else:
        allowed_basis = {b.upper() for b in approval_cfg["allowed_basis"]}
        if basis.upper() not in allowed_basis:
            warnings.append(f"基组 '{basis}' 不在允许列表（需人工确认）")

    for tok in base_tokens[1:]:
        if not _EXTRA_OK.match(tok):
            warnings.append(f"未知路由关键字 '{tok}'（需人工确认）")

    needs_approval = bool(warnings)
    upper_route = " ".join(base_tokens).upper()
    if any(k.upper() in upper_route for k in approval_cfg["expensive_keywords"]):
        needs_approval = True
    return {"ok": not problems, "problems": problems, "warnings": warnings,
            "needs_approval": needs_approval, "method": method, "basis": basis}


def compose_route(base_tokens: list[str], stage: str,
                  solvent: str | None = None,
                  extra_keywords: list[str] | None = None,
                  retry_patch: str | None = None) -> str:
    """Build the final '#p ...' route for a stage."""
    tokens = list(base_tokens)
    stage_patch = None
    if retry_patch and retry_patch.lower().startswith(("opt=", "freq=")):
        # e.g. "opt=(maxcycles=200,calcfc)" replaces the plain stage keyword
        stage_patch = retry_patch
    if stage in ("opt", "preopt"):
        tokens.append(stage_patch or "opt")
    elif stage == "freq":
        tokens.append(stage_patch or "freq")
    elif stage == "sp":
        pass
    else:
        raise ValueError(f"未知 stage: {stage}")
    if solvent:
        tokens.append(f"scrf=(smd,solvent={solvent})")
    if extra_keywords:
        tokens.extend(extra_keywords)
    if retry_patch and not stage_patch:
        tokens.append(retry_patch)
    return "#p " + " ".join(tokens)


def build_gjf(
    *,
    route_line: str,
    title: str,
    charge: int,
    multiplicity: int,
    coords: list[tuple[str, float, float, float]] | None = None,
    chk: str = "job.chk",
    oldchk: str | None = None,
    mem: str = "8GB",
    nproc: int = 8,
    geom_check: bool = False,
    comment: str | None = None,
) -> str:
    """Render a complete .gjf text."""
    lines: list[str] = []
    lines.append(f"%chk={chk}")
    if oldchk:
        lines.append(f"%oldchk={oldchk}")
    lines.append(f"%mem={mem}")
    lines.append(f"%nproc={nproc}")
    lines.append("")
    lines.append(route_line)
    lines.append("")
    lines.append(title or "gaussflow job")
    if comment:
        lines.append(comment)
    lines.append("")
    lines.append(f"{int(charge)} {int(multiplicity)}")
    if geom_check:
        # geometry and wavefunction are read from %oldchk
        lines.append("")
    else:
        if not coords:
            raise ValueError("非 geom=check 任务必须提供坐标")
        for s, x, y, z in coords:
            lines.append(f"{s:<2s} {x:14.6f} {y:14.6f} {z:14.6f}")
        lines.append("")
    return "\n".join(lines) + "\n"


def describe_input(coords, charge: int, multiplicity: int, route_line: str) -> dict:
    return {
        "natoms": len(coords) if coords else None,
        "formula": formula_of_coords(coords) if coords else None,
        "charge": charge,
        "multiplicity": multiplicity,
        "route": route_line,
        "elements": sorted({s for s, *_ in coords}) if coords else [],
    }


def build_counterpoise_gjf(
    *,
    route_line: str,
    title: str,
    coords: list[tuple[str, float, float, float]],
    fragments: list[dict],
    total_charge: int,
    total_multiplicity: int,
    chk: str = "job.chk",
    mem: str = "8GB",
    nproc: int = 8,
) -> str:
    """Counterpoise (BSSE) single point for a complex.

    Verified against Gaussian 16W on this machine: the charge/multiplicity line
    is **total first, then one q/m pair per fragment**, e.g. for Li⁺···water:
    ``1 1 0 1 1 1`` (total +1/1, fragment1 0/1, fragment2 +1/1).  Every atom
    carries a ``Fragment=N`` attribute.
    """
    total = sum(int(f["n"]) for f in fragments)
    if total != len(coords):
        raise ValueError(f"fragment atom count {total} != coordinates {len(coords)}")
    lines = [f"%chk={chk}", f"%mem={mem}", f"%nproc={nproc}", "",
             route_line, "", title, ""]
    qm: list[str] = [str(int(total_charge)), str(int(total_multiplicity))]
    for frag in fragments:
        qm += [str(int(frag["q"])), str(int(frag["m"]))]
    lines.append(" ".join(qm))
    idx = 0
    for fi, frag in enumerate(fragments, start=1):
        for _ in range(int(frag["n"])):
            s, x, y, z = coords[idx]
            idx += 1
            lines.append(f"{s}(Fragment={fi}) {x:14.6f} {y:14.6f} {z:14.6f}")
    lines.append("")
    return "\n".join(lines) + "\n"
