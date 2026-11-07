"""Deterministic F-substitution site enumeration (symmetry-unique).

Answers the question "F 加在哪里": instead of letting the LLM guess positions,
this module replaces H atoms by F systematically and de-duplicates symmetry-
equivalent results via canonical SMILES.

For a molecule such as ethanol ``CCO`` the three symmetry-unique mono-fluoro
isomers are F on the terminal CH3, F on the CH2, and F on the OH.  For an
ethyl group attached to a larger skeleton (``CH3-CH2-R``) the same logic finds
both carbon positions and reports them as ``C1 (CH3 端基)`` / ``C2 (CH2)``, so
the GUI can show exactly where each F went.  Multi-F enumeration prioritises
per-site (e.g. CF3) products first, then mixed sites.
"""
from __future__ import annotations

import itertools
import random
import re

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def _canonical(mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def _formula(mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def _site_type(atom) -> str:
    """Human-readable local environment, e.g. ``CH3 端基`` / ``CH2`` / ``OH``."""
    symbol = atom.GetSymbol()
    hs = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
    if atom.GetIsAromatic():
        return f"芳环 {symbol}H" if hs else f"芳环 {symbol}"
    if symbol == "C":
        return {3: "CH3 端基", 2: "CH2", 1: "CH", 0: "C"}.get(hs, f"C(H{hs})")
    if symbol == "O":
        return "OH" if hs else "O"
    if symbol == "N":
        return {2: "NH2", 1: "NH"}.get(hs, "N")
    if symbol == "S":
        return "SH" if hs else "S"
    if symbol == "P":
        return "PH" if hs else "P"
    return symbol + (f"H{hs}" if hs else "")


def _site_desc(atom) -> str:
    heavy = sorted({n.GetSymbol() for n in atom.GetNeighbors()
                    if n.GetAtomicNum() > 1})
    kind = _site_type(atom)
    return f"{kind}，连 {'、'.join(heavy) or '—'}"


def _h_sites(mol_h) -> list[dict]:
    """All H atoms with the heavy atom they are attached to."""
    out = []
    for atom in mol_h.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        neighbours = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        if not neighbours:
            continue
        heavy = neighbours[0]
        out.append({
            "h_index": atom.GetIdx(),
            "atom_index": heavy.GetIdx(),
            "site": f"{heavy.GetSymbol()}{heavy.GetIdx() + 1}",
            "site_type": _site_type(heavy),
            "desc": _site_desc(heavy),
        })
    return out


def _site_order(site: str):
    m = re.fullmatch(r"([A-Za-z]+)(\d+)", site or "")
    if m:
        return (m.group(1), int(m.group(2)), "")
    return (site or "", 0, site or "")


def _combo_label(combo, by_h: dict[int, dict]) -> tuple[str, list[str], dict]:
    """Return (label, unique sites, site→count) for one F combination.

    Repeated substitutions on the same atom are collapsed into ``C1×2`` so
    labels stay readable (e.g. ``F 取代在 C1×2 + O3``).
    """
    counts: dict[str, int] = {}
    desc: dict[str, str] = {}
    for h in combo:
        info = by_h[h]
        counts[info["site"]] = counts.get(info["site"], 0) + 1
        desc[info["site"]] = info["desc"]
    sites = sorted(counts, key=_site_order)
    parts = []
    for s in sites:
        kind = desc[s].split("，")[0]
        atom = f"{s}×{counts[s]}" if counts[s] > 1 else s
        parts.append(f"{atom}（{kind}）")
    label = "F 取代在 " + " + ".join(parts)
    return label, sites, counts


def describe_fluorination_sites(smiles: str, max_sites: int = 40) -> list[dict]:
    """Unique mono-fluorination positions of ``smiles`` (symmetry-de-duplicated).

    Returns a list of dicts::

        {"site": "C1", "atom_index": 0, "site_type": "CH3 端基",
         "smiles": "<canonical>", "label": "F 取代在 C1（CH3 端基…）",
         "n_F": 1, "formula": "C2H5FO"}
    """
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smiles!r}")
    mol_h = Chem.AddHs(mol)
    seen: dict[str, dict] = {}
    for info in _h_sites(mol_h):
        rw = Chem.RWMol(mol_h)
        atom = rw.GetAtomWithIdx(info["h_index"])
        atom.SetAtomicNum(9)
        atom.SetNoImplicit(True)
        try:
            candidate = rw.GetMol()
            Chem.SanitizeMol(candidate)
        except Exception:
            continue
        can = _canonical(candidate)
        if can in seen:
            continue
        seen[can] = {
            "site": info["site"],
            "atom_index": info["atom_index"],
            "site_type": info["site_type"],
            "smiles": can,
            "label": f"F 取代在 {info['site']}（{info['desc']}）",
            "n_F": 1,
            "formula": _formula(Chem.MolFromSmiles(can)),
        }
        if len(seen) >= max_sites:
            break
    return sorted(seen.values(), key=lambda v: _site_order(v["site"]))


def enumerate_fluorinated_variants(smiles: str, n_f: int = 1,
                                   max_variants: int = 20,
                                   seed: int = 42) -> list[dict]:
    """All symmetry-unique ``n_f``-fold fluorination products of ``smiles``.

    For ``n_f == 1`` the enumeration is complete (every symmetry-unique H
    position).  For ``n_f > 1`` products with all F on the same site (e.g.
    CF3) are considered first, then mixed sites; at most ``max_variants``
    unique structures are returned.  The selection is deterministic.
    """
    if n_f < 1:
        raise ValueError("n_f 必须 ≥ 1")
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smiles!r}")
    mol_h = Chem.AddHs(mol)
    h_info = _h_sites(mol_h)
    if len(h_info) < n_f:
        raise ValueError(f"分子只有 {len(h_info)} 个 H，无法替换 {n_f} 个 F")
    by_h = {info["h_index"]: info for info in h_info}

    def make(combo) -> dict | None:
        rw = Chem.RWMol(mol_h)
        for h in combo:
            a = rw.GetAtomWithIdx(h)
            a.SetAtomicNum(9)
            a.SetNoImplicit(True)
        try:
            candidate = rw.GetMol()
            Chem.SanitizeMol(candidate)
        except Exception:
            return None
        can = _canonical(candidate)
        if not can:
            return None
        label, sites, counts = _combo_label(combo, by_h)
        return {
            "smiles": can,
            "n_F": n_f,
            "sites": sites,
            "site_counts": counts,
            "label": label,
            "formula": _formula(Chem.MolFromSmiles(can)),
        }

    seen: dict[str, dict] = {}
    # A single H site can only accept n_f F atoms when it has n_f H atoms.
    if n_f == 1:
        combos = ((info["h_index"],) for info in h_info)
        examine_limit = None
    else:
        def _all_combos():
            # 1) per-site products first: all F on the same carbon/heteroatom
            for info in h_info:
                same = [x["h_index"] for x in h_info
                        if x["atom_index"] == info["atom_index"]]
                if len(same) >= n_f:
                    for combo in itertools.combinations(same, n_f):
                        yield combo
            # 2) then mixed sites in a deterministic (seeded) order
            order = [x["h_index"] for x in h_info]
            random.Random(seed).shuffle(order)
            for combo in itertools.combinations(order, n_f):
                yield combo
        combos = _all_combos()
        examine_limit = 200_000
    examined = 0
    for combo in combos:
        examined += 1
        if examine_limit is not None and examined > examine_limit:
            break
        product = make(combo)
        if product is None:
            continue
        if product["smiles"] in seen:
            continue
        seen[product["smiles"]] = product
        if len(seen) >= max_variants:
            break
    out = sorted(seen.values(), key=lambda v: v["smiles"])
    for i, v in enumerate(out, 1):
        v["variant_index"] = i
        v["name_suffix"] = f"F{n_f}-{i}"
    return out


def molecules_from_variants(base_name: str, base_smiles: str, variants: list[dict],
                            charge: int = 0, multiplicity: int = 1) -> list[dict]:
    """Convert enumerated variants into the molecule dicts used by the pipeline."""
    base = Chem.MolFromSmiles(base_smiles)
    base_formula = rdMolDescriptors.CalcMolFormula(base) if base else ""
    out = []
    for v in variants:
        mol = Chem.MolFromSmiles(v["smiles"])
        if mol is None:
            continue
        out.append({
            "name": f"{base_name}-{v['name_suffix']}",
            "smiles": v["smiles"],
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "n_F": v["n_F"],
            "natoms": mol.GetNumAtoms() + sum(a.GetTotalNumHs() for a in mol.GetAtoms()),
            "charge": int(charge),
            "multiplicity": int(multiplicity),
            "sites": v.get("sites") or [],
            "site_counts": v.get("site_counts") or {},
            "site_label": v.get("label"),
            "base_smiles": base_smiles,
            "base_formula": base_formula,
            "note": v.get("label") or "",
        })
    return out


def molecules_for_f_counts(base: dict, counts: list[int],
                            max_variants: int = 20) -> tuple[list[dict], list[str]]:
    """Build the F=0,1,2,… series for one AI-generated base molecule.

    ``counts`` may contain 0: it yields the untouched parent (F0) so a
    ``0-3`` series has four complete structures.  Returns ``(molecules, notes)``.
    """
    name = base.get("name") or "mol"
    smiles = base.get("smiles") or ""
    formula = base.get("formula") or ""
    charge = int(base.get("charge") or 0)
    multiplicity = int(base.get("multiplicity") or 1)
    out: list[dict] = []
    notes: list[str] = []
    for n in counts:
        if int(n) <= 0:
            out.append({
                "name": f"{name}-F0-0",
                "smiles": smiles,
                "formula": formula,
                "n_F": 0,
                "natoms": base.get("natoms"),
                "charge": charge,
                "multiplicity": multiplicity,
                "sites": [],
                "site_counts": {},
                "site_label": "未取代基准（F0）",
                "base_smiles": smiles,
                "base_formula": formula,
                "note": "未取代基准（F0）",
            })
            continue
        try:
            variants = enumerate_fluorinated_variants(
                smiles, n_f=int(n), max_variants=max_variants)
        except ValueError as exc:
            notes.append(f"{name}: {exc}")
            continue
        out.extend(molecules_from_variants(
            name, smiles, variants, charge=charge, multiplicity=multiplicity))
    return out, notes
