"""Deterministic parser for Gaussian 16 output (.out).

Primary extraction uses cclib; a regex fallback covers truncated / unusual
outputs.  The parser never *decides* whether a result is scientifically
acceptable -- that is :mod:`gaussflow.validate`.  It only extracts facts and
normalises units into JSON-safe Python types.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

HARTREE2KCAL = 627.5094740631
HARTREE2EV = 27.211386245988
HARTREE2KJ = 2625.4996394799
HARTREE2CAL = 627509.4740631          # for cclib entropy (hartree/(mol K))

_SCI = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"


def _f(x) -> float:
    return float(str(x).replace("D", "E").replace("d", "e"))


def _last(seq):
    return seq[-1] if seq else None


def _enum(seq):
    return list(seq) if seq else []


# ---------------------------------------------------------------------------
# geometry fallback
# ---------------------------------------------------------------------------
def _parse_standard_orientations(text: str) -> list[list[tuple[str, float, float, float]]]:
    from .structure import element_symbol
    blocks = []
    pattern = re.compile(
        r"Standard orientation:.*?-{10,}\s*\n.*?-{10,}\s*\n(.*?)-{10,}",
        re.S)
    for m in pattern.finditer(text):
        rows = []
        for line in m.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 6:
                try:
                    z = int(parts[1])
                    rows.append((element_symbol(z), _f(parts[3]), _f(parts[4]), _f(parts[5])))
                except (ValueError, IndexError):
                    continue
        if rows:
            blocks.append(rows)
    return blocks


def _parse_input_orientation(text: str) -> list[tuple[str, float, float, float]]:
    from .structure import element_symbol
    m = None
    for m2 in re.finditer(
            r"Input orientation:.*?-{10,}\s*\n.*?-{10,}\s*\n(.*?)-{10,}", text, re.S):
        m = m2
    if not m:
        return []
    rows = []
    for line in m.group(1).strip().splitlines():
        parts = line.split()
        if len(parts) >= 6:
            try:
                rows.append((element_symbol(int(parts[1])), _f(parts[3]), _f(parts[4]), _f(parts[5])))
            except (ValueError, IndexError):
                continue
    return rows


# ---------------------------------------------------------------------------
# main entry points
# ---------------------------------------------------------------------------
def parse_output(path: str | Path, cclib_attempt: bool = True) -> dict:
    path = Path(path)
    if not path.exists():
        return {"output_file": str(path), "exists": False,
                "normal_termination": False, "error_termination": True,
                "error_lines": ["输出文件不存在"], "warnings": [],
                "size_bytes": 0}
    text = path.read_text(encoding="utf-8", errors="replace")
    return parse_text(text, path=str(path), cclib_attempt=cclib_attempt)


def parse_text(text: str, path: str | None = None, cclib_attempt: bool = False) -> dict:
    r: dict = {
        "output_file": path,
        "exists": True,
        "size_bytes": len(text.encode("utf-8", errors="replace")),
        "cclib_ok": False,
        "cclib_error": None,
    }
    norm = list(re.finditer(r"Normal termination of Gaussian \d+ at (.+)", text))
    err_lines = [ln.strip() for ln in text.splitlines()
                 if re.search(r"Error termination|Convergence failure|"
                              r"Erroneous write|End of file in|QPErr|"
                              r"Illegal IType|Syntax error", ln)][:12]
    r["normal_termination"] = bool(norm)
    r["termination_time"] = norm[-1].group(1).strip() if norm else None
    r["error_termination"] = bool(err_lines) or not r["normal_termination"]
    r["error_lines"] = err_lines

    m = re.search(r"Gaussian 16, Revision ([A-Za-z0-9.\-]+)", text)
    r["version"] = f"Gaussian 16 Rev {m.group(1)}" if m else None

    # route line
    m = re.search(r"^\s*#\s*[pPnNtT]?\s*(.+)$", text, re.M)
    r["route"] = m.group(1).strip() if m else None

    # atom count / charge / multiplicity
    m = _last(re.findall(r"NAtoms=\s*(\d+)", text))
    r["natoms"] = int(m) if m else None
    chargemult = re.findall(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text)
    if chargemult:
        r["charge"], r["multiplicity"] = int(chargemult[-1][0]), int(chargemult[-1][1])
    else:
        r["charge"], r["multiplicity"] = None, None

    # SCF energies
    scf = re.findall(
        r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(" + _SCI + r")\s+A\.U\.\s+after\s+(\d+) cycles",
        text)
    r["scf_energies"] = [{"energy_hartree": _f(e), "cycles": int(c)} for _, e, c in scf]
    r["final_energy_hartree"] = _f(scf[-1][1]) if scf else None
    r["scf_cycles"] = int(scf[-1][2]) if scf else None
    r["scf_converged"] = bool(scf) and ("Convergence failure" not in text)

    # geometry optimisation
    is_opt = bool(re.search(r"Optimization completed|Stationary point found|"
                            r"Step number\s+\d+|Optimization stopped", text))
    r["is_opt"] = is_opt
    if is_opt:
        stationary_found = ("-- Stationary point found." in text)
        stationary_not = ("Stationary point NOT found" in text)
        completed = ("Optimization completed." in text)
        stopped = ("Optimization stopped." in text)
        steps = re.findall(r"Step number\s+(\d+)", text)
        criteria = {}
        for item in ("Maximum Force", "RMS     Force", "RMS Force",
                     "Maximum Displacement", "RMS     Displacement", "RMS Displacement"):
            matches = re.findall(
                rf"{re.escape(item)}\s+({_SCI})\s+({_SCI})\s+(YES|NO)", text)
            if matches:
                v, thr, ok = matches[-1]
                key = item.replace("     ", " ")
                criteria[key] = {"value": _f(v), "threshold": _f(thr),
                                 "converged": ok == "YES"}
        r["opt"] = {
            "stationary_point": stationary_found and not stationary_not,
            "completed": completed,
            "stopped": stopped,
            "nsteps": int(steps[-1]) if steps else None,
            "criteria": criteria,
            "criteria_all_yes": (all(c["converged"] for c in criteria.values())
                                 if criteria else None),
            # Gaussian's own verdict is authoritative: it may report
            # "Optimization completed"/"Stationary point found" while one of the
            # printed displacement criteria is a hair above its threshold.
            "converged": bool(completed and stationary_found and not stationary_not
                              and not stopped),
        }
    else:
        r["opt"] = None

    # frequencies + thermochemistry
    freqs: list[float] = []
    vib_source = "regex"
    zpve = enthalpy = freeenergy = entropy = None
    dipole = None
    last_geom: list[tuple[str, float, float, float]] = []

    cclib_data = None
    if cclib_attempt and path:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                from cclib.io import ccread
                cclib_data = ccread(str(path))
            r["cclib_ok"] = cclib_data is not None
        except Exception as exc:  # pragma: no cover - depends on output
            r["cclib_error"] = f"{type(exc).__name__}: {exc}"
            cclib_data = None

    if cclib_data is not None:
        vib = getattr(cclib_data, "vibfreqs", None)
        if vib is not None:
            freqs = [float(x) for x in vib]
            vib_source = "cclib"
        for attr in ("zpve", "enthalpy", "freeenergy", "entropy"):
            val = getattr(cclib_data, attr, None)
            if val is not None:
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    val = None
            if attr == "zpve":
                zpve = val
            elif attr == "enthalpy":
                enthalpy = val
            elif attr == "freeenergy":
                freeenergy = val
            else:
                # cclib reports entropy in hartree/(mol K); normalise to cal/(mol K)
                if val is not None and abs(val) < 1.0:
                    val = val * HARTREE2CAL
                entropy = val
        diffs = getattr(cclib_data, "moments", None)
        if diffs is not None and len(diffs) > 0:
            try:
                m0 = [float(x) for x in diffs[0][:3]]
                if any(abs(x) > 1e-8 for x in m0):
                    dipole = sum(x * x for x in m0) ** 0.5
            except (TypeError, ValueError):
                dipole = None
        coords = getattr(cclib_data, "atomcoords", None)
        atomnos = getattr(cclib_data, "atomnos", None)
        if coords is not None and atomnos is not None and len(coords):
            from .structure import element_symbol
            last_geom = [(element_symbol(int(z)), float(p[0]), float(p[1]), float(p[2]))
                         for z, p in zip(atomnos, coords[-1])]

    if not freqs:
        for m in re.finditer(r"Frequencies --\s+(.+)", text):
            vals = re.findall(_SCI, m.group(1))
            freqs.extend(_f(v) for v in vals)
    if not last_geom:
        blocks = _parse_standard_orientations(text)
        if blocks:
            last_geom = blocks[-1]
        else:
            last_geom = _parse_input_orientation(text)

    # regex thermochemistry fallback
    def grab(pattern: str):
        m = re.search(pattern, text)
        return _f(m.group(1)) if m and m.group(1) else None

    if zpve is None:
        zpve = grab(r"Sum of electronic and zero-point Energies=\s+(" + _SCI + ")")
    if enthalpy is None:
        enthalpy = grab(r"Sum of electronic and thermal Enthalpies=\s+(" + _SCI + ")")
    if freeenergy is None:
        freeenergy = grab(r"Sum of electronic and thermal Free Energies=\s+(" + _SCI + ")")
    zpc = grab(r"Zero-point correction=\s+(" + _SCI + ")")
    thermal_corr = {
        "zero_point_hartree": zpc,
        "thermal_energy": grab(r"Thermal correction to Energy=\s+(" + _SCI + ")"),
        "thermal_enthalpy": grab(r"Thermal correction to Enthalpy=\s+(" + _SCI + ")"),
        "thermal_gibbs": grab(r"Thermal correction to Gibbs Free Energy=\s+(" + _SCI + ")"),
    }
    if entropy is None:
        m = re.search(r"Total\s+(" + _SCI + r")\s+(" + _SCI + r")\s+(" + _SCI + r")\s*$",
                      text, re.M)
        if m:
            entropy = _f(m.group(3))

    n_imag = sum(1 for x in freqs if x < 0)
    r["freq"] = {
        "has_freq": bool(freqs),
        "n_freqs": len(freqs),
        "n_imag": n_imag,
        "imag_freqs": [x for x in freqs if x < 0],
        "frequencies_cm1": freqs,
        "source": vib_source,
        "zpve_hartree": zpve,
        "enthalpy_hartree": enthalpy,
        "freeenergy_hartree": freeenergy,
        "entropy_cal_mol_k": entropy,
        "thermal_corrections": thermal_corr,
        "lowest_modes": sorted(freqs)[:6],
    }

    # dipole: prefer the final archive line "Dipole=x,y,z" (always present),
    # then cclib, then the explicit "Dipole moment (...)" block
    if dipole is None or dipole == 0.0:
        dip_lines = re.findall(
            r"Dipole=\s*(" + _SCI + r"),\s*(" + _SCI + r"),\s*(" + _SCI + r")", text)
        if dip_lines:
            x, y, z = (_f(v) for v in dip_lines[-1])
            dipole = (x * x + y * y + z * z) ** 0.5
    if not dipole:
        dm = re.search(r"Dipole moment \(field-independent basis, Debye\):\s*\n", text)
        if dm:
            seg = text[dm.end():dm.end() + 400]
            x = re.search(r"X=\s*(" + _SCI + ")", seg)
            y = re.search(r"Y=\s*(" + _SCI + ")", seg)
            z = re.search(r"Z=\s*(" + _SCI + ")", seg)
            if x and y and z:
                dipole = (_f(x.group(1)) ** 2 + _f(y.group(1)) ** 2
                          + _f(z.group(1)) ** 2) ** 0.5
    r["dipole_debye"] = dipole

    r["final_geometry"] = [{"symbol": s, "x": x, "y": y, "z": z}
                           for s, x, y, z in last_geom]

    # timings
    det = re.search(r"Job cpu time:\s*(.+)", text)
    r["cpu_time"] = det.group(1).strip() if det else None
    det = re.search(r"Elapsed time:\s*(.+)", text)
    r["elapsed_time"] = det.group(1).strip() if det else None

    # warnings (Gaussian's legal boilerplate is not a scientific warning)
    warn_lines = []
    for ln in text.splitlines():
        if re.search(r"Warning|WARNING|Small basis|imaginary frequencies|"
                     r"Basis set superposition", ln):
            s = ln.strip()
            if "This program may not be used" in s:
                continue
            if s and s not in warn_lines and len(s) < 300:
                warn_lines.append(s)
    r["warnings"] = warn_lines[:20]

    # --- frontier orbitals (HOMO/LUMO) ---------------------------------
    orbitals = _parse_orbitals(text)
    r["orbitals"] = orbitals
    if orbitals.get("homo_ev") is not None and orbitals.get("lumo_ev") is not None:
        orbitals["gap_ev"] = orbitals["lumo_ev"] - orbitals["homo_ev"]

    # --- counterpoise / BSSE (for Li+ binding energies) -----------------
    r["counterpoise"] = _parse_counterpoise(text)

    # unit conversions (deterministic, for export layer)
    e = r["final_energy_hartree"]
    r["final_energy_kcal_mol"] = e * HARTREE2KCAL if e is not None else None
    r["final_energy_ev"] = e * HARTREE2EV if e is not None else None
    return r


def _nums(text: str) -> list[float]:
    return [_f(x) for x in re.findall(_SCI, text or "")]


def _parse_orbitals(text: str) -> dict:
    """Extract HOMO/LUMO (Hartree + eV) from the MO eigenvalue list."""
    occ: list[float] = []
    virt: list[float] = []
    current = None
    for line in text.splitlines():
        s = line.strip()
        if re.search(r"(?i)alpha\s+occ\.\s*eigenvalues", s):
            current = occ
            occ.extend(_nums(s.split("--")[-1] if "--" in s else ""))
            continue
        if re.search(r"(?i)alpha\s+virt\.\s*eigenvalues", s):
            current = virt
            virt.extend(_nums(s.split("--")[-1] if "--" in s else ""))
            continue
        if current is not None:
            if not s:
                current = None
                continue
            if re.fullmatch(r"[-+0-9.EeDd\s]+", s):
                current.extend(_nums(s))
            else:
                current = None
    out = {"homo_hartree": None, "homo_ev": None,
           "lumo_hartree": None, "lumo_ev": None,
           "n_occupied_basis": len(occ), "n_virtual_basis": len(virt)}
    if occ:
        out["homo_hartree"] = occ[-1]
        out["homo_ev"] = occ[-1] * HARTREE2EV
    if virt:
        out["lumo_hartree"] = virt[0]
        out["lumo_ev"] = virt[0] * HARTREE2EV
    return out


def _parse_counterpoise(text: str) -> dict | None:
    """Parse Gaussian counterpoise=2 output (BSSE, fragment energies).

    Verified output layout (Gaussian 16W):
      Counterpoise: doing full-system calculation.
      SCF Done: ...                                  <- complex energy
      Counterpoise: doing calculation for fragment 1 using the basis set of the
        full-system ... SCF Done: ...                <- fragment in complex basis
      Counterpoise: doing calculation for fragment 1 using the basis set of the
        fragment. ... SCF Done: ...                  <- fragment own basis
      ...
      Counterpoise corrected energy = ...
                        BSSE energy = ...
                   sum of fragments = ...
                complexation energy =  -6.42 kcal/mole (raw)
                complexation energy =  -2.82 kcal/mole (corrected)
    Note: BSSE is positive and *added*: corrected = raw + BSSE.
    """
    if "Counterpoise" not in text and "BSSE" not in text:
        return None
    out: dict = {
        "bsse_energy_hartree": None,
        "corrected_energy_hartree": None,
        "sum_fragments_hartree": None,
        "complex_energy_hartree": None,
        "complexation_raw_kcal": None,
        "complexation_corrected_kcal": None,
        "fragment_complex_basis": {},
        "fragment_own_basis": {},
        "raw_lines": [],
    }
    for m in re.finditer(r"Counterpoise corrected energy\s*=\s*(" + _SCI + ")", text):
        out["corrected_energy_hartree"] = _f(m.group(1))
    for m in re.finditer(r"BSSE energy\s*=\s*(" + _SCI + ")", text):
        out["bsse_energy_hartree"] = _f(m.group(1))
    for m in re.finditer(r"sum of fragments\s*=\s*(" + _SCI + ")", text):
        out["sum_fragments_hartree"] = _f(m.group(1))
    for m in re.finditer(
            r"complexation energy\s*=\s*(" + _SCI + r")\s*kcal/mole\s*\((raw|corrected)\)",
            text, re.I):
        key = (out["complexation_raw_kcal"] if m.group(2).lower() == "raw"
               else out["complexation_corrected_kcal"])
        out[f"complexation_{m.group(2).lower()}_kcal"] = _f(m.group(1))
    m = re.search(
        r"Counterpoise:\s*doing full-system calculation\.(.*?)"
        r"SCF Done:\s+E\([^)]*\)\s*=\s*(" + _SCI + ")", text, re.S)
    if m:
        out["complex_energy_hartree"] = _f(m.group(2))
    for m in re.finditer(
            r"doing calculation for fragment\s+(\d+)\s+using the basis set of the\s+"
            r"(full-system|fragment)\.?(.*?)SCF Done:\s+E\([^)]*\)\s*=\s*(" + _SCI + ")",
            text, re.S | re.I):
        key = ("fragment_complex_basis" if "full-system" in m.group(2).lower()
               else "fragment_own_basis")
        out[key][int(m.group(1))] = _f(m.group(4))
    out["raw_lines"] = [ln.strip() for ln in text.splitlines()
                        if any(k in ln for k in ("Counterpoise", "BSSE",
                                                 "sum of fragments",
                                                 "complexation energy"))][:24]
    return out


def output_tail(path: str | Path, n_chars: int = 4000) -> str:
    """Small tail excerpt used by diagnose_failure; never the whole log."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    return text[-n_chars:]
