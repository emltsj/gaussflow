"""System-level analysis: Li⁺ binding energies (with BSSE), redox potentials,
frontier orbitals and thermochemistry, aggregated over a calculation batch.

Energy bookkeeping (all in Hartree unless stated):

* ``optfreq`` level E  : SCF energy of the opt+freq job in the recipe
* ``sp_gas``  level E  : gas-phase single point in the recipe
* ``sp_sol``  level E  : solvent-phase single point in the recipe
* G from a level      : E(level SP) + [G(optfreq) - E(optfreq)]
  (thermal correction from the frequency job, applied to the higher-level
  electronic energy -- the usual rigid-rotor/harmonic-oscillator protocol)

Binding energy::

    Eb = E(complex) - E(monomer) - E(Li+)          (+ BSSE correction)

With Gaussian's counterpoise convention ``corrected = raw + BSSE`` (BSSE > 0),
so ``Eb(BSSE) = Eb(raw) + BSSE``.

Redox::
    dG_ox(sol) = G_sol(cation) - G_sol(neutral)     [eV after conversion]
    E_ox       = -dG_ox(sol) + reference_offset    [V]
"""
from __future__ import annotations

import json
from pathlib import Path

from .parse import HARTREE2EV, HARTREE2KCAL
from .pipeline import study_summary

F_CONST = 96485.33212  # C/mol (only needed if absolute potentials are wanted)

LEVELS = ["optfreq", "sp_gas", "sp_sol"]
LEVEL_LABELS = {
    "optfreq": "opt+freq",
    "sp_gas": "gas-phase SP",
    "sp_sol": "solvent-phase SP",
}
LI_STAGE = {"optfreq": "li_opt", "sp_gas": "li_gas", "sp_sol": "li_sol"}
BSSE_STAGE = {"optfreq": "bsse_optfreq", "sp_gas": "bsse_sp_gas",
              "sp_sol": "bsse_sp_sol"}


# ---------------------------------------------------------------------------
def _stage(summary: dict, name: str) -> dict:
    return (summary.get("stage_results") or {}).get(name) or {}


def _energy(summary: dict, name: str) -> float | None:
    return _stage(summary, name).get("energy_hartree")


def _G_of_level(summary: dict, level: str) -> float | None:
    """G at a level: E(level) + thermal correction from the opt+freq job."""
    opt = _stage(summary, "optfreq")
    e_opt = opt.get("energy_hartree")
    g_opt = opt.get("freeenergy_hartree")
    e = _energy(summary, level)
    if e is None:
        return None
    if e_opt is not None and g_opt is not None:
        return e + (g_opt - e_opt)
    return e


def _bsse(summary: dict, level: str) -> float | None:
    cp = _stage(summary, BSSE_STAGE[level]).get("counterpoise") or {}
    return cp.get("bsse_energy_hartree")


def _counterpoise(summary: dict, level: str) -> dict:
    return _stage(summary, BSSE_STAGE[level]).get("counterpoise") or {}


def _role(summary: dict) -> str:
    return str(summary.get("role") or "")


def _is_li_complex(role: str) -> bool:
    return (role.startswith("li_") and role != "li_ref"
            and "redox" not in role)


def _is_redox(role: str) -> bool:
    return "redox" in role


# ---------------------------------------------------------------------------
def system_summary(store, cfg, batch_id: str) -> dict:
    """Aggregate one calculation system (batch of studies)."""
    studies = store.list_studies(batch_id=batch_id, limit=500)
    if not studies:
        return {"ok": False, "error": f"未知 batch: {batch_id}"}
    summaries = {s["id"]: study_summary(store, cfg, s["id"]) for s in studies}
    neutral = next((v for v in summaries.values() if _role(v) == "neutral"), None)
    li_ref = next((v for v in summaries.values() if _role(v) == "li_ref"), None)
    li_complexes = [v for v in summaries.values() if _is_li_complex(_role(v))]
    redox = [v for v in summaries.values() if _is_redox(_role(v))]
    base_name = (studies[0].get("name") or "").split("-")[0]

    # ---- overview -----------------------------------------------------
    overview = []
    for sid, s in summaries.items():
        sr = s.get("stage_results") or {}
        overview.append({
            "study_id": sid,
            "name": s.get("name"),
            "role": _role(s),
            "formula": s.get("formula"),
            "natoms": s.get("natoms"),
            "charge": s.get("charge"),
            "multiplicity": s.get("multiplicity"),
            "status": s.get("status"),
            "validation": s.get("validation_overall"),
            "E_optfreq": _energy(s, "optfreq"),
            "E_sp_gas": _energy(s, "sp_gas"),
            "E_sp_sol": _energy(s, "sp_sol"),
            "G_optfreq": (sr.get("optfreq") or {}).get("freeenergy_hartree"),
            "homo_ev": s.get("orbitals", {}).get("homo_ev"),
            "lumo_ev": s.get("orbitals", {}).get("lumo_ev"),
            "gap_ev": s.get("orbitals", {}).get("gap_ev"),
            "n_imag": s.get("n_imag"),
            "walltime_s": s.get("walltime_total_s"),
        })

    # ---- Li+ binding energies -----------------------------------------
    binding: list[dict] = []
    if neutral is not None and li_complexes:
        mono = {lv: _energy(neutral, lv) for lv in LEVELS}
        li_e = {lv: (_energy(li_ref, LI_STAGE[lv]) if li_ref else None)
                for lv in LEVELS}
        for comp in li_complexes:
            site = comp.get("li_site") or {}
            for lv in LEVELS:
                ec = _energy(comp, lv)
                em, el = mono.get(lv), li_e.get(lv)
                ebind = (ec - em - el) if None not in (ec, em, el) else None
                bsse = _bsse(comp, lv)
                ebind_bsse = (ebind + bsse) if (ebind is not None and bsse is not None) else None
                cp = _counterpoise(comp, lv)
                binding.append({
                    "complex": comp.get("name"),
                    "complex_role": _role(comp),
                    "site": f"{site.get('element') or ''}{site.get('site_no') or ''}",
                    "level": lv,
                    "level_label": LEVEL_LABELS[lv],
                    "E_complex_hartree": ec,
                    "E_monomer_hartree": em,
                    "E_li_hartree": el,
                    "E_binding_hartree": ebind,
                    "E_binding_kcal": ebind * HARTREE2KCAL if ebind is not None else None,
                    "E_binding_ev": ebind * HARTREE2EV if ebind is not None else None,
                    "BSSE_hartree": bsse,
                    "BSSE_kcal": bsse * HARTREE2KCAL if bsse is not None else None,
                    "E_binding_BSSE_hartree": ebind_bsse,
                    "E_binding_BSSE_kcal": (ebind_bsse * HARTREE2KCAL
                                            if ebind_bsse is not None else None),
                    "E_binding_BSSE_ev": (ebind_bsse * HARTREE2EV
                                          if ebind_bsse is not None else None),
                    "vertical_CP_corrected_kcal": cp.get("complexation_corrected_kcal"),
                    "vertical_CP_raw_kcal": cp.get("complexation_raw_kcal"),
                    "counterpoise_ok": bool(cp),
                })

    # ---- redox ---------------------------------------------------------
    redox_rows: list[dict] = []
    offset = float(cfg.get("redox.reference_offset_v") or 0.0)
    for state_sum in redox:
        role = _role(state_sum)
        state = "cation" if role.endswith("cation") else (
            "anion" if role.endswith("anion") else "unknown")
        tag = role.split("_redox_")[0] if "_redox_" in role else ""
        if tag and tag != "redox":
            parent = next((v for v in summaries.values()
                           if _role(v) == "li_" + tag[2:]), None)
        else:
            parent = neutral
        if parent is None:
            continue
        gg_state = _G_of_level(state_sum, "sp_gas")
        gs_state = _G_of_level(state_sum, "sp_sol")
        gg_parent = _G_of_level(parent, "sp_gas")
        gs_parent = _G_of_level(parent, "sp_sol")
        dg_gas = (gg_state - gg_parent) if None not in (gg_state, gg_parent) else None
        dg_sol = (gs_state - gs_parent) if None not in (gs_state, gs_parent) else None
        dg_sol_ev = dg_sol * HARTREE2EV if dg_sol is not None else None
        dg_gas_ev = dg_gas * HARTREE2EV if dg_gas is not None else None
        potential = None
        if dg_sol_ev is not None:
            potential = -dg_sol_ev + offset
        redox_rows.append({
            "species": state_sum.get("name"),
            "parent": parent.get("name"),
            "state": state,
            "level_pair": "sp_sol (PCM/DMSO) vs sp_gas",
            "G_gas_parent_hartree": gg_parent,
            "G_gas_state_hartree": gg_state,
            "G_sol_parent_hartree": gs_parent,
            "G_sol_state_hartree": gs_state,
            "dG_gas_hartree": dg_gas,
            "dG_gas_ev": dg_gas_ev,
            "dG_sol_hartree": dg_sol,
            "dG_sol_ev": dg_sol_ev,
            "dG_sol_kcal": dg_sol * HARTREE2KCAL if dg_sol is not None else None,
            "E_potential_v": potential,
            "reference_offset_v": offset,
            "note": ("E = -ΔG_sol(eV) + offset  (相对参考；offset 需按文献参考电极设置)"
                     if state in ("cation", "anion") else ""),
        })

    # ---- frontier orbitals / thermochemistry ---------------------------
    orbitals = []
    thermo = []
    for s in summaries.values():
        o = s.get("orbitals") or {}
        if o.get("homo_ev") is not None:
            orbitals.append({
                "study": s.get("name"), "role": _role(s),
                "HOMO_eV": o.get("homo_ev"), "LUMO_eV": o.get("lumo_ev"),
                "gap_eV": o.get("gap_ev"),
                "HOMO_hartree": o.get("homo_hartree"),
                "LUMO_hartree": o.get("lumo_hartree"),
            })
        sr = s.get("stage_results") or {}
        opt = sr.get("optfreq") or {}
        shermo = {}
        if opt.get("job_id"):
            res = store.get_result(opt["job_id"]) or {}
            shermo = ((res.get("parsed") or {}).get("shermo") or {})
        thermo.append({
            "study": s.get("name"), "role": _role(s),
            "ZPE_hartree": opt.get("zpe_hartree"),
            "H_298_hartree": opt.get("enthalpy_hartree"),
            "G_298_hartree": opt.get("freeenergy_hartree"),
            "S_cal_mol_K": opt.get("entropy_cal_mol_k"),
            "Shermo_G_hartree": shermo.get("G_hartree") if shermo.get("ok") else None,
            "Shermo_available": bool(shermo.get("ok")),
            "Shermo_error": None if shermo.get("ok") else shermo.get("error"),
        })

    # ---- warnings / notes ----------------------------------------------
    warnings: list[str] = []
    notes: list[str] = []
    for s in summaries.values():
        if s.get("status") == "failed":
            diag = s.get("last_diagnosis") or {}
            warnings.append(f"{s.get('name')}: 失败 ({diag.get('class')})")
        elif s.get("validation_overall") == "warn":
            warnings.append(f"{s.get('name')}: 验收有警告")
        for note in (store.get_study(s["study_id"]) or {}).get("meta", {}).get("notes") or []:
            if note not in notes:
                notes.append(note)
    return {
        "ok": True,
        "batch_id": batch_id,
        "system_name": base_name,
        "recipe": next((s.get("recipe") for s in summaries.values() if s.get("recipe")), None),
        "solvent": next((s.get("solvent") for s in studies if s.get("solvent")), None),
        "n_studies": len(studies),
        "studies": overview,
        "binding": binding,
        "redox": redox_rows,
        "orbitals": orbitals,
        "thermo": thermo,
        "warnings": warnings,
        "notes": notes,
        "reference_offset_v": offset,
    }


# ---------------------------------------------------------------------------
def export_system(store, cfg, batch_id: str, out_dir: str | Path | None = None) -> dict:
    """Write batch artifacts + an analysis workbook (binding/redox/orbitals)."""
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    from . import export as base_export

    analysis = system_summary(store, cfg, batch_id)
    if not analysis.get("ok"):
        return analysis
    if out_dir is None:
        out_dir = cfg.exports_dir / f"system_{batch_id}"
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # standard batch artifacts (results.csv / summary.xlsx / xyz / report)
    base = base_export.export_batch(store, cfg, batch_id, out_dir=str(out / "batch"))

    (out / "analysis.json").write_text(
        json.dumps(analysis, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8")

    def write_csv(name: str, rows: list[dict]) -> str | None:
        if not rows:
            return None
        import csv
        path = out / name
        keys = list(rows[0].keys())
        with open(path, "w", encoding="utf-8-sig", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=keys)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: ("" if row.get(k) is None else row.get(k))
                                 for k in keys})
        return str(path)

    csvs = {
        "binding_csv": write_csv("binding_energies.csv", analysis["binding"]),
        "redox_csv": write_csv("redox_potentials.csv", analysis["redox"]),
        "orbitals_csv": write_csv("frontier_orbitals.csv", analysis["orbitals"]),
    }

    wb = Workbook()
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="1F4E78")

    def sheet(ws, rows, columns):
        ws.append([label for _, label in columns])
        for cell in ws[1]:
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center")
        for row in rows:
            ws.append(["" if row.get(k) is None else row.get(k) for k, _ in columns])
        ws.freeze_panes = "A2"
        for idx, (_, label) in enumerate(columns, start=1):
            ws.column_dimensions[get_column_letter(idx)].width = max(
                12, min(38, len(label) + 6))
        return ws

    ws = wb.active
    ws.title = "体系总览"
    sheet(ws, analysis["studies"], [
        ("name", "研究"), ("role", "角色"), ("formula", "分子式"),
        ("charge", "电荷"), ("multiplicity", "多重度"), ("status", "状态"),
        ("validation", "验收"), ("E_optfreq", "E(opt) / Ha"),
        ("E_sp_gas", "E(gas SP) / Ha"), ("E_sp_sol", "E(solv SP) / Ha"),
        ("homo_ev", "HOMO / eV"), ("lumo_ev", "LUMO / eV"),
        ("gap_ev", "Gap / eV"), ("n_imag", "NImag"),
    ])
    for row in ws.iter_rows(min_row=2, min_col=8, max_col=10):
        for cell in row:
            cell.number_format = "0.0000000000"

    sheet(wb.create_sheet("Li结合能"), analysis["binding"], [
        ("complex", "配合物"), ("site", "位点"), ("level_label", "理论级别"),
        ("E_complex_hartree", "E(complex)/Ha"), ("E_monomer_hartree", "E(monomer)/Ha"),
        ("E_li_hartree", "E(Li+)/Ha"),
        ("E_binding_kcal", "Eb(未校正)/(kcal/mol)"),
        ("BSSE_kcal", "BSSE/(kcal/mol)"),
        ("E_binding_BSSE_kcal", "Eb(BSSE校正)/(kcal/mol)"),
        ("E_binding_BSSE_ev", "Eb(BSSE校正)/eV"),
        ("vertical_CP_corrected_kcal", "垂直CP校正/(kcal/mol)"),
    ])
    sheet(wb.create_sheet("氧化还原"), analysis["redox"], [
        ("species", "物种"), ("parent", "母体"), ("state", "状态"),
        ("G_sol_state_hartree", "G(sol) 物种/Ha"),
        ("G_sol_parent_hartree", "G(sol) 母体/Ha"),
        ("dG_sol_ev", "ΔG_sol / eV"),
        ("dG_sol_kcal", "ΔG_sol / kcal/mol"),
        ("dG_gas_ev", "ΔG_gas / eV"),
        ("E_potential_v", "E / V (相对)"),
        ("note", "说明"),
    ])
    sheet(wb.create_sheet("前线轨道"), analysis["orbitals"], [
        ("study", "研究"), ("role", "角色"), ("HOMO_eV", "HOMO / eV"),
        ("LUMO_eV", "LUMO / eV"), ("gap_eV", "Gap / eV"),
    ])
    sheet(wb.create_sheet("热化学"), analysis["thermo"], [
        ("study", "研究"), ("role", "角色"), ("ZPE_hartree", "ZPE / Ha"),
        ("H_298_hartree", "H(298) / Ha"), ("G_298_hartree", "G(298) / Ha"),
        ("S_cal_mol_K", "S / cal/mol/K"),
        ("Shermo_G_hartree", "Shermo G / Ha"),
        ("Shermo_available", "Shermo 可用"),
    ])
    xlsx = out / "analysis.xlsx"
    wb.save(xlsx)

    report = out / "analysis_report.md"
    report.write_text(_render_report(analysis), encoding="utf-8")
    return {"ok": True, "batch_id": batch_id, "out_dir": str(out),
            "analysis_xlsx": str(xlsx), "analysis_json": str(out / "analysis.json"),
            "report_md": str(report), **csvs,
            "batch_artifacts": base}


def _render_report(a: dict) -> str:
    lines = [f"# 体系分析报告 — {a.get('system_name')}", "",
             f"- Batch: `{a['batch_id']}`",
             f"- 计算方案: {a.get('recipe')}",
             f"- 溶剂: {a.get('solvent')}",
             f"- 研究数: {a['n_studies']}",
             f"- 参考电位平移: {a.get('reference_offset_v')} V", ""]
    if a["binding"]:
        lines += ["## Li⁺ 结合能", "",
                  "| 配合物 | 位点 | 理论级别 | Eb(BSSE校正) kcal/mol | BSSE kcal/mol |",
                  "|---|---|---|---|---|"]
        for row in a["binding"]:
            if row.get("E_binding_BSSE_kcal") is None:
                continue
            lines.append(
                f"| {row['complex']} | {row['site']} | {row['level_label']} "
                f"| {row['E_binding_BSSE_kcal']:.2f} | "
                f"{(row.get('BSSE_kcal') or 0):.2f} |")
        lines += ["", "> 结合能定义：Eb = E(complex) − E(monomer) − E(Li⁺) + BSSE；"
                      "BSSE 由 counterpoise=2 单点给出（Gaussian 约定 corrected = raw + BSSE）。", ""]
    if a["redox"]:
        lines += ["## 氧化还原（相对）", "",
                  "| 物种 | 状态 | ΔG_sol eV | E / V (相对) |", "|---|---|---|---|"]
        for row in a["redox"]:
            lines.append(f"| {row['species']} | {row['state']} "
                         f"| {row.get('dG_sol_ev') or 0:.3f} "
                         f"| {row.get('E_potential_v') or 0:.3f} |")
        lines += ["", "> E = -ΔG_sol(eV) + offset；绝对电位需按文献参考电极设置 offset。", ""]
    if a["orbitals"]:
        lines += ["## 前线轨道", "", "| 研究 | HOMO / eV | LUMO / eV | Gap / eV |",
                  "|---|---|---|---|"]
        for row in a["orbitals"]:
            lines.append(f"| {row['study']} | {row.get('HOMO_eV') or 0:.3f} "
                         f"| {row.get('LUMO_eV') or 0:.3f} | {row.get('gap_eV') or 0:.3f} |")
        lines.append("")
    if a["warnings"]:
        lines += ["## 警告", ""] + [f"- {w}" for w in a["warnings"]] + [""]
    return "\n".join(lines) + "\n"
