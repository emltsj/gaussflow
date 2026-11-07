"""Deterministic aggregation and export.

The Agent never transcribes numbers: everything here is produced from the
parsed/validated results in the store.  Outputs per batch:

* ``results.csv``        one row per variant
* ``summary.xlsx``       Summary / Energies / Jobs / Frequencies / Failures / Metadata
* ``structures/*.xyz``   final optimised geometry of each variant
* ``batch_summary.json`` machine-readable aggregate
* ``report.md``          human-readable report with explicit validation caveats
"""
from __future__ import annotations

import csv
import json
import re
from datetime import datetime
from pathlib import Path

from .parse import HARTREE2KCAL
from .pipeline import study_summary
from .store import Store

SUMMARY_COLUMNS = [
    ("study_id", "Study ID"),
    ("name", "Name"),
    ("variant_label", "Variant"),
    ("formula", "Formula"),
    ("natoms", "NAtoms"),
    ("charge", "Charge"),
    ("multiplicity", "Mult"),
    ("route", "Route"),
    ("workflow", "Workflow"),
    ("status", "Status"),
    ("validation_overall", "Validation"),
    ("energy_hartree", "E (Hartree)"),
    ("energy_kcal_mol", "E (kcal/mol)"),
    ("zpe_hartree", "ZPE (Hartree)"),
    ("enthalpy_hartree", "H (Hartree)"),
    ("freeenergy_hartree", "G (Hartree)"),
    ("entropy_cal_mol_k", "S (cal/mol/K)"),
    ("n_imag", "NImag"),
    ("opt_converged", "OptConv"),
    ("opt_steps", "OptSteps"),
    ("scf_cycles", "SCFcycles"),
    ("walltime_total_s", "Walltime(s)"),
    ("final_stage", "FinalStage"),
    ("n_jobs", "Jobs"),
    ("last_error", "Error"),
]

FREQ_COLUMNS = [
    ("study_id", "Study ID"), ("name", "Name"), ("variant_label", "Variant"),
    ("n_imag", "NImag"), ("imag_freqs", "ImagFreqs"),
    ("lowest_modes", "Lowest modes (cm^-1)"),
]

JOB_COLUMNS = [
    ("study_id", "Study ID"), ("job_id", "Job ID"), ("stage", "Stage"),
    ("attempt", "Attempt"), ("status", "Status"), ("route", "Route"),
    ("walltime_s", "Walltime(s)"), ("validation", "Validation"),
    ("normal_termination", "NormalTerm"), ("energy_hartree", "E (Hartree)"),
    ("n_imag", "NImag"), ("opt_converged", "OptConv"),
    ("error_class", "ErrorClass"), ("output", "Output"),
]


def _safe(name: str) -> str:
    s = re.sub(r"[^0-9A-Za-z_.\-]+", "_", str(name or "study")).strip("_")
    return s[:60] or "study"


def _variant_label(s: dict) -> str:
    v = s.get("variant") or {}
    if v.get("type") == "ion":
        return v.get("ion") or "ion"
    return "base"


def _row(summary: dict) -> dict:
    s = dict(summary)
    s["variant_label"] = _variant_label(s)
    err = s.get("last_diagnosis") or {}
    s["last_error"] = f"{err.get('class')}: {err.get('title')}" if err else None
    return s


def _fmt(v, key: str | None = None):
    """Text formatting for CSV / Markdown (columns get their own precision)."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "YES" if v else "NO"
    if isinstance(v, float):
        if key in ("energy_hartree", "zpe_hartree", "enthalpy_hartree",
                   "freeenergy_hartree"):
            return f"{v:.10f}"
        if key in ("energy_kcal_mol", "relative_kcal_mol"):
            return f"{v:.2f}"
        if key == "entropy_cal_mol_k":
            return f"{v:.3f}"
        if key in ("walltime_s", "walltime_total_s"):
            return f"{v:.1f}"
        return f"{v:.6g}"
    if isinstance(v, (list, tuple)):
        return ", ".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, ensure_ascii=False)
    return v


def _excel_value(v):
    """Keep numbers numeric in Excel so charts/aggregations work."""
    if v is None:
        return ""
    if isinstance(v, bool):
        return "YES" if v else "NO"
    if isinstance(v, (int, float)):
        return v
    if isinstance(v, (list, tuple)):
        return ", ".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, ensure_ascii=False)
    return v


# ---------------------------------------------------------------------------
def export_study(store: Store, cfg, study_id: str) -> dict:
    summary = study_summary(store, cfg, study_id)
    if not summary.get("ok"):
        return summary
    run_dir = cfg.work_root / study_id
    run_dir.mkdir(parents=True, exist_ok=True)
    result_json = run_dir / "result.json"
    result_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                      default=str), encoding="utf-8")
    xyz_path = None
    geom = summary.get("final_geometry")
    if geom:
        xyz = f"{len(geom)}\n{summary['name']} {summary['study_id']}\n" + "\n".join(
            f"{g['symbol']:<2s} {g['x']:14.6f} {g['y']:14.6f} {g['z']:14.6f}"
            for g in geom) + "\n"
        xyz_path = run_dir / "final.xyz"
        xyz_path.write_text(xyz, encoding="utf-8")
    return {"ok": True, "study_id": study_id, "result_json": str(result_json),
            "xyz": str(xyz_path) if xyz_path else None, "summary": summary}


def export_batch(store: Store, cfg, batch_id: str, out_dir: str | None = None) -> dict:
    studies = store.list_studies(batch_id=batch_id)
    if not studies:
        return {"ok": False, "error": f"未知 batch: {batch_id}"}
    out = Path(out_dir) if out_dir else cfg.exports_dir / f"batch_{batch_id}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "structures").mkdir(exist_ok=True)

    summaries = [_row(study_summary(store, cfg, s["id"])) for s in studies]

    # --- CSV ----------------------------------------------------------
    csv_path = out / "results.csv"
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([label for _, label in SUMMARY_COLUMNS])
        for s in summaries:
            writer.writerow([_fmt(s.get(k), k) for k, _ in SUMMARY_COLUMNS])

    # --- XYZ ----------------------------------------------------------
    for s in summaries:
        geom = s.get("final_geometry")
        if not geom:
            continue
        xyz = f"{len(geom)}\n{s['name']} {s['study_id']}\n" + "\n".join(
            f"{g['symbol']:<2s} {g['x']:14.6f} {g['y']:14.6f} {g['z']:14.6f}"
            for g in geom) + "\n"
        (out / "structures" / f"{_safe(s['name'])}.xyz").write_text(xyz, encoding="utf-8")

    # --- relative energies -------------------------------------------
    # Only meaningful between variants with identical formula + charge
    # (e.g. the same molecule with different ions is NOT comparable).
    by_group: dict = {}
    for s in summaries:
        e = s.get("energy_hartree")
        s["relative_kcal_mol"] = None
        s["comparable_group"] = None
        if e is None:
            continue
        key = (s.get("formula"), s.get("charge"))
        s["comparable_group"] = f"{key[0]}/{key[1]:+d}" if key[1] is not None else str(key[0])
        by_group.setdefault(key, []).append(s)
    for key, grp in by_group.items():
        if len(grp) < 2:
            continue
        emin = min(x["energy_hartree"] for x in grp)
        for s in grp:
            s["relative_kcal_mol"] = (s["energy_hartree"] - emin) * HARTREE2KCAL

    # --- Excel --------------------------------------------------------
    xlsx_path = out / "summary.xlsx"
    _write_excel(xlsx_path, summaries, cfg, batch_id)

    # --- JSON + Markdown ---------------------------------------------
    batch_json = out / "batch_summary.json"
    batch_json.write_text(json.dumps(
        {"batch_id": batch_id, "exported_at": datetime.now().isoformat(timespec="seconds"),
         "n_studies": len(summaries), "studies": summaries},
        ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    report = out / "report.md"
    report.write_text(_render_report(batch_id, summaries), encoding="utf-8")

    return {"ok": True, "batch_id": batch_id, "out_dir": str(out),
            "results_csv": str(csv_path), "summary_xlsx": str(xlsx_path),
            "batch_summary_json": str(batch_json), "report_md": str(report),
            "n_studies": len(summaries)}


# ---------------------------------------------------------------------------
def _write_excel(path: Path, summaries: list[dict], cfg, batch_id: str) -> None:
    from openpyxl import Workbook
    from openpyxl.chart import BarChart, Reference
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="1F4E78")
    warn_fill = PatternFill("solid", fgColor="FFF2CC")
    fail_fill = PatternFill("solid", fgColor="F8CBAD")

    def fill_sheet(ws, columns, rows, freeze="A2"):
        ws.append([label for _, label in columns])
        for cell in ws[1]:
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center")
        for r in rows:
            ws.append([_excel_value(r.get(k)) for k, _ in columns])
        ws.freeze_panes = freeze
        ws.auto_filter.ref = ws.dimensions
        formats = {
            "E (Hartree)": "0.0000000000", "ZPE (Hartree)": "0.0000000000",
            "H (Hartree)": "0.0000000000", "G (Hartree)": "0.0000000000",
            "E (kcal/mol)": "0.00", "S (cal/mol/K)": "0.000",
            "Walltime(s)": "0.0", "Walltime (s)": "0.0",
        }
        for idx, (key, label) in enumerate(columns, start=1):
            width = max(10, min(34, len(label) + 4))
            ws.column_dimensions[get_column_letter(idx)].width = width
            fmt = formats.get(label)
            if fmt:
                for (cell,) in ws.iter_rows(min_row=2, min_col=idx, max_col=idx):
                    cell.number_format = fmt
        return ws

    # Summary
    ws = wb.active
    ws.title = "Summary"
    fill_sheet(ws, SUMMARY_COLUMNS, summaries)
    for row in ws.iter_rows(min_row=2):
        status = row[9].value
        val = row[10].value
        if val == "fail" or status == "failed":
            for c in row:
                c.fill = fail_fill
        elif val == "warn":
            for c in row:
                c.fill = warn_fill

    # Energies + bar chart
    ws_e = wb.create_sheet("Energies")
    ws_e.append(["Name", "Variant", "Group", "E (Hartree)",
                 "Relative (kcal/mol)"])
    for cell in ws_e[1]:
        cell.font = header_font
        cell.fill = header_fill
    for s in summaries:
        ws_e.append([s.get("name"), _variant_label(s),
                     s.get("comparable_group") or "",
                     s.get("energy_hartree"),
                     s.get("relative_kcal_mol")])
    ws_e.freeze_panes = "A2"
    for col, width in (("A", 28), ("B", 14), ("C", 18), ("D", 20), ("E", 20)):
        ws_e.column_dimensions[col].width = width
    for (cell,) in ws_e.iter_rows(min_row=2, min_col=4, max_col=4):
        cell.number_format = "0.0000000000"
    for (cell,) in ws_e.iter_rows(min_row=2, min_col=5, max_col=5):
        cell.number_format = "0.00"
    if summaries:
        chart = BarChart()
        chart.type = "col"
        chart.title = "Relative energy (kcal/mol, within same composition)"
        chart.y_axis.title = "kcal/mol"
        data = Reference(ws_e, min_col=5, min_row=1, max_row=len(summaries) + 1)
        cats = Reference(ws_e, min_col=1, min_row=2, max_row=len(summaries) + 1)
        chart.add_data(data, titles_from_data=True)
        chart.set_categories(cats)
        chart.height = 8
        chart.width = 18
        ws_e.add_chart(chart, "G2")

    # Frequencies
    ws_f = wb.create_sheet("Frequencies")
    fill_sheet(ws_f, FREQ_COLUMNS, summaries)
    for row in ws_f.iter_rows(min_row=2):
        try:
            n_imag = int(row[3].value or 0)
        except (TypeError, ValueError):
            n_imag = 0
        if n_imag > 0:
            for c in row:
                c.fill = fail_fill

    # Jobs
    ws_j = wb.create_sheet("Jobs")
    job_rows = []
    for s in summaries:
        for j in s.get("jobs") or []:
            job_rows.append({**j, "study_id": s["study_id"]})
    fill_sheet(ws_j, JOB_COLUMNS, job_rows)

    # Failures
    ws_x = wb.create_sheet("Failures")
    ws_x.append(["Study ID", "Name", "Status", "Error class", "Error title",
                 "Needs human", "Suggestions", "Validation reasons"])
    for cell in ws_x[1]:
        cell.font = header_font
        cell.fill = header_fill
    for s in summaries:
        err = s.get("last_diagnosis") or {}
        if s.get("status") != "failed" and s.get("validation_overall") != "fail":
            continue
        ws_x.append([s.get("study_id"), s.get("name"), s.get("status"),
                     err.get("class"), err.get("title"),
                     err.get("needs_human"), "; ".join(err.get("suggestions") or []),
                     "; ".join(s.get("validation_reasons") or [])])
    for col, width in (("A", 22), ("B", 24), ("C", 12), ("D", 18), ("E", 26),
                       ("F", 12), ("G", 40), ("H", 40)):
        ws_x.column_dimensions[col].width = width

    # Metadata
    ws_m = wb.create_sheet("Metadata")
    meta = {
        "Batch": batch_id,
        "Exported at": datetime.now().isoformat(timespec="seconds"),
        "G16 exe": str(cfg.g16_exe),
        "Default route": cfg.chemistry.get("route"),
        "Pre-opt route": cfg.chemistry.get("preopt_route"),
        "Resources": json.dumps(cfg.resources, ensure_ascii=False),
        "Validation rules": json.dumps(cfg.validation, ensure_ascii=False),
        "Note": "Energies are electronic energies in Hartree; relative values are in "
                "kcal/mol vs. the lowest variant in this batch. Frequency checks are "
                "only valid for the level of theory of the final stage.",
    }
    ws_m.append(["Key", "Value"])
    for cell in ws_m[1]:
        cell.font = header_font
        cell.fill = header_fill
    for k, v in meta.items():
        ws_m.append([k, v])
    ws_m.column_dimensions["A"].width = 22
    ws_m.column_dimensions["B"].width = 90

    wb.save(path)


def _render_report(batch_id: str, summaries: list[dict]) -> str:
    lines = [f"# Batch report `{batch_id}`", "",
             f"Exported: {datetime.now().isoformat(timespec='seconds')}", ""]
    ok = sum(1 for s in summaries if s.get("status") == "completed"
             and s.get("validation_overall") in ("pass", "warn"))
    lines.append(f"- Variants: **{len(summaries)}**, accepted: **{ok}**")
    groups: dict[str, list[str]] = {}
    for s in summaries:
        if s.get("comparable_group"):
            groups.setdefault(s["comparable_group"], []).append(s.get("name") or "")
    if groups:
        lines.append("- Comparable groups (relative energies only valid inside a group): "
                     + "; ".join(f"`{k}`: {', '.join(v)}" for k, v in groups.items()))
    lines.append("")
    lines.append("| Name | Variant | Formula | Charge | E (Hartree) | E_rel (kcal/mol) "
                 "| NImag | OptConv | Status |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for s in sorted(summaries, key=lambda x: (x.get("energy_hartree") is None,
                                              x.get("energy_hartree") or 0)):
        rel = s.get("relative_kcal_mol")
        lines.append(
            f"| {s.get('name')} | {_variant_label(s)} | {s.get('formula') or ''} "
            f"| {s.get('charge')} | {_fmt(s.get('energy_hartree'), 'energy_hartree')} "
            f"| {'' if rel is None else _fmt(rel, 'relative_kcal_mol')} "
            f"| {_fmt(s.get('n_imag'))} | {s.get('opt_converged')} | {s.get('status')} |")
    lines += ["", "## Validation notes", "",
              "- Only results with normal termination, converged SCF and converged "
              "optimisation are marked accepted.",
              "- For minimum structures the frequency job must show 0 imaginary modes.",
              "- Final electronic energy is taken from the last stage of the "
              "opt -> freq chain (same level of theory).",
              "- Relative energies are only computed between variants with identical "
              "formula and charge (same composition); otherwise left blank. They are "
              "electronic-energy differences, not free energies -- use G (Hartree) for "
              "thermochemical comparisons.",
              "- Every failed variant keeps its diagnosis in the Failures sheet."]
    return "\n".join(lines) + "\n"
