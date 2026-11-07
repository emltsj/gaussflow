"""``gflow`` — the controlled tool interface for the Agent.

Each subcommand prints exactly one JSON object when ``--json`` is used, so the
Agent consumes structured data only.  The Agent never writes Gaussian input,
never calls g16.exe directly and never edits output summaries by hand.

Recommended Agent loop:
    gflow doctor --json
    gflow prepare ... --json          # structure + route + warnings, no compute
    gflow submit --study <id> --json  # async, returns job ids
    gflow status --study <id> --json  # poll
    gflow result --study <id> --json  # structured summary
    gflow diagnose --job <id> --json  # only when something failed
    gflow export --batch <id> --json  # CSV / Excel / XYZ / report
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
def _ensure_utf8() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass


def _out(obj, as_json: bool = True) -> None:
    if as_json:
        print(json.dumps(obj, ensure_ascii=False, indent=2, default=str))
        return
    _print_human(obj)


def _print_human(obj) -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, (dict, list)):
                print(f"{k}:")
                print(json.dumps(v, ensure_ascii=False, indent=2, default=str))
            else:
                print(f"{k}: {v}")
    else:
        print(json.dumps(obj, ensure_ascii=False, indent=2, default=str))


def _table(rows: list[dict], columns: list[tuple[str, str]]) -> str:
    headers = [label for _, label in columns]
    data = []
    for r in rows:
        line = []
        for key, _ in columns:
            v = r.get(key)
            if isinstance(v, float):
                v = f"{v:.6f}"
            elif isinstance(v, (list, dict)):
                v = json.dumps(v, ensure_ascii=False)
            line.append("" if v is None else v)
        data.append(line)
    try:
        from tabulate import tabulate
        return tabulate(data, headers=headers, tablefmt="simple", missingval="")
    except Exception:  # pragma: no cover
        out = ["\t".join(headers)]
        out += ["\t".join(map(str, row)) for row in data]
        return "\n".join(out)


def _cfg(args):
    from .config import Config
    return Config.load(getattr(args, "config", None))


def _apply_model(cfg, args, *, save: bool = True) -> str:
    """Apply an explicit --model (v4pro / v4flash / full id) to the config."""
    from . import llm
    name = getattr(args, "model", None)
    if name:
        cfg.set("llm.model", llm.resolve_model(name, cfg))
        if save:
            cfg.save()
    return str(cfg.get("llm.model") or llm.MODEL_FLASH)


def _store(cfg):
    from .store import Store
    return Store(cfg.db_path)


def _auto_name(args) -> str:
    if getattr(args, "name", None):
        return args.name
    if getattr(args, "smiles", None):
        base = re.sub(r"[^0-9A-Za-z]+", "_", args.smiles).strip("_")
        return base[:32] or "molecule"
    if getattr(args, "file", None):
        return Path(args.file).stem
    return "molecule"


def _split_ions(text: str | None) -> list[str]:
    if not text:
        return []
    return [x.strip() for x in re.split(r"[,\s]+", text) if x.strip()]


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------
def cmd_init(args) -> int:
    cfg = _cfg(args)
    store = _store(cfg)
    store.close()
    _out({"ok": True, "config": str(cfg.path), **cfg.describe()}, args.json)
    return 0


def cmd_doctor(args) -> int:
    cfg = _cfg(args)
    checks: dict = {}
    problems: list[str] = []

    problems += cfg.require_binaries()
    for label, p in (("g16_exe", cfg.g16_exe), ("gview_exe", cfg.gview_exe),
                     ("formchk_exe", cfg.formchk_exe)):
        checks[label] = {"path": str(p), "exists": Path(p).exists()}

    import shutil
    free_gb = round(shutil.disk_usage(str(cfg.work_root)).free / 1e9, 1)
    checks["disk_free_gb"] = free_gb
    if free_gb < 5:
        problems.append(f"磁盘剩余空间不足: {free_gb} GB")

    for mod, label in (("rdkit", "rdkit"), ("cclib", "cclib"),
                       ("openpyxl", "openpyxl"), ("numpy", "numpy")):
        try:
            m = __import__(mod)
            checks[label] = getattr(m, "__version__", "ok")
        except Exception as exc:
            checks[label] = f"MISSING ({exc})"
            problems.append(f"Python 依赖缺失: {label}")

    try:
        store = _store(cfg)
        store.close()
        checks["db_writable"] = True
    except Exception as exc:
        checks["db_writable"] = False
        problems.append(f"数据库不可写: {exc}")

    from .daemon import daemon_status
    checks["worker"] = daemon_status(cfg)

    if args.test_run:
        checks["g16_test"] = _doctor_test_run(cfg, problems)

    _out({"ok": not problems, "problems": problems, "checks": checks,
          "config": str(cfg.path)}, args.json)
    return 0 if not problems else 1


def _doctor_test_run(cfg, problems: list[str]) -> dict:
    from .parse import parse_output
    from .runner import run_job
    d = cfg.work_root / "_doctor" / f"t{int(__import__('time').time())}"
    d.mkdir(parents=True, exist_ok=True)
    (d / "doctor.gjf").write_text(
        "%chk=doctor.chk\n%mem=1GB\n%nproc=2\n#p HF/STO-3G sp\n\n"
        "gaussflow doctor\n\n0 1\n"
        "O   0.000000   0.000000   0.117300\n"
        "H   0.000000   0.757200  -0.469200\n"
        "H   0.000000  -0.757200  -0.469200\n\n",
        encoding="ascii")
    job = {"dir": str(d), "input_path": str(d / "doctor.gjf"),
           "output_path": str(d / "doctor.out")}
    try:
        outcome = run_job(cfg, job)
        parsed = parse_output(outcome["output_path"])
        ok = bool(parsed.get("normal_termination"))
        if not ok:
            problems.append("Gaussian 测试任务未正常结束（检查许可/安装）")
        return {"ok": ok, "exit_code": outcome.get("exit_code"),
                "walltime_s": outcome.get("walltime"),
                "energy_hartree": parsed.get("final_energy_hartree"),
                "version": parsed.get("version"),
                "output": outcome.get("output_path")}
    except Exception as exc:
        problems.append(f"Gaussian 测试任务异常: {exc}")
        return {"ok": False, "error": str(exc)}


def cmd_prepare(args) -> int:
    from . import pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    preopt = None
    if args.no_preopt:
        preopt = False
    study = pipeline.create_study(
        store, cfg,
        name=_auto_name(args),
        smiles=args.smiles, file=args.file,
        route=args.route, workflow=args.workflow,
        charge=args.charge, multiplicity=args.multiplicity,
        solvent=args.solvent, preopt=preopt, seed=args.seed,
        n_conformers=args.conformers, ion=args.ion, sp_route=args.sp_route,
    )
    meta = study.get("meta") or {}
    result = {
        "ok": True,
        "study_id": study["id"],
        "name": study["name"],
        "status": study["status"],
        "natoms": meta.get("natoms"),
        "formula": meta.get("formula"),
        "charge": study["charge"],
        "multiplicity": study["multiplicity"],
        "workflow": study["workflow"],
        "preopt": bool(study["preopt"]),
        "route": study.get("route"),
        "route_preview": meta.get("route_preview"),
        "route_check": meta.get("route_check"),
        "needs_approval": meta.get("needs_approval"),
        "approval_reasons": meta.get("approval_reasons"),
        "run_dir": str(cfg.work_root / study["id"]),
        "gaussview_preview_file": str(cfg.work_root / study["id"] / "initial.gjf"),
        "next": f"gflow submit --study {study['id']}",
    }
    _out(result, args.json)
    return 0 if not meta.get("needs_approval") else 3


def cmd_submit(args) -> int:
    from . import daemon, pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    res = pipeline.submit_study(store, cfg, args.study, approve=args.approve)
    if res.get("status") == "submitted" and not args.no_daemon:
        res["worker"] = daemon.start_daemon(cfg)
    _out(res, args.json)
    return 0 if res.get("ok") and res.get("status") == "submitted" else (
        3 if res.get("status") == "needs_approval" else 1)


def cmd_submit_batch(args) -> int:
    from . import daemon, pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    res = pipeline.submit_batch(store, cfg, args.batch, approve=args.approve)
    if not args.no_daemon:
        res["worker"] = daemon.start_daemon(cfg)
    _out(res, args.json)
    return 0 if res.get("ok") else 1


def cmd_status(args) -> int:
    from .daemon import daemon_status
    cfg = _cfg(args)
    store = _store(cfg)
    if args.study:
        s = store.get_study(args.study)
        if s is None:
            _out({"ok": False, "error": f"未知 study: {args.study}"}, args.json)
            return 1
        studies = [s]
    elif args.batch:
        studies = store.list_studies(batch_id=args.batch)
    else:
        studies = store.list_studies(limit=args.limit)

    rows = []
    for s in studies:
        jobs = store.study_jobs(s["id"])
        running = [j for j in jobs if j["status"] == "running"]
        pending = [j for j in jobs if j["status"] == "pending"]
        last_done = [j for j in jobs if j["status"] in ("completed", "failed")]
        current = (running or pending or last_done[-1:] or [None])[0]
        res = store.get_result(current["id"]) if current else None
        parsed = (res or {}).get("parsed") or {}
        validation = (res or {}).get("validation") or {}
        rows.append({
            "study_id": s["id"],
            "name": s["name"],
            "status": s["status"],
            "variant": ((s.get("meta") or {}).get("variant") or {}).get("label", "base"),
            "stage": current["stage"] if current else None,
            "job_status": current["status"] if current else None,
            "energy_hartree": parsed.get("final_energy_hartree"),
            "validation": validation.get("overall"),
            "n_imag": (parsed.get("freq") or {}).get("n_imag"),
            "updated_at": s.get("updated_at"),
        })
    result = {"ok": True, "count": len(rows), "studies": rows,
              "worker": daemon_status(cfg),
              "pending_jobs": store.pending_count()}
    if args.json:
        _out(result, True)
    else:
        print(_table(rows, [("study_id", "Study ID"), ("name", "Name"),
                            ("status", "Status"), ("variant", "Variant"),
                            ("stage", "Stage"), ("job_status", "Job"),
                            ("energy_hartree", "E (Hartree)"),
                            ("validation", "Valid"),
                            ("n_imag", "NImag")]))
        print(f"\nworker: {result['worker']}  pending jobs: {result['pending_jobs']}")
    return 0


def cmd_result(args) -> int:
    from . import pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    summary = pipeline.study_summary(store, cfg, args.study)
    _out(summary, args.json)
    return 0 if summary.get("ok") else 1


def cmd_diagnose(args) -> int:
    from .diagnose import diagnose
    from .parse import parse_output
    from .validate import validate
    cfg = _cfg(args)
    store = _store(cfg)
    job = store.get_job(args.job)
    if job is None:
        _out({"ok": False, "error": f"未知 job: {args.job}"}, args.json)
        return 1
    res = store.get_result(job["id"]) or {}
    parsed = res.get("parsed")
    if not parsed:
        parsed = parse_output(job.get("output_path") or "")
    validation = res.get("validation") or validate(parsed, job, cfg)
    diag = diagnose(parsed, job, validation, cfg,
                    output_path=job.get("output_path"))
    _out(diag, args.json)
    return 0


def cmd_export(args) -> int:
    from . import export as exporter
    cfg = _cfg(args)
    store = _store(cfg)
    if args.batch:
        res = exporter.export_batch(store, cfg, args.batch, out_dir=args.out)
    elif args.study:
        res = exporter.export_study(store, cfg, args.study)
    else:
        studies = store.list_studies(limit=args.limit)
        out = Path(args.out) if args.out else cfg.exports_dir / "all"
        out.mkdir(parents=True, exist_ok=True)
        res = {"ok": True, "out_dir": str(out), "studies": []}
        for s in studies:
            res["studies"].append(exporter.export_study(store, cfg, s["id"]))
    _out(res, args.json)
    return 0 if res.get("ok") else 1


def cmd_iterate(args) -> int:
    from . import daemon, pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    res = pipeline.create_batch(
        store, cfg, base_name=_auto_name(args),
        smiles=args.smiles, file=args.file, n=args.n,
        ions=_split_ions(args.ions), seed=args.seed,
        route=args.route, workflow=args.workflow, solvent=args.solvent,
        sp_route=args.sp_route,
        preopt=False if args.no_preopt else None,
        include_base=not args.no_base,
    )
    if args.submit:
        res["submit"] = pipeline.submit_batch(store, cfg, res["batch_id"],
                                              approve=args.approve)
        res["worker"] = daemon.start_daemon(cfg)
    res["next"] = (f"gflow submit-batch --batch {res['batch_id']}" if not args.submit
                   else f"gflow status --batch {res['batch_id']}")
    _out(res, args.json)
    return 0


def cmd_daemon(args) -> int:
    from . import daemon
    cfg = _cfg(args)
    if args.action == "start":
        res = daemon.start_daemon(cfg)
    elif args.action == "stop":
        res = daemon.stop_daemon(cfg)
    else:
        res = daemon.daemon_status(cfg)
    _out(res, args.json)
    return 0 if res.get("ok", True) else 1


def cmd_run(args) -> int:
    from .daemon import Worker
    cfg = _cfg(args)
    worker = Worker(cfg)
    if args.once:
        did = worker.process_one()
        _out({"ok": True, "processed": bool(did),
              "pending": worker.store.pending_count()}, args.json)
        return 0
    worker.loop()
    return 0


def cmd_preview(args) -> int:
    cfg = _cfg(args)
    store = _store(cfg)
    if args.job:
        job = store.get_job(args.job)
        if job is None:
            _out({"ok": False, "error": f"未知 job: {args.job}"}, args.json)
            return 1
        target = Path(args.final and job.get("output_path") or job.get("input_path") or "")
    else:
        study = store.get_study(args.study)
        if study is None:
            _out({"ok": False, "error": f"未知 study: {args.study}"}, args.json)
            return 1
        if args.final:
            from .pipeline import study_summary
            summary = study_summary(store, cfg, args.study)
            fid = summary.get("final_job_id")
            job = store.get_job(fid) if fid else None
            target = Path(job.get("output_path")) if job else None
        else:
            target = cfg.work_root / args.study / "initial.gjf"
    if target is None or not Path(target).exists():
        _out({"ok": False, "error": f"文件不存在: {target}"}, args.json)
        return 1
    if not cfg.gview_exe:
        _out({"ok": False, "error": "GaussView is not configured; run `gflow detect` "
              "or set gview_exe in config.json."}, args.json)
        return 1
    proc = subprocess.Popen([str(cfg.gview_exe), str(target)])
    _out({"ok": True, "launched": "GaussView", "pid": proc.pid,
          "file": str(target),
          "note": "GaussView 是交互式工具，仅用于人工检查；流水线不依赖它"},
         args.json)
    return 0


def cmd_list(args) -> int:
    cfg = _cfg(args)
    store = _store(cfg)
    studies = store.list_studies(limit=args.limit)
    rows = []
    for s in studies:
        meta = s.get("meta") or {}
        rows.append({
            "study_id": s["id"], "name": s["name"], "status": s["status"],
            "batch_id": s.get("batch_id"),
            "variant": (meta.get("variant") or {}).get("label", "base"),
            "charge": s.get("charge"), "route": s.get("route"),
            "created_at": s.get("created_at"),
        })
    if args.json:
        _out({"ok": True, "studies": rows}, True)
    else:
        print(_table(rows, [("study_id", "Study ID"), ("name", "Name"),
                            ("status", "Status"), ("batch_id", "Batch"),
                            ("variant", "Variant"), ("charge", "Q"),
                            ("route", "Route"), ("created_at", "Created")]))
    return 0


def cmd_jobs(args) -> int:
    cfg = _cfg(args)
    store = _store(cfg)
    jobs = store.study_jobs(args.study)
    rows = []
    for j in jobs:
        res = store.get_result(j["id"]) or {}
        parsed = res.get("parsed") or {}
        validation = res.get("validation") or {}
        diag = (j.get("meta") or {}).get("diagnosis") or {}
        rows.append({
            "job_id": j["id"], "stage": j["stage"], "attempt": j.get("attempt"),
            "status": j["status"], "walltime_s": j.get("walltime"),
            "validation": validation.get("overall"),
            "error_class": j.get("error_class") or diag.get("class"),
            "energy_hartree": parsed.get("final_energy_hartree"),
            "output": j.get("output_path"),
        })
    if args.json:
        _out({"ok": True, "jobs": rows}, True)
    else:
        print(_table(rows, [("job_id", "Job ID"), ("stage", "Stage"),
                            ("attempt", "Try"), ("status", "Status"),
                            ("walltime_s", "Wall time"), ("validation", "Valid"),
                            ("error_class", "Error"),
                            ("energy_hartree", "E (Hartree)")]))
    return 0


def cmd_detect(args) -> int:
    from .locate import detect_all
    found = detect_all({
        "g16_exe": None, "gview_exe": None, "shermo_exe": None,
    })
    _out({"ok": True, "found": found["found"],
          "gaussian": found["gaussian"]["candidates"],
          "gaussview": found["gaussview"]["candidates"],
          "shermo": found["shermo"]["candidates"]}, args.json)
    return 0


def cmd_recipe(args) -> int:
    from . import daemon, pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    options = {
        "li": None if args.li is None else bool(args.li),
        "li_elements": _split_ions(args.li_elements) or None,
        "redox": None if args.redox is None else bool(args.redox),
        "preopt": None if args.preopt is None else bool(args.preopt),
        "solvent": args.solvent,
    }
    res = pipeline.create_recipe_system(
        store, cfg, name=_auto_name(args),
        smiles=args.smiles, file=args.file,
        recipe_name=args.recipe, options=options,
        charge=args.charge, multiplicity=args.multiplicity)
    if args.submit:
        res["submit"] = pipeline.submit_batch(store, cfg, res["batch_id"],
                                              approve=args.approve)
        res["worker"] = daemon.start_daemon(cfg)
    else:
        res["next"] = f"gflow submit-batch --batch {res['batch_id']}"
    _out(res, args.json)
    return 0


def cmd_analysis(args) -> int:
    from . import analysis
    cfg = _cfg(args)
    store = _store(cfg)
    _out(analysis.system_summary(store, cfg, args.batch), args.json)
    return 0


def cmd_export_system(args) -> int:
    from . import analysis
    cfg = _cfg(args)
    store = _store(cfg)
    res = analysis.export_system(store, cfg, args.batch, out_dir=args.out)
    _out(res, args.json)
    return 0 if res.get("ok") else 1


def cmd_list_recipes(args) -> int:
    from .recipes import list_recipes
    cfg = _cfg(args)
    recipes = cfg.get("recipes") or {}
    items = [{"key": k, "name": v.get("name", k),
              "description": v.get("description", ""),
              "stages": [s.get("name") for s in v.get("core_stages", [])],
              "li": bool((v.get("li") or {}).get("on")),
              "redox": bool((v.get("redox") or {}).get("on")),
              "ai_generated": bool(v.get("ai_generated"))}
             for k, v in recipes.items()]
    _out({"ok": True, "recipes": items}, args.json)
    return 0


def cmd_ai(args) -> int:
    from . import daemon, fluorinate, llm, pipeline
    cfg = _cfg(args)
    _apply_model(cfg, args)
    store = _store(cfg)
    res = llm.generate_molecules(
        cfg, args.prompt, f_range=args.f_range,
        n_molecules=getattr(args, "n_molecules", 1) or 1)
    enum_counts = llm.parse_f_range(getattr(args, "enum_f", None))
    if enum_counts:
        molecules = []
        maxv = int(cfg.get("llm.max_fluorination_variants") or 20)
        notes = []
        for base in res["molecules"]:
            built, base_notes = fluorinate.molecules_for_f_counts(
                base, enum_counts, max_variants=maxv)
            molecules.extend(built)
            notes.extend(base_notes)
        if not molecules:
            _out({"ok": False, "error": "按位点枚举 F 没有得到结构",
                  "notes": notes}, args.json)
            return 1
        res["molecules"] = molecules
        res["enumeration"] = {"n_f_values": enum_counts, "max_variants": maxv,
                              "notes": notes}
    preview = pipeline.export_gjf_previews(
        cfg, res["molecules"], recipe_name=args.recipe,
        name_prefix=args.name or "AI")
    out = {
        "ok": True,
        "n_molecules": len(res["molecules"]),
        "molecules": res["molecules"],
        "explanation": res.get("explanation"),
        "missing_f_counts": res.get("missing_f_counts"),
        "enumeration": res.get("enumeration"),
        "model": res.get("model"),
        "gjf_dir": preview["gjf_dir"],
        "files": preview["files"],
    }
    if args.submit:
        out["calc"] = pipeline.create_ai_systems(
            store, cfg, res["molecules"], recipe_name=args.recipe,
            options={"li": None if args.li is None else bool(args.li),
                     "redox": None if args.redox is None else bool(args.redox),
                     "preopt": None if args.preopt is None else bool(args.preopt),
                     "solvent": args.solvent},
            name_prefix=args.name or "AI", submit=True,
            approve=args.approve)
        out["worker"] = daemon.start_daemon(cfg)
    _out(out, args.json)
    return 0


def cmd_fluorinate(args) -> int:
    """Enumerate symmetry-unique F substitution positions (deterministic)."""
    from . import daemon, pipeline
    from .fluorinate import (describe_fluorination_sites,
                             enumerate_fluorinated_variants,
                             molecules_from_variants)
    cfg = _cfg(args)
    store = _store(cfg)
    if args.describe:
        sites = describe_fluorination_sites(args.smiles)
        _out({"ok": True, "smiles": args.smiles, "sites": sites}, args.json)
        return 0
    n = int(args.n_f or 1)
    maxv = int(args.max_variants or cfg.get("llm.max_fluorination_variants") or 20)
    variants = enumerate_fluorinated_variants(args.smiles, n_f=n, max_variants=maxv)
    molecules = molecules_from_variants(
        args.name or "mol", args.smiles, variants,
        charge=int(args.charge or 0), multiplicity=int(args.multiplicity or 1))
    preview = pipeline.export_gjf_previews(
        cfg, molecules, recipe_name=args.recipe, name_prefix=args.name or "F-enum")
    out = {"ok": True, "smiles": args.smiles, "n_f": n,
           "n_variants": len(molecules), "molecules": molecules,
           "gjf_dir": preview["gjf_dir"], "files": preview["files"]}
    if args.submit:
        out["calc"] = pipeline.create_ai_systems(
            store, cfg, molecules, recipe_name=args.recipe,
            options={"li": None if args.li is None else bool(args.li),
                     "redox": None if args.redox is None else bool(args.redox)},
            name_prefix=args.name or "F-enum", submit=True,
            approve=args.approve)
        out["worker"] = daemon.start_daemon(cfg)
    _out(out, args.json)
    return 0


def cmd_ai_workflow(args) -> int:
    """Generate a calculation-chain recipe with DeepSeek and store it."""
    from . import llm, pipeline
    cfg = _cfg(args)
    _apply_model(cfg, args)
    res = llm.generate_workflow(cfg, args.prompt)
    reg = pipeline.register_ai_recipe(cfg, res["recipe"], key=args.key)
    out = {"ok": bool(reg.get("ok")), "key": reg.get("key"),
           "recipe": reg.get("recipe"), "problems": reg.get("problems"),
           "explanation": res.get("explanation"),
           "model": res.get("model"), "usage": res.get("usage")}
    _out(out, args.json)
    return 0 if reg.get("ok") else 1


def cmd_recover(args) -> int:
    """Restore the result view from the real-time JSONL result log."""
    from . import resultlog
    cfg = _cfg(args)
    store = _store(cfg)
    res = resultlog.recover(store, cfg, batch_id=args.batch)
    res["log"] = resultlog.summary(cfg)
    _out(res, args.json)
    return 0


def cmd_logs(args) -> int:
    """Show where the real-time result log is and what it contains."""
    from . import resultlog
    cfg = _cfg(args)
    _out(resultlog.summary(cfg), args.json)
    return 0


def cmd_gjf_batch(args) -> int:
    from . import daemon, pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    paths = [p.strip().strip('"') for p in re.split(r"[;,\n]+", args.files or "")
             if p.strip()]
    if args.dir:
        paths += [str(p) for p in Path(args.dir).glob("*.gjf")]
        paths += [str(p) for p in Path(args.dir).glob("*.com")]
    res = pipeline.create_batch_from_gjfs(
        store, cfg, paths, name_prefix=args.name, submit=args.submit)
    if args.submit:
        res["worker"] = daemon.start_daemon(cfg)
    _out(res, args.json)
    return 0


def cmd_llm_test(args) -> int:
    from .llm import chat
    cfg = _cfg(args)
    _apply_model(cfg, args)
    try:
        res = chat(cfg, [{"role": "user",
                          "content": '只输出 JSON：{"ok": true, "msg": "pong"}'}])
    except Exception as exc:
        _out({"ok": False, "error": f"{type(exc).__name__}: {exc}"}, args.json)
        return 1
    _out({"ok": True, "model": res.get("model"),
          "content": (res.get("content") or "")[:500],
          "usage": res.get("usage")}, args.json)
    return 0


def cmd_stop(args) -> int:
    """Force-stop running Gaussian jobs and cancel the queue."""
    from . import daemon, pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    if args.batch:
        res = pipeline.cancel_jobs(store, cfg, batch_ids=[args.batch])
    elif args.all:
        res = pipeline.cancel_jobs(store, cfg, all_jobs=True)
    else:
        _out({"ok": False, "error": "需要 --batch <id> 或 --all"}, args.json)
        return 1
    if args.daemon or args.all:
        res["daemon"] = daemon.stop_daemon(cfg, timeout=8)
    _out(res, args.json)
    return 0


def cmd_purge(args) -> int:
    """Remove cancelled studies from the task list."""
    from . import pipeline
    cfg = _cfg(args)
    store = _store(cfg)
    batch_ids = [args.batch] if args.batch else None
    res = pipeline.purge_cancelled(store, cfg, batch_ids=batch_ids,
                                   delete_files=args.delete_files)
    _out(res, args.json)
    return 0


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gflow",
        description="Agent-driven Gaussian 16 pipeline (controlled tool interface)")
    p.add_argument("--config", help="config.json 路径（默认 GFLOW_CONFIG 或项目目录）")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_json(sp):
        sp.add_argument("--json", action="store_true", help="输出单个 JSON 对象")

    sp = sub.add_parser("init", help="初始化配置与数据库")
    add_json(sp); sp.set_defaults(func=cmd_init)

    sp = sub.add_parser("doctor", help="环境自检（可含 Gaussian 试算）")
    sp.add_argument("--test-run", action="store_true", help="运行一个水分子测试任务")
    add_json(sp); sp.set_defaults(func=cmd_doctor)

    sp = sub.add_parser("prepare", help="生成结构+输入预览（不提交计算）")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles")
    g.add_argument("--file", help=".sdf/.mol/.xyz/.gjf 结构文件")
    sp.add_argument("--name")
    sp.add_argument("--route", help="主计算级别，如 'B3LYP/6-31G(d)'")
    sp.add_argument("--workflow", choices=["opt", "opt_freq", "opt_freq_sp", "sp", "freq"])
    sp.add_argument("--charge", type=int)
    sp.add_argument("--multiplicity", type=int)
    sp.add_argument("--solvent", help="隐式溶剂，如 water")
    sp.add_argument("--sp-route", help="opt_freq_sp 工作流中单点能级别，如 'B3LYP/6-311+G(d,p)'")
    sp.add_argument("--ion", help="给分子加一个离子，如 Na+ / Cl-")
    sp.add_argument("--no-preopt", action="store_true", help="跳过 HF/3-21G 预优化")
    sp.add_argument("--conformers", type=int, default=1, help="构象搜索数量")
    sp.add_argument("--seed", type=int, default=42)
    add_json(sp); sp.set_defaults(func=cmd_prepare)

    sp = sub.add_parser("submit", help="提交 study（异步）")
    sp.add_argument("--study", required=True)
    sp.add_argument("--approve", action="store_true", help="批准高风险/昂贵任务")
    sp.add_argument("--no-daemon", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_submit)

    sp = sub.add_parser("submit-batch", help="提交整个 batch")
    sp.add_argument("--batch", required=True)
    sp.add_argument("--approve", action="store_true")
    sp.add_argument("--no-daemon", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_submit_batch)

    sp = sub.add_parser("status", help="查询状态")
    g = sp.add_mutually_exclusive_group()
    g.add_argument("--study")
    g.add_argument("--batch")
    sp.add_argument("--limit", type=int, default=30)
    add_json(sp); sp.set_defaults(func=cmd_status)

    sp = sub.add_parser("result", help="结构化结果摘要")
    sp.add_argument("--study", required=True)
    add_json(sp); sp.set_defaults(func=cmd_result)

    sp = sub.add_parser("diagnose", help="失败诊断 + 白名单重试建议")
    sp.add_argument("--job", required=True)
    add_json(sp); sp.set_defaults(func=cmd_diagnose)

    sp = sub.add_parser("export", help="导出 CSV/Excel/XYZ/报告")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--batch")
    g.add_argument("--study")
    g.add_argument("--all", action="store_true")
    sp.add_argument("--out")
    sp.add_argument("--limit", type=int, default=200)
    add_json(sp); sp.set_defaults(func=cmd_export)

    sp = sub.add_parser("iterate", help="随机迭代：生成加离子等变体 batch")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles")
    g.add_argument("--file")
    sp.add_argument("--name")
    sp.add_argument("--n", type=int, default=3, help="变体数量")
    sp.add_argument("--ions", help="候选离子，逗号分隔，如 'Na+,Cl-,K+'")
    sp.add_argument("--no-base", action="store_true", help="不包含原始分子")
    sp.add_argument("--route")
    sp.add_argument("--sp-route", help="opt_freq_sp 的单点能级别")
    sp.add_argument("--workflow", choices=["opt", "opt_freq", "opt_freq_sp", "sp", "freq"])
    sp.add_argument("--solvent")
    sp.add_argument("--no-preopt", action="store_true")
    sp.add_argument("--seed", type=int, default=42)
    sp.add_argument("--submit", action="store_true", help="创建后立即提交")
    sp.add_argument("--approve", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_iterate)

    sp = sub.add_parser("daemon", help="后台 worker 管理")
    sp.add_argument("action", choices=["start", "stop", "status"])
    add_json(sp); sp.set_defaults(func=cmd_daemon)

    sp = sub.add_parser("run", help="前台运行 worker（调试用）")
    sp.add_argument("--once", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("preview", help="用 GaussView 打开结构/结果（人工检查）")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--study")
    g.add_argument("--job")
    sp.add_argument("--final", action="store_true", help="打开最终结果而非初始结构")
    add_json(sp); sp.set_defaults(func=cmd_preview)

    sp = sub.add_parser("list", help="列出 studies")
    sp.add_argument("--limit", type=int, default=50)
    add_json(sp); sp.set_defaults(func=cmd_list)

    sp = sub.add_parser("jobs", help="列出某个 study 的所有 job")
    sp.add_argument("--study", required=True)
    add_json(sp); sp.set_defaults(func=cmd_jobs)

    # ---- calculation-system / GUI-oriented commands --------------------
    sp = sub.add_parser("detect", help="自动检测本机 Gaussian / GaussView / Shermo")
    add_json(sp); sp.set_defaults(func=cmd_detect)

    sp = sub.add_parser("recipes", help="列出预设计算方案")
    add_json(sp); sp.set_defaults(func=cmd_list_recipes)

    sp = sub.add_parser("recipe", help="创建计算体系（opt+freq / SP / Li⁺+BSSE / 氧化还原）")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles")
    g.add_argument("--file")
    sp.add_argument("--name")
    sp.add_argument("--recipe", default="opt_freq",
                    help="方案 key（见 gflow recipes）")
    sp.add_argument("--li", dest="li", action="store_true", default=None,
                    help="计算 Li⁺ 配合物")
    sp.add_argument("--no-li", dest="li", action="store_false")
    sp.add_argument("--li-elements", default="O,F", help="Li⁺ 配位元素，如 O,F")
    sp.add_argument("--redox", dest="redox", action="store_true", default=None,
                    help="计算氧化态/还原态")
    sp.add_argument("--no-redox", dest="redox", action="store_false")
    sp.add_argument("--preopt", dest="preopt", action="store_true", default=None)
    sp.add_argument("--no-preopt", dest="preopt", action="store_false")
    sp.add_argument("--solvent")
    sp.add_argument("--charge", type=int)
    sp.add_argument("--multiplicity", type=int)
    sp.add_argument("--submit", action="store_true")
    sp.add_argument("--approve", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_recipe)

    sp = sub.add_parser("analysis", help="体系级分析：Li⁺结合能/氧化还原/HOMO-LUMO/热化学")
    sp.add_argument("--batch", required=True)
    add_json(sp); sp.set_defaults(func=cmd_analysis)

    sp = sub.add_parser("export-system", help="导出体系分析 Excel/CSV/报告")
    sp.add_argument("--batch", required=True)
    sp.add_argument("--out")
    add_json(sp); sp.set_defaults(func=cmd_export_system)

    # ---- AI generation / multi-gjf queue -------------------------------
    sp = sub.add_parser("ai", help="用 DeepSeek 按自然语言生成分子结构并导出 gjf")
    sp.add_argument("--prompt", required=True, help="分子描述，例如 '生成 TEP 溶剂，加 1 个 F'")
    sp.add_argument("--f-range", help="F 原子个数范围，如 '0-3'（生成系列分子）")
    sp.add_argument("--n-molecules", type=int, default=1,
                    help="一次生成多少个不同候选分子（默认 1）")
    sp.add_argument("--enum-f", default=None,
                    help="按对称不等价位点枚举 F 取代，支持个数或范围：1 / 0-3")
    sp.add_argument("--name", help="体系名称前缀")
    sp.add_argument("--recipe", default="opt_freq")
    sp.add_argument("--li", dest="li", action="store_true", default=None)
    sp.add_argument("--no-li", dest="li", action="store_false")
    sp.add_argument("--redox", dest="redox", action="store_true", default=None)
    sp.add_argument("--no-redox", dest="redox", action="store_false")
    sp.add_argument("--preopt", dest="preopt", action="store_true", default=None)
    sp.add_argument("--no-preopt", dest="preopt", action="store_false")
    sp.add_argument("--solvent")
    sp.add_argument("--model", help="DeepSeek 模型：v4flash（默认）/ v4pro / 完整 id")
    sp.add_argument("--submit", action="store_true")
    sp.add_argument("--approve", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_ai)

    sp = sub.add_parser("fluorinate", help="按对称不等价位点枚举 F 取代（RDKit，确定性）")
    sp.add_argument("--smiles", required=True, help="基准分子 SMILES")
    sp.add_argument("--n-f", type=int, default=1, help="替换的 F 个数（默认 1）")
    sp.add_argument("--max-variants", type=int, default=20)
    sp.add_argument("--name", help="名称前缀")
    sp.add_argument("--charge", type=int, default=0)
    sp.add_argument("--multiplicity", type=int, default=1)
    sp.add_argument("--describe", action="store_true",
                    help="只列出对称不等价的单 F 取代位点")
    sp.add_argument("--recipe", default="opt_freq")
    sp.add_argument("--li", dest="li", action="store_true", default=None)
    sp.add_argument("--no-li", dest="li", action="store_false")
    sp.add_argument("--redox", dest="redox", action="store_true", default=None)
    sp.add_argument("--no-redox", dest="redox", action="store_false")
    sp.add_argument("--submit", action="store_true")
    sp.add_argument("--approve", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_fluorinate)

    sp = sub.add_parser("ai-workflow", help="用 DeepSeek 生成自定义计算链并保存为方案")
    sp.add_argument("--prompt", required=True, help="计算链描述（方法/基组/溶剂/Li+/氧化还原）")
    sp.add_argument("--key", help="保存的方案 key（默认自动生成 ai_*）")
    sp.add_argument("--model", help="DeepSeek 模型：v4flash（默认）/ v4pro / 完整 id")
    add_json(sp); sp.set_defaults(func=cmd_ai_workflow)

    sp = sub.add_parser("recover", help="从实时结果日志恢复结果显示")
    sp.add_argument("--batch", help="只恢复该批次")
    add_json(sp); sp.set_defaults(func=cmd_recover)

    sp = sub.add_parser("logs", help="查看实时结果日志位置与内容统计")
    add_json(sp); sp.set_defaults(func=cmd_logs)

    sp = sub.add_parser("gjf-batch", help="批量导入 gjf 并排队计算（可多选/整目录）")
    sp.add_argument("--files", help="gjf 路径，逗号/分号分隔")
    sp.add_argument("--dir", help="目录（导入其中所有 .gjf/.com）")
    sp.add_argument("--name", help="名称前缀")
    sp.add_argument("--submit", action="store_true")
    add_json(sp); sp.set_defaults(func=cmd_gjf_batch)

    sp = sub.add_parser("llm-test", help="测试 DeepSeek API 连接")
    sp.add_argument("--model", help="DeepSeek 模型：v4flash（默认）/ v4pro / 完整 id")
    add_json(sp); sp.set_defaults(func=cmd_llm_test)

    sp = sub.add_parser("stop", help="强制终止正在运行的 Gaussian 并取消排队任务")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--batch", help="只停止该批次")
    g.add_argument("--all", action="store_true", help="停止所有批次")
    sp.add_argument("--daemon", action="store_true", help="同时停止后台 worker")
    add_json(sp); sp.set_defaults(func=cmd_stop)

    sp = sub.add_parser("purge", help="从任务列表清除已取消的研究记录")
    g = sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--batch", help="只清理该批次")
    g.add_argument("--all", action="store_true", help="清理所有批次的取消记录")
    sp.add_argument("--delete-files", action="store_true",
                    help="同时删除 runs/ 下对应任务目录（默认保留文件）")
    add_json(sp); sp.set_defaults(func=cmd_purge)

    return p


def main(argv=None) -> int:
    _ensure_utf8()
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130
    except Exception as exc:
        import traceback
        print(json.dumps({"ok": False, "error": f"{type(exc).__name__}: {exc}",
                          "traceback": traceback.format_exc()[-1500:]},
                         ensure_ascii=False, indent=2), file=sys.stdout)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
