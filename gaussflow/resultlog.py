"""Append-only real-time result log (JSONL) + recovery into the task DB.

Every finished job and every completed study is appended immediately (with
``flush + fsync``) to ``<home>/logs/results/results_YYYYMMDD.jsonl`` so that a
crash / forced GUI close does not lose computed values.  ``recover`` rebuilds
the result view from those logs (studies are re-created with their summary in
``meta['result']``; analysis falls back to it).
"""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

from .store import Store, now


def log_dir(cfg) -> Path:
    d = Path(cfg.home) / "logs" / "results"
    d.mkdir(parents=True, exist_ok=True)
    return d


def log_path(cfg, day: str | None = None) -> Path:
    day = day or datetime.now().strftime("%Y%m%d")
    return log_dir(cfg) / f"results_{day}.jsonl"


def append_record(cfg, record: dict) -> Path:
    """Append one JSON line and fsync it (real-time durability)."""
    rec = {"ts": now(), **record}
    path = log_path(cfg)
    line = json.dumps(rec, ensure_ascii=False, default=str)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
        fh.flush()
        try:
            os.fsync(fh.fileno())
        except OSError:
            pass
    return path


def iter_records(cfg):
    d = log_dir(cfg)
    for path in sorted(d.glob("results_*.jsonl")):
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue


def summary(cfg) -> dict:
    counts: dict[str, int] = {}
    days = set()
    latest = None
    for rec in iter_records(cfg):
        typ = rec.get("type") or "unknown"
        counts[typ] = counts.get(typ, 0) + 1
        days.add((rec.get("ts") or "")[:10])
        latest = rec.get("ts") or latest
    return {"ok": True, "log_dir": str(log_dir(cfg)),
            "days": sorted(d for d in days if d), "counts": counts,
            "latest": latest}


def logged_study_records(cfg, batch_id: str | None = None) -> dict[str, dict]:
    """Latest ``study_result`` record per study id (later records win)."""
    latest: dict[str, dict] = {}
    for rec in iter_records(cfg):
        if rec.get("type") != "study_result":
            continue
        if batch_id and rec.get("batch_id") != batch_id:
            continue
        sid = rec.get("study_id")
        if sid:
            latest[sid] = rec
    return latest


def missing_studies(store: Store, cfg, batch_id: str | None = None) -> dict:
    """Logged completed studies that are absent from the task database."""
    logged = logged_study_records(cfg, batch_id=batch_id)
    missing = [sid for sid in logged if store.get_study(sid) is None]
    return {"ok": True, "logged": len(logged), "missing": missing,
            "count": len(missing), "log_dir": str(log_dir(cfg))}


def recover(store: Store, cfg, batch_id: str | None = None) -> dict:
    """Recreate studies from study_result records that are missing in the DB."""
    latest = logged_study_records(cfg, batch_id=batch_id)
    recovered, skipped, failed = [], 0, []
    for sid, rec in latest.items():
        if store.get_study(sid) is not None:
            skipped += 1
            continue
        summary = rec.get("summary") or {}
        try:
            variant = summary.get("variant") or rec.get("variant") or {}
            study = store.create_study(
                name=rec.get("name") or summary.get("name") or sid,
                smiles=rec.get("smiles") or summary.get("smiles") or "",
                charge=int(summary.get("charge") or rec.get("charge") or 0),
                multiplicity=int(summary.get("multiplicity") or 1),
                route=summary.get("route") or rec.get("route") or "HF/3-21G",
                workflow=summary.get("workflow") or "opt_freq",
                solvent=summary.get("solvent"),
                preopt=bool(summary.get("preopt")),
                batch_id=rec.get("batch_id") or summary.get("batch_id"),
                status="completed",
                variant=variant or None,
                meta={
                    "result": summary,
                    "recovered": True,
                    "recovered_at": now(),
                    "role": (summary.get("variant") or {}).get("role")
                            or rec.get("role"),
                    "formula": summary.get("formula"),
                    "natoms": summary.get("natoms"),
                    "variant": summary.get("variant"),
                    "li_site": summary.get("li_site"),
                    "redox": summary.get("redox"),
                },
            )
            # align the real study id: keep the log id so result lookups work
            if study["id"] != sid:
                store.conn.execute(
                    "update studies set id=?, batch_id=? where id=?",
                    (sid, study["batch_id"], study["id"]))
                store.conn.commit()
            recovered.append(sid)
        except Exception as exc:  # pragma: no cover - defensive
            failed.append({"study_id": sid, "error": f"{type(exc).__name__}: {exc}"})
    if recovered:
        store.add_event(
            f"recovered {len(recovered)} studies from real-time result logs",
            level="warning")
    return {"ok": True, "recovered": recovered, "skipped": skipped,
            "failed": failed, "count": len(recovered),
            "log_dir": str(log_dir(cfg))}
