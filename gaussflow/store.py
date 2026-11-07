"""SQLite-backed study/job store.

One database per gaussflow home (WAL mode) so that the CLI, the worker
daemon and ad-hoc analysis scripts can share state safely.  All scientific
payloads (parsed results, validation reports) are stored as JSON text so a
new parser/validator version can be applied by simply re-parsing outputs,
without a schema migration of every field.
"""
from __future__ import annotations

import json
import secrets
import sqlite3
from datetime import datetime
from pathlib import Path

_SCHEMA = """
CREATE TABLE IF NOT EXISTS studies (
    id           TEXT PRIMARY KEY,
    batch_id     TEXT,
    name         TEXT NOT NULL,
    smiles       TEXT,
    variant      TEXT,
    charge       INTEGER,
    multiplicity INTEGER,
    route        TEXT,
    workflow     TEXT,
    solvent      TEXT,
    preopt       INTEGER DEFAULT 1,
    status       TEXT NOT NULL,
    spec         TEXT,
    meta         TEXT,
    created_at   TEXT NOT NULL,
    updated_at   TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_studies_batch ON studies(batch_id);
CREATE INDEX IF NOT EXISTS idx_studies_status ON studies(status);

CREATE TABLE IF NOT EXISTS jobs (
    id           TEXT PRIMARY KEY,
    study_id     TEXT NOT NULL,
    stage        TEXT NOT NULL,
    seq          INTEGER NOT NULL,
    attempt      INTEGER NOT NULL DEFAULT 1,
    lane         INTEGER NOT NULL DEFAULT 0,
    status       TEXT NOT NULL,
    dir          TEXT,
    route        TEXT,
    charge       INTEGER,
    multiplicity INTEGER,
    input_path   TEXT,
    output_path  TEXT,
    chk_path     TEXT,
    pid          INTEGER,
    exit_code    INTEGER,
    timed_out    INTEGER DEFAULT 0,
    walltime     REAL,
    error_class  TEXT,
    parent_job   TEXT,
    retry_of     TEXT,
    meta         TEXT,
    created_at   TEXT NOT NULL,
    started_at   TEXT,
    finished_at  TEXT
);
CREATE INDEX IF NOT EXISTS idx_jobs_study ON jobs(study_id);
CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);

CREATE TABLE IF NOT EXISTS results (
    job_id       TEXT PRIMARY KEY,
    study_id     TEXT NOT NULL,
    stage        TEXT,
    parsed       TEXT,
    validation   TEXT,
    created_at   TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_results_study ON results(study_id);

CREATE TABLE IF NOT EXISTS events (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    ts           TEXT NOT NULL,
    study_id     TEXT,
    job_id       TEXT,
    level        TEXT,
    message      TEXT
);
CREATE INDEX IF NOT EXISTS idx_events_study ON events(study_id);
"""


def now() -> str:
    """Local ISO timestamp (seconds)."""
    return datetime.now().astimezone().isoformat(timespec="seconds")


def new_id(prefix: str) -> str:
    return f"{prefix}{datetime.now():%Y%m%d%H%M%S}_{secrets.token_hex(3)}"


def _dumps(obj) -> str | None:
    if obj is None:
        return None
    return json.dumps(obj, ensure_ascii=False, default=str)


def _loads(txt: str | None):
    if not txt:
        return None
    try:
        return json.loads(txt)
    except (TypeError, ValueError):
        return txt


class Store:
    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.db_path), timeout=60)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.executescript(_SCHEMA)
        # migration: add the lane column on older databases
        cols = {r[1] for r in self.conn.execute("PRAGMA table_info(jobs)")}
        if "lane" not in cols:
            self.conn.execute("ALTER TABLE jobs ADD COLUMN lane INTEGER NOT NULL DEFAULT 0")
        self.conn.commit()

    # -- low level ------------------------------------------------------
    def close(self) -> None:
        try:
            self.conn.close()
        except sqlite3.Error:
            pass

    def _exec(self, sql: str, params=()):
        cur = self.conn.execute(sql, params)
        self.conn.commit()
        return cur

    @staticmethod
    def _row(row) -> dict | None:
        if row is None:
            return None
        d = dict(row)
        for k in ("spec", "meta", "variant", "parsed", "validation"):
            if k in d:
                d[k] = _loads(d[k])
        return d

    # -- studies --------------------------------------------------------
    def create_study(self, **kw) -> dict:
        sid = kw.get("id") or new_id("st_")
        ts = now()
        self._exec(
            """INSERT INTO studies
               (id,batch_id,name,smiles,variant,charge,multiplicity,route,workflow,
                solvent,preopt,status,spec,meta,created_at,updated_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                sid, kw.get("batch_id"), kw["name"], kw.get("smiles"),
                _dumps(kw.get("variant")), kw.get("charge"), kw.get("multiplicity"),
                kw.get("route"), kw.get("workflow"), kw.get("solvent"),
                int(bool(kw.get("preopt", True))), kw.get("status", "prepared"),
                _dumps(kw.get("spec")), _dumps(kw.get("meta")), ts, ts,
            ),
        )
        return self.get_study(sid)

    def get_study(self, study_id: str) -> dict | None:
        row = self.conn.execute("SELECT * FROM studies WHERE id=?", (study_id,)).fetchone()
        return self._row(row)

    def update_study(self, study_id: str, **fields) -> None:
        allowed = {
            "batch_id", "name", "smiles", "variant", "charge", "multiplicity",
            "route", "workflow", "solvent", "preopt", "status", "spec", "meta",
        }
        sets, vals = [], []
        for k, v in fields.items():
            if k not in allowed:
                continue
            sets.append(f"{k}=?")
            if k in ("variant", "spec", "meta"):
                vals.append(_dumps(v))
            elif k == "preopt":
                vals.append(int(bool(v)))
            else:
                vals.append(v)
        if not sets:
            return
        sets.append("updated_at=?")
        vals.extend([now(), study_id])
        self._exec(f"UPDATE studies SET {', '.join(sets)} WHERE id=?", vals)

    def list_studies(self, batch_id: str | None = None, status: str | None = None,
                     limit: int = 200) -> list[dict]:
        sql, params = "SELECT * FROM studies WHERE 1=1", []
        if batch_id:
            sql += " AND batch_id=?"
            params.append(batch_id)
        if status:
            sql += " AND status=?"
            params.append(status)
        sql += " ORDER BY created_at LIMIT ?"
        params.append(limit)
        return [self._row(r) for r in self.conn.execute(sql, params)]

    # -- jobs -----------------------------------------------------------
    def create_job(self, **kw) -> dict:
        jid = kw.get("id") or new_id("jb_")
        ts = now()
        self._exec(
            """INSERT INTO jobs
               (id,study_id,stage,seq,attempt,lane,status,dir,route,charge,multiplicity,
                input_path,output_path,chk_path,pid,error_class,parent_job,retry_of,
                meta,created_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                jid, kw["study_id"], kw["stage"], kw.get("seq", 0),
                kw.get("attempt", 1), kw.get("lane", 0),
                kw.get("status", "pending"), kw.get("dir"),
                kw.get("route"), kw.get("charge"), kw.get("multiplicity"),
                kw.get("input_path"), kw.get("output_path"), kw.get("chk_path"),
                kw.get("pid"), kw.get("error_class"), kw.get("parent_job"),
                kw.get("retry_of"), _dumps(kw.get("meta")), ts,
            ),
        )
        return self.get_job(jid)

    def get_job(self, job_id: str) -> dict | None:
        row = self.conn.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
        return self._row(row)

    def update_job(self, job_id: str, **fields) -> None:
        allowed = {
            "status", "dir", "route", "charge", "multiplicity", "input_path",
            "output_path", "chk_path", "pid", "exit_code", "timed_out",
            "walltime", "error_class", "parent_job", "retry_of", "meta",
            "started_at", "finished_at", "attempt", "seq", "stage",
        }
        sets, vals = [], []
        for k, v in fields.items():
            if k not in allowed:
                continue
            sets.append(f"{k}=?")
            vals.append(_dumps(v) if k == "meta" else v)
        if not sets:
            return
        vals.append(job_id)
        self._exec(f"UPDATE jobs SET {', '.join(sets)} WHERE id=?", vals)

    def list_jobs(self, study_id: str | None = None, status: str | None = None,
                  limit: int = 500) -> list[dict]:
        sql, params = "SELECT * FROM jobs WHERE 1=1", []
        if study_id:
            sql += " AND study_id=?"
            params.append(study_id)
        if status:
            sql += " AND status=?"
            params.append(status)
        sql += " ORDER BY created_at, seq LIMIT ?"
        params.append(limit)
        return [self._row(r) for r in self.conn.execute(sql, params)]

    def study_jobs(self, study_id: str) -> list[dict]:
        rows = self.conn.execute(
            "SELECT * FROM jobs WHERE study_id=? ORDER BY seq, created_at", (study_id,)
        ).fetchall()
        return [self._row(r) for r in rows]

    # -- queue ----------------------------------------------------------
    def claim_next_job(self) -> dict | None:
        """Atomically claim the oldest pending job, keeping one study's
        stages contiguous (lane ordering)."""
        self.conn.execute("BEGIN IMMEDIATE")
        try:
            row = self.conn.execute(
                "SELECT * FROM jobs WHERE status='pending' "
                "ORDER BY lane, created_at, seq LIMIT 1"
            ).fetchone()
            if row is None:
                self.conn.commit()
                return None
            self.conn.execute(
                "UPDATE jobs SET status='running', started_at=? WHERE id=? AND status='pending'",
                (now(), row["id"]),
            )
            self.conn.commit()
        except Exception:
            self.conn.rollback()
            raise
        job = self.get_job(row["id"])
        return job

    def next_lane(self) -> int:
        row = self.conn.execute("SELECT COALESCE(MAX(lane), 0) + 1 AS n FROM jobs").fetchone()
        return int(row["n"]) if row else 1

    def running_jobs(self) -> list[dict]:
        return [self._row(r) for r in self.conn.execute(
            "SELECT * FROM jobs WHERE status='running'")]

    def pending_count(self) -> int:
        row = self.conn.execute(
            "SELECT COUNT(*) AS n FROM jobs WHERE status IN ('pending','running')").fetchone()
        return int(row["n"]) if row else 0

    # -- results --------------------------------------------------------
    def save_result(self, job_id: str, study_id: str, stage: str,
                    parsed: dict, validation: dict) -> None:
        self._exec(
            """INSERT OR REPLACE INTO results (job_id,study_id,stage,parsed,validation,created_at)
               VALUES (?,?,?,?,?,?)""",
            (job_id, study_id, stage, _dumps(parsed), _dumps(validation), now()),
        )

    def get_result(self, job_id: str) -> dict | None:
        row = self.conn.execute("SELECT * FROM results WHERE job_id=?", (job_id,)).fetchone()
        if row is None:
            return None
        d = dict(row)
        d["parsed"] = _loads(d.pop("parsed"))
        d["validation"] = _loads(d.pop("validation"))
        return d

    def study_results(self, study_id: str) -> dict[str, dict]:
        rows = self.conn.execute(
            "SELECT * FROM results WHERE study_id=?", (study_id,)).fetchall()
        out = {}
        for r in rows:
            d = dict(r)
            d["parsed"] = _loads(d.pop("parsed"))
            d["validation"] = _loads(d.pop("validation"))
            out[d["job_id"]] = d
        return out

    # -- events ---------------------------------------------------------
    def add_event(self, message: str, study_id: str | None = None,
                  job_id: str | None = None, level: str = "info") -> None:
        self._exec(
            "INSERT INTO events (ts,study_id,job_id,level,message) VALUES (?,?,?,?,?)",
            (now(), study_id, job_id, level, message),
        )

    def list_events(self, study_id: str | None = None, limit: int = 200) -> list[dict]:
        sql, params = "SELECT * FROM events WHERE 1=1", []
        if study_id:
            sql += " AND study_id=?"
            params.append(study_id)
        sql += " ORDER BY id DESC LIMIT ?"
        params.append(limit)
        return [dict(r) for r in self.conn.execute(sql, params)]
