"""Background worker daemon: claim -> run -> parse -> validate -> advance.

Long Gaussian jobs must be asynchronous: the Agent submits, gets job ids and
polls ``gflow status``; it never blocks on a job.  The daemon is a single
detached Python process (one job at a time, so the resource budget stays
predictable).
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

from .config import Config
from .pipeline import after_job
from .runner import _pid_alive, run_job
from .store import Store, now

_CREATE_NO_WINDOW = 0x08000000
_CREATE_NEW_PROCESS_GROUP = 0x00000200
_DETACHED = 0x00000008


# ---------------------------------------------------------------------------
def pid_file(cfg: Config) -> Path:
    return cfg.home / "worker.pid"


def stop_file(cfg: Config) -> Path:
    return cfg.home / "worker.stop"


def daemon_status(cfg: Config) -> dict:
    pf = pid_file(cfg)
    if not pf.exists():
        return {"running": False, "pid": None}
    try:
        pid = int(pf.read_text(encoding="utf-8").strip())
    except (OSError, ValueError):
        return {"running": False, "pid": None, "stale": True}
    alive = _pid_alive(pid)
    return {"running": alive, "pid": pid, "stale": not alive}


def start_daemon(cfg: Config) -> dict:
    status = daemon_status(cfg)
    if status.get("running"):
        return {"ok": True, "already_running": True, **status}
    pf = pid_file(cfg)
    if pf.exists():
        pf.unlink()
    sf = stop_file(cfg)
    if sf.exists():
        sf.unlink()

    log_path = cfg.log_dir / "worker.log"
    log = open(log_path, "a", encoding="utf-8", errors="replace")
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root) + os.pathsep + env.get("PYTHONPATH", "")
    env["GFLOW_CONFIG"] = str(cfg.path)
    flags = 0
    if os.name == "nt":
        flags = _DETACHED | _CREATE_NEW_PROCESS_GROUP | _CREATE_NO_WINDOW
    proc = subprocess.Popen(
        [sys.executable, "-m", "gaussflow.daemon"],
        cwd=str(root),
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=log,
        stderr=subprocess.STDOUT,
        creationflags=flags,
        close_fds=True,
    )
    log.close()
    # wait until the daemon publishes its pid file
    deadline = time.time() + 10
    while time.time() < deadline:
        st = daemon_status(cfg)
        if st.get("running"):
            return {"ok": True, "already_running": False, **st, "log": str(log_path)}
        if proc.poll() is not None:
            return {"ok": False, "error": "daemon 启动后立即退出",
                    "log": str(log_path), "exit_code": proc.returncode}
        time.sleep(0.2)
    return {"ok": False, "error": "daemon 启动超时（见 worker.log）",
            "log": str(log_path)}


def stop_daemon(cfg: Config, timeout: float = 30.0) -> dict:
    status = daemon_status(cfg)
    if not status.get("running"):
        pf = pid_file(cfg)
        if pf.exists():
            pf.unlink()
        return {"ok": True, "was_running": False}
    stop_file(cfg).write_text("stop", encoding="utf-8")
    pid = status["pid"]
    deadline = time.time() + timeout
    while time.time() < deadline:
        if not _pid_alive(pid):
            break
        time.sleep(0.3)
    if _pid_alive(pid):
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(pid)],
                           capture_output=True, creationflags=_CREATE_NO_WINDOW)
        else:
            try:
                os.kill(pid, 9)
            except OSError:
                pass
    for f in (pid_file(cfg), stop_file(cfg)):
        if f.exists():
            try:
                f.unlink()
            except OSError:
                pass
    return {"ok": True, "was_running": True, "pid": pid}


# ---------------------------------------------------------------------------
class Worker:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.store = Store(cfg.db_path)
        self.running = True
        self.stop = stop_file(cfg)

    # -- recovery -------------------------------------------------------
    def recover(self) -> None:
        """Requeue/mark jobs left 'running' by a crashed daemon."""
        for job in self.store.running_jobs():
            pid = job.get("pid")
            if pid and _pid_alive(int(pid)):
                continue
            output = job.get("output_path")
            outcome = {"exit_code": None, "timed_out": 0, "walltime": None,
                       "output_path": output, "output_exists": bool(
                           output and Path(output).exists())}
            self.store.update_job(job["id"], status="failed", pid=None,
                                  error_class="interrupted")
            self.store.add_event(
                "worker recovered interrupted job", study_id=job["study_id"],
                job_id=job["id"], level="warning")
            try:
                after_job(self.store, self.cfg, job, outcome)
            except Exception:
                self.store.add_event(
                    "recovery failed: " + traceback.format_exc()[:800],
                    study_id=job["study_id"], job_id=job["id"], level="error")

    # -- single job -----------------------------------------------------
    def process_one(self) -> bool:
        job = self.store.claim_next_job()
        if job is None:
            return False
        self.store.update_job(job["id"], status="running")
        study = self.store.get_study(job["study_id"])
        if study and study["status"] == "submitted":
            self.store.update_study(job["study_id"], status="running")
        self.store.add_event(f"job {job['stage']} started", study_id=job["study_id"],
                             job_id=job["id"])
        # let the previous Gaussian process fully release files/AV locks
        delay = float(self.cfg.limits.get("inter_job_delay_seconds", 1.0))
        if delay > 0:
            time.sleep(delay)
        try:
            outcome = run_job(self.cfg, job,
                              on_start=lambda pid: self.store.update_job(
                                  job["id"], pid=pid))
            self.store.update_job(job["id"], pid=None,
                                  walltime=outcome.get("walltime"))
            if outcome.get("cancelled"):
                self.store.update_job(
                    job["id"], status="cancelled",
                    error_class="stopped_by_user", finished_at=now())
                self.store.update_study(job["study_id"], status="cancelled")
                self.store.add_event(
                    f"job {job['stage']} cancelled by user (Gaussian terminated)",
                    study_id=job["study_id"], job_id=job["id"], level="warning")
                return True
            result = after_job(self.store, self.cfg,
                               self.store.get_job(job["id"]), outcome)
            self.store.add_event(
                f"job finished -> {result.get('action')}",
                study_id=job["study_id"], job_id=job["id"])
        except Exception:
            self.store.update_job(job["id"], status="failed", pid=None,
                                  error_class="worker_exception")
            self.store.add_event("worker exception: " + traceback.format_exc()[:1200],
                                 study_id=job["study_id"], job_id=job["id"],
                                 level="error")
            self.store.update_study(job["study_id"], status="failed")
        return True

    # -- loop -----------------------------------------------------------
    def loop(self) -> None:
        self.recover()
        poll = float(self.cfg.limits.get("poll_seconds", 2.0))
        while self.running:
            if self.stop.exists():
                break
            try:
                if self.process_one():
                    continue
            except Exception:
                self.store.add_event("worker loop error: " + traceback.format_exc()[:1200],
                                     level="error")
            time.sleep(poll)

    def shutdown(self) -> None:
        self.running = False


def main() -> int:
    cfg = Config.load()
    pf = pid_file(cfg)
    pf.write_text(str(os.getpid()), encoding="utf-8")
    sf = stop_file(cfg)
    if sf.exists():
        sf.unlink()
    worker = Worker(cfg)
    try:
        worker.loop()
    finally:
        try:
            pf.unlink()
        except OSError:
            pass
        worker.store.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
