"""Gaussian job execution with resource limits and a wall-clock watchdog."""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

_CREATE_NO_WINDOW = 0x08000000
_CREATE_NEW_PROCESS_GROUP = 0x00000200

# job_id -> pid of the currently running g16.exe (for forced cancellation)
_RUNNING_PIDS: dict[str, int] = {}
# job ids for which the user pressed "stop"; run_job reports outcome.cancelled
_CANCEL_REQUESTED: set[str] = set()


def running_pids() -> dict[str, int]:
    """Snapshot of job_id -> g16 pid for all currently running jobs."""
    return dict(_RUNNING_PIDS)


def is_cancel_requested(job_id: str) -> bool:
    return str(job_id) in _CANCEL_REQUESTED


def kill_pid(pid: int) -> None:
    _kill_tree(int(pid))


def request_cancel(job_id: str) -> int | None:
    """Mark a job as cancelled and kill its Gaussian process tree immediately."""
    job_id = str(job_id)
    _CANCEL_REQUESTED.add(job_id)
    pid = _RUNNING_PIDS.get(job_id)
    if pid:
        _kill_tree(int(pid))
        return int(pid)
    return None


def _kill_tree(pid: int) -> None:
    if os.name == "nt":
        try:
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(pid)],
                capture_output=True, timeout=30,
                creationflags=_CREATE_NO_WINDOW,
            )
            return
        except Exception:
            pass
    try:
        os.kill(pid, 9)
    except Exception:
        pass


def _pid_alive(pid: int) -> bool:
    if not pid:
        return False
    try:
        import psutil
        return psutil.pid_exists(pid)
    except Exception:
        try:
            os.kill(pid, 0)
            return True
        except Exception:
            return False


def gaussian_env(cfg) -> dict:
    env = os.environ.copy()
    if cfg.g16_home:
        env["GAUSS_EXEDIR"] = str(cfg.g16_home)
        env["PATH"] = str(cfg.g16_home) + os.pathsep + env.get("PATH", "")
    env["GAUSS_SCRDIR"] = ""  # replaced per job
    env.pop("GAUSS_SCRDIR", None)
    return env


def run_job(cfg, job: dict, on_start=None) -> dict:
    """Run one Gaussian job synchronously (called from the worker)."""
    if not cfg.g16_exe:
        raise RuntimeError(
            "Gaussian g16.exe is not configured. Run `python gflow.py detect` "
            "or set g16_exe in config.json / the Settings tab.")
    job_dir = Path(job["dir"])
    job_dir.mkdir(parents=True, exist_ok=True)
    scratch = job_dir / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    input_path = Path(job["input_path"])
    output_path = Path(job["output_path"])

    env = gaussian_env(cfg)
    # NOTE (verified on this machine): Gaussian 16W ignores GAUSS_SCRDIR and
    # always writes scratch into the *current directory*.  We therefore run
    # with cwd=job_dir and sweep Gau-* leftovers afterwards.
    env["GAUSS_SCRDIR"] = str(scratch)

    console_path = job_dir / "console.log"
    max_wall = float(cfg.limits.get("max_walltime_seconds", 7200))
    t0 = time.time()
    rc = None
    timed_out = False
    job_id = str(job.get("id") or "")
    cancelled = False

    with open(console_path, "w", encoding="utf-8", errors="replace") as log:
        flags = _CREATE_NO_WINDOW
        # The input MUST be passed as a bare relative filename in cwd:
        # absolute paths break the Windows driver (truncation at spaces).
        proc = subprocess.Popen(
            [str(cfg.g16_exe), str(input_path.name)],
            cwd=str(job_dir),
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            creationflags=flags if os.name == "nt" else 0,
        )
        if job_id:
            _RUNNING_PIDS[job_id] = proc.pid
        if on_start:
            try:
                on_start(proc.pid)
            except Exception:
                pass
        try:
            while True:
                try:
                    rc = proc.wait(timeout=2.0)
                    break
                except subprocess.TimeoutExpired:
                    if job_id and job_id in _CANCEL_REQUESTED:
                        cancelled = True
                        _kill_tree(proc.pid)
                        try:
                            rc = proc.wait(timeout=30)
                        except subprocess.TimeoutExpired:
                            rc = -9
                        break
                    if time.time() - t0 > max_wall:
                        timed_out = True
                        _kill_tree(proc.pid)
                        try:
                            rc = proc.wait(timeout=30)
                        except subprocess.TimeoutExpired:
                            rc = -9
                        break
        except KeyboardInterrupt:  # pragma: no cover
            _kill_tree(proc.pid)
            rc = -2
            raise
        finally:
            walltime = time.time() - t0
            if job_id:
                _RUNNING_PIDS.pop(job_id, None)
            if job_id and job_id in _CANCEL_REQUESTED:
                cancelled = True
                _CANCEL_REQUESTED.discard(job_id)

    # Gaussian 16W writes <input_stem>.out next to the input
    if not output_path.exists():
        candidate = job_dir / (input_path.stem + ".out")
        if candidate.exists():
            output_path = candidate

    # sweep scratch leftovers (Gaussian deletes them itself on success only)
    leftover_bytes = 0
    for p in list(job_dir.glob("Gau-*")):
        try:
            leftover_bytes += p.stat().st_size
            p.unlink()
        except OSError:
            pass
    try:
        if scratch.exists() and not any(scratch.iterdir()):
            scratch.rmdir()
    except OSError:
        pass

    return {
        "exit_code": rc,
        "timed_out": timed_out,
        "cancelled": cancelled,
        "walltime": round(walltime, 2),
        "output_path": str(output_path),
        "output_exists": output_path.exists(),
        "console_log": str(console_path),
        "scratch_leftover_bytes": leftover_bytes,
    }
