"""Optional Shermo thermochemistry integration.

Shermo (Tian Lu) post-processes a Gaussian frequency output.  It is NOT
bundled with this application (it is a separate free download); when the exe
is configured / found we run it non-interactively and parse the standard
``========== Total ==========`` block:

    Shermo.exe  freq.out  -T 298.15  -noset

Everything goes to stdout, so we redirect it to ``shermo_thermo.txt`` next to
the Gaussian output and store the parsed G/H/U/S values in the result JSON.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

_CACHE: dict = {}

_NUM = r"[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"


def _num(x) -> float:
    return float(str(x).replace("D", "E").replace("d", "e"))


def shermo_path(cfg) -> str | None:
    """Locate Shermo: config first, then the regular auto-detection scan."""
    configured = cfg.get("shermo.exe")
    if configured and Path(configured).is_file():
        return str(configured)
    key = str(cfg.home)
    if key in _CACHE:
        return _CACHE[key]
    exe = None
    try:
        from .locate import find_software
        res = find_software("shermo", configured)
        if res.get("best"):
            exe = res["best"]["exe"]
    except Exception:
        exe = None
    _CACHE[key] = exe
    return exe


def parse_shermo_output(text: str) -> dict:
    """Parse the Shermo stdout summary into structured values."""
    out: dict = {"ok": True, "U_hartree": None, "H_hartree": None,
                 "G_hartree": None, "S_cal_mol_k": None,
                 "zpe_energy_hartree": None, "electronic_energy_hartree": None,
                 "raw": []}
    m = re.search(r"Sum of electronic energy and ZPE[^:]*:\s*(" + _NUM + r")\s*a\.u\.", text)
    if m:
        out["zpe_energy_hartree"] = _num(m.group(1))
    for key, label in (("U_hartree", "U"), ("H_hartree", "H"), ("G_hartree", "G")):
        m = re.search(r"Sum of electronic energy and thermal correction to "
                      + label + r":\s*(" + _NUM + r")\s*a\.u\.", text)
        if m:
            out[key] = _num(m.group(1))
    m = re.search(r"Total S:\s+" + _NUM + r"\s*J/mol/K\s+(" + _NUM + r")\s*cal/mol/K", text)
    if m:
        out["S_cal_mol_k"] = _num(m.group(1))
    m = re.search(r"Electronic energy:\s*(" + _NUM + r")", text)
    if m:
        out["electronic_energy_hartree"] = _num(m.group(1))
    out["raw"] = [ln.strip() for ln in text.splitlines()
                  if any(k in ln for k in ("Total S:", "Sum of electronic energy",
                                           "Electronic energy", "ZPE"))][:12]
    if out["G_hartree"] is None and out["H_hartree"] is None:
        out["ok"] = False
        out["error"] = "未能从 Shermo 输出中解析热化学量"
    return out


def run_shermo(cfg, output_path: str | Path,
               energy_hartree: float | None = None,
               temperature: float | None = None) -> dict:
    """Run Shermo on a Gaussian frequency output. Never raises."""
    out_path = Path(output_path)
    if not out_path.exists():
        return {"ok": False, "error": f"输出文件不存在: {out_path}"}
    exe = shermo_path(cfg)
    if not exe:
        return {"ok": False, "error": "未找到 Shermo.exe（可在设置中手动指定）",
                "available": False}
    T = float(temperature or cfg.get("redox.temperature_k") or 298.15)
    args = [exe, str(out_path), "-T", f"{T:g}", "-noset"]
    if energy_hartree is not None:
        args += ["-E", f"{float(energy_hartree):.10f}"]
    timeout = int(cfg.get("shermo.timeout_seconds") or 300)
    try:
        proc = subprocess.run(args, cwd=str(out_path.parent), capture_output=True,
                              text=True, encoding="utf-8", errors="replace",
                              timeout=timeout,
                              creationflags=0x08000000 if __import__("os").name == "nt" else 0)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"Shermo 超时 ({timeout}s)", "exe": exe}
    except OSError as exc:
        return {"ok": False, "error": f"Shermo 启动失败: {exc}", "exe": exe}
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    log_file = out_path.parent / "shermo_thermo.txt"
    try:
        log_file.write_text(text, encoding="utf-8", errors="replace")
    except OSError:
        pass
    parsed = parse_shermo_output(text)
    parsed.update({"exe": exe, "returncode": proc.returncode,
                   "temperature_k": T, "stdout_file": str(log_file)})
    if proc.returncode != 0 and not parsed.get("ok"):
        parsed["error"] = f"Shermo 退出码 {proc.returncode}"
    return parsed
