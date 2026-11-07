"""Locate Gaussian 16 / GaussView 6 / Shermo on an arbitrary Windows machine.

The packaged .exe must work on other computers, where these programs may live
in very different places, so detection combines:

1. explicit paths already stored in the config,
2. environment variables (GAUSS_EXEDIR, PATH),
3. well-known installation directories,
4. registry uninstall entries,
5. a bounded keyword scan of all fixed drives,
6. running processes.

All results are scored; the best candidate is returned together with the full
candidate list so the UI can offer a dropdown / manual override.
"""
from __future__ import annotations

import os
import re
import string
from pathlib import Path

_KEYWORDS = re.compile(
    r"(gauss|g16|g09|gview|gaussview|shermo|gaussian16|gaussian\s*16)", re.I)

_COMMON_DIRS = [
    r"C:\G16W", r"C:\G16w", r"C:\Gaussian16", r"C:\Gaussian 16",
    r"C:\Program Files\Gaussian 16", r"C:\Program Files\Gaussian",
    r"C:\Program Files (x86)\Gaussian 16", r"C:\Program Files (x86)\Gaussian",
    r"D:\G16W", r"D:\g16w", r"D:\Gaussian16", r"D:\Gaussian 16",
    r"D:\Program Files\Gaussian 16", r"E:\G16W", r"E:\Gaussian16",
    r"C:\GaussView 6", r"C:\GaussView6", r"C:\Program Files\GaussView 6",
    r"C:\Program Files (x86)\GaussView 6", r"D:\GaussView 6", r"D:\GaussView6",
    r"C:\Shermo", r"D:\Shermo", r"C:\Program Files\Shermo",
]

_SKIP_DIRS = {
    "windows", "winsxs", "system32", "syswow64", "temp", "tmp",
    "$recycle.bin", "system volume information", "node_modules", ".git",
    "appdata", "recovery", "perflogs", "programdata",
}

# binaries that identify a Gaussian installation
_GAUSS_MARKERS = ("g16.exe", "g16w.exe", "l9999.exe", "formchk.exe")
_GVIEW_MARKERS = ("gview.exe",)


def _fixed_drives() -> list[Path]:
    out = []
    for letter in string.ascii_uppercase:
        root = Path(f"{letter}:\\")
        try:
            if root.exists():
                out.append(root)
        except OSError:
            continue
    return out


def _dir_has(path: Path, markers) -> bool:
    try:
        names = {p.name.lower() for p in path.iterdir() if p.is_file()}
    except OSError:
        return False
    return all(m.lower() in names for m in markers)


def _find_exe(path: Path, names) -> Path | None:
    for name in names:
        p = path / name
        if p.is_file():
            return p
    return None


def gaussian_version(home: Path) -> str | None:
    """Infer the Gaussian revision from the release notes filenames."""
    try:
        for p in home.iterdir():
            m = re.match(r"release_a(\d+)\.txt", p.name, re.I)
            if m:
                return f"16 Rev A.{m.group(1)}"
    except OSError:
        pass
    return None


def _candidate(path: Path, kind: str) -> dict | None:
    path = Path(path)
    if not path.exists():
        return None
    home = path if path.is_dir() else path.parent
    exe = None
    if kind == "gaussian":
        exe = _find_exe(home, ("g16.exe", "g16w.exe"))
        if exe is None:
            sub = home / "g16w"
            if sub.is_dir():
                exe = _find_exe(sub, ("g16.exe",))
                if exe:
                    home = sub
    elif kind == "gaussview":
        exe = _find_exe(home, ("gview.exe",))
    elif kind == "shermo":
        exe = _find_exe(home, ("shermo.exe", "Shermo.exe"))
    if exe is None:
        return None
    info = {"kind": kind, "home": str(home), "exe": str(exe)}
    if kind == "gaussian":
        info["version"] = gaussian_version(home)
        info["has_gview"] = (home / "gview.exe").is_file()
        try:
            info["n_link_exes"] = len([p for p in home.glob("l[0-9]*.exe")])
        except OSError:
            info["n_link_exes"] = 0
    if kind == "gaussview":
        info["version"] = _gview_version(home)
    return info


def _gview_version(home: Path) -> str | None:
    for name in ("fullversionname.txt", "GaussView.url"):
        p = home / name
        if p.is_file():
            try:
                txt = p.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            m = re.search(r"(\d+\.\d+\.\d+)", txt)
            if m:
                return m.group(1)
    return None


def _registry_candidates() -> list[Path]:
    paths: list[Path] = []
    try:
        import winreg
    except ImportError:  # pragma: no cover - non Windows
        return paths
    keys = [
        (winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall"),
        (winreg.HKEY_LOCAL_MACHINE,
         r"SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall"),
        (winreg.HKEY_CURRENT_USER, r"SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall"),
    ]
    for hive, sub in keys:
        try:
            root = winreg.OpenKey(hive, sub)
        except OSError:
            continue
        for i in range(256):
            try:
                name = winreg.EnumKey(root, i)
            except OSError:
                break
            try:
                with winreg.OpenKey(root, name) as k:
                    try:
                        display = winreg.QueryValueEx(k, "DisplayName")[0]
                    except OSError:
                        continue
                    if not _KEYWORDS.search(str(display)):
                        continue
                    for val in ("InstallLocation", "DisplayIcon", "UninstallString"):
                        try:
                            raw = str(winreg.QueryValueEx(k, val)[0])
                        except OSError:
                            continue
                        raw = raw.strip('"').split(",")[0]
                        p = Path(raw)
                        if p.suffix.lower() == ".exe":
                            p = p.parent
                        if p.exists():
                            paths.append(p)
            except OSError:
                continue
    # Gaussian also writes a plain SOFTWARE\Gaussian key on some installs
    for hive in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
        for sub in (r"SOFTWARE\Gaussian", r"SOFTWARE\Gaussian 16"):
            try:
                with winreg.OpenKey(hive, sub) as k:
                    for i in range(64):
                        try:
                            n, v, _ = winreg.EnumValue(k, i)
                        except OSError:
                            break
                        if isinstance(v, str) and v:
                            p = Path(v.split(",")[0].strip('"'))
                            if p.suffix.lower() == ".exe":
                                p = p.parent
                            if p.exists() and ("gauss" in n.lower() or p.name.lower().startswith("g")):
                                paths.append(p)
            except OSError:
                continue
    return paths


def _scan_dirs(roots: list[Path], max_depth: int = 3, max_visited: int = 4000) -> list[Path]:
    """Bounded breadth-first scan for directories whose name matches keywords."""
    hits: list[Path] = []
    seen = set()
    queue = [(r, 0) for r in roots if r.exists()]
    visited = 0
    while queue and visited < max_visited:
        current, depth = queue.pop(0)
        try:
            entries = list(os.scandir(current))
        except OSError:
            continue
        for e in entries:
            try:
                if not e.is_dir(follow_symlinks=False):
                    continue
            except OSError:
                continue
            name = e.name
            if name.lower() in _SKIP_DIRS:
                continue
            p = Path(e.path)
            key = str(p).lower()
            if key in seen:
                continue
            seen.add(key)
            visited += 1
            if _KEYWORDS.search(name):
                hits.append(p)
            # descend into keyword dirs, plus shallow levels and standard roots
            if depth < max_depth and (depth < 2 or _KEYWORDS.search(name)):
                queue.append((p, depth + 1))
    return hits


def find_software(kind: str, explicit: str | None = None) -> dict:
    """Return {'best': info|None, 'candidates': [info...]} for a program kind."""
    candidates: list[dict] = []
    sources: dict[str, str] = {}

    def add(path, source):
        if not path:
            return
        info = _candidate(Path(path), kind)
        if info is None:
            return
        key = info["exe"].lower()
        if key in sources:
            return
        sources[key] = source
        info["source"] = source
        candidates.append(info)

    # 1) explicit config path
    if explicit:
        add(Path(explicit), "config")
        add(Path(explicit).parent if not Path(explicit).is_dir() else Path(explicit), "config")

    # 2) environment
    if kind == "gaussian":
        exedir = os.environ.get("GAUSS_EXEDIR")
        if exedir:
            add(Path(exedir), "GAUSS_EXEDIR")
    for d in os.environ.get("PATH", "").split(os.pathsep):
        try:
            p = Path(d)
            if p.is_dir() and (p / "g16.exe").is_file():
                add(p, "PATH")
            elif p.is_dir() and (p / "gview.exe").is_file():
                add(p, "PATH")
            elif p.is_dir() and (p / "shermo.exe").is_file():
                add(p, "PATH")
        except OSError:
            continue

    # 3) well-known dirs
    for d in _COMMON_DIRS:
        add(Path(d), "common-dir")

    # 4) registry
    for p in _registry_candidates():
        add(p, "registry")

    # 5) bounded drive scan
    roots: list[Path] = []
    for drive in _fixed_drives():
        roots.append(drive)
        for sub in ("Program Files", "Program Files (x86)", "Software", "Apps",
                    "Programs", "tools", "Tools", "chem", "Chem"):
            roots.append(drive / sub)
    for p in _scan_dirs(roots):
        add(p, "scan")

    # 6) running processes
    try:
        import psutil
        wanted = {"gaussian": {"g16.exe", "g16w.exe"}, "gaussview": {"gview.exe"},
                  "shermo": {"shermo.exe"}}[kind]
        for proc in psutil.process_iter(["name", "exe"]):
            try:
                if (proc.info.get("name") or "").lower() in wanted and proc.info.get("exe"):
                    add(Path(proc.info["exe"]).parent, "process")
            except (psutil.Error, OSError):
                continue
    except Exception:
        pass

    # rank: explicit/env/common/registry > scan/process; prefer gaussview alongside Gaussian
    rank = {"config": 0, "GAUSS_EXEDIR": 1, "PATH": 1, "common-dir": 2,
            "registry": 3, "process": 4, "scan": 5}
    candidates.sort(key=lambda c: (rank.get(c.get("source"), 9),
                                   -int(c.get("n_link_exes") or 0)))
    return {"best": candidates[0] if candidates else None, "candidates": candidates}


def detect_all(explicit: dict | None = None) -> dict:
    """Detect all supported programs; returns a JSON-friendly dict."""
    explicit = explicit or {}
    g = find_software("gaussian", explicit.get("g16_exe"))
    v = find_software("gaussview", explicit.get("gview_exe"))
    s = find_software("shermo", explicit.get("shermo_exe"))

    # If GaussView wasn't found standalone but lives in the Gaussian folder
    if v["best"] is None and g["best"]:
        home = Path(g["best"]["home"])
        cand = _candidate(home, "gaussview")
        if cand:
            cand["source"] = "gaussian-home"
            v = {"best": cand, "candidates": [cand]}

    return {
        "gaussian": g,
        "gaussview": v,
        "shermo": s,
        "found": {
            "gaussian": g["best"] is not None,
            "gaussview": v["best"] is not None,
            "shermo": s["best"] is not None,
        },
    }
