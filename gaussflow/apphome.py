"""Application home / data directory bootstrap (works for source and frozen exe).

* Running from source: data stays inside the project folder.
* Running as a packaged .exe: data goes to ``%LOCALAPPDATA%\\GaussFlow`` so the
  exe can live in Program Files (read-only) or on a USB stick.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def app_home(portable: bool = False) -> Path:
    if getattr(sys, "frozen", False):
        if portable:
            return Path(sys.executable).resolve().parent
        base = os.environ.get("LOCALAPPDATA") or str(Path.home())
        return Path(base) / "GaussFlow"
    return Path(__file__).resolve().parents[1]


def bootstrap_home(portable: bool = False) -> Path:
    home = app_home(portable=portable)
    home.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("GFLOW_HOME", str(home))
    return home


def resource_path(*parts: str) -> Path:
    """Path to a bundled resource (PyInstaller ``_MEIPASS`` aware)."""
    base = getattr(sys, "_MEIPASS", None)
    if base:
        return Path(base).joinpath(*parts)
    return Path(__file__).resolve().parents[1].joinpath(*parts)
