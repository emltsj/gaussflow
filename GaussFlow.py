#!/usr/bin/env python
"""GaussFlow desktop entry point (source or PyInstaller-frozen).

    python GaussFlow.py              # open the GUI
    python GaussFlow.py --selftest out.json
"""
import sys
from pathlib import Path

if not getattr(sys, "frozen", False):
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from gaussflow.gui import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
