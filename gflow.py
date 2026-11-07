#!/usr/bin/env python
"""Entry point: ``python gflow.py <command>`` (or via gflow.cmd)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gaussflow.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
