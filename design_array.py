#!/usr/bin/env python3
"""Convenience entry point: ``python design_array.py --help``."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from thz_opt.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
