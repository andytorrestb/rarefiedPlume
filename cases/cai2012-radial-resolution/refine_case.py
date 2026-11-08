#!/usr/bin/env python3
"""Create a fresh revision only when selected measured criteria require it."""
from pathlib import Path
import sys

for root in Path(__file__).resolve().parents:
    if (root/'plumetools').is_dir():
        sys.path.insert(0,str(root))
        break
from plumetools.radialstudy.__main__ import main

if __name__ == '__main__':
    raise SystemExit(main(['refine',*sys.argv[1:]]))
