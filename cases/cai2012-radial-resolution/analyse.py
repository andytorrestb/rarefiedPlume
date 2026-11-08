#!/usr/bin/env python3
"""Inventory the study and report measured gaps and available comparisons."""
from pathlib import Path
import sys

for root in Path(__file__).resolve().parents:
    if (root/'plumetools').is_dir():
        sys.path.insert(0,str(root))
        break
from plumetools.radialstudy.__main__ import main

if __name__ == '__main__':
    raise SystemExit(main(['analyse',str(Path(__file__).with_name('Cases-adaptive')),*sys.argv[1:]]))
