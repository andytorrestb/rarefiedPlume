#!/usr/bin/env python3
"""Generate the prediction-based radial cases; use --dry-run to inspect costs."""
from pathlib import Path
import sys

for root in Path(__file__).resolve().parents:
    if (root/'plumetools').is_dir():
        sys.path.insert(0,str(root))
        break
from plumetools.radialstudy.__main__ import main

if __name__ == '__main__':
    raise SystemExit(main(['generate','--study',str(Path(__file__).with_name('study.yaml')),
                          '--output',str(Path(__file__).with_name('Cases')),*sys.argv[1:]]))
