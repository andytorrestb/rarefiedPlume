#!/usr/bin/env python3
"""Generate separate100-cell local patches; never launch or modify existing runs."""
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from plumetools.radialstudy.__main__ import main

if __name__=='__main__':
    raise SystemExit(main(['generate','--study',str(Path(__file__).with_name('study-local.yaml')),
                          '--output',str(Path(__file__).with_name('Cases-local')),*sys.argv[1:]]))
