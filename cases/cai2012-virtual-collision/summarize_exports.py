#!/usr/bin/env python3
"""Refresh the index/overview from completed per-case exports; no model rerun."""
from pathlib import Path
import sys
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from plumetools.virtualcollision.export_summary import summarize
if __name__ == '__main__':
    print(f'Summarized {summarize(ROOT/"cases"/"cai2012-kn-extension")} completed cases.')
