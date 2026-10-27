#!/usr/bin/env python3
"""Inventory, fit and report empirical whole-field Cai coordinate maps."""
from pathlib import Path
import os
import sys

# BLAS reads these when NumPy is imported, before the configuration is loaded.
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from plumetools.cai2012ext.scaling_workflow import main

if __name__ == "__main__":
    raise SystemExit(main(study_root=Path(__file__).resolve().parent))
