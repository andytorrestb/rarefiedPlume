#!/usr/bin/env python3
"""Inventory, fit and report the empirical Cai scaling for the uniGasFoam cases.

    python Allscale.py inventory          results/scaling/inventory.{json,csv}
    python Allscale.py fit --full         fits, with the uniGasFoam stratum apart
    python Allscale.py report             README and figures from assessment.json

The kn-extension's workflow functions, on this family's inventory: every row
carries solver uniGasFoam, a stratum of its own, and the post-run verdict as its
suitability (plumetools/cai2012ugf/scaling.py). Nothing here pools solvers.
"""
from pathlib import Path
import os
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from plumetools.cai2012ugf.scaling import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(study_root=Path(__file__).resolve().parent))
