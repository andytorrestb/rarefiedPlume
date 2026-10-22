#!/usr/bin/env python3
"""Overlap validation: this truncated run against the full-domain benchmark.

    ./Alloverlap.py                         baseline from interface.yaml
    ./Alloverlap.py --baseline <case>       another copy of it (read-only)

Interface-VALIDATION cases only. Reads the predeclared region, targets and
uncertainty settings from ``interface.yaml`` (``overlap:``), the boundary frames
of both runs one at a time, and writes:

    results/overlap.yaml        E_n, E_speed, E_T with intervals and verdicts,
                                u_R, u_theta, excluded and low-signal bins, the
                                definitions and the region
    results/overlap_bins.csv    every bin: both runs, errors, standard errors

The baseline is never written to. See plumetools/cai2012ext/overlap.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012ext import overlap
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


if __name__ == "__main__":
    raise SystemExit(overlap.main([str(HERE), *sys.argv[1:]]))
