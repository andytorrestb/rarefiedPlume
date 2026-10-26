#!/usr/bin/env python3
"""Run the virtual model for existing per-case DSMC/Cai CSV comparisons."""
import argparse
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from plumetools.virtualcollision import NumericalOptions
from plumetools.virtualcollision.exported import run_exports


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', type=Path, default=ROOT/'cases'/'cai2012-kn-extension')
    p.add_argument('--samples', type=int, default=2048)
    p.add_argument('--hazard-steps', type=int, default=24)
    p.add_argument('--grid-refinement', type=int, choices=(0, 1, 2), default=1)
    p.add_argument('--kn', nargs='+', type=float,
                   help='Only these Kn values; extension cases are 10, 1, 0.001 and 0.0001')
    p.add_argument('--no-plots', action='store_true')
    args = p.parse_args()
    options = NumericalOptions(samples=args.samples, batches=8, hazard_steps=args.hazard_steps,
        partner_radial=4, partner_azimuthal=16, partner_speed=16,
        domain_length_D=20, domain_radius_D=16, chunk_size=128)
    processed, skipped = run_exports(ROOT, args.study, options,
        refinement=args.grid_refinement, plots=not args.no_plots, selected_kn=args.kn)
    if processed and not args.no_plots:
        from plumetools.virtualcollision.export_summary import summarize
        summarize(args.study)
    print(f'Completed {len(processed)} cases; {len(skipped)} missing/unsupported records.')
    return 0 if processed else 1


if __name__ == '__main__':
    raise SystemExit(main())
