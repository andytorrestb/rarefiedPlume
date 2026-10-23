#!/usr/bin/env python3
"""Rarefaction and mesh-resolution diagnostics -- the nozzle-plume family's.

    ./Alldiagnose.py                  the latest time
    ./Alldiagnose.py --time 0.063     a specific time directory
    ./Alldiagnose.py --no-fields      the summary only; write nothing into <time>/

Run after ./Allpost, which writes the cell centres and volumes. What R-4D's
./Allpost step 4 does (python -m plumetools.nozzleplume diagnostics), fed from
this cai2012 case.yaml and from ./diagnostics.yaml:

    <time>/mfp, KnD, KnCellMax, KnCellV, cellOverMfpMax, cellOverMfpV,
           cellSizeMax, cellSizeV, cellSizeTarget, KnGradN, KnGradT,
           mfpResolved, parcelsPerStep, sampledParcels, diagValid,
           diagValidGradN, diagValidGradT, diagStatus
                               VOLUME FIELDS beside rhoNMean -- open open.foam
                               in ParaView. -1 where there is no value
    results/diagnostics_summary.yaml
                               validity, Delta/lambda statistics, the worst
                               cells, the lambda/k mesh estimate, and this
                               family's exit_plane check
    results/centreline_diagnostics.csv, centreline-diagnostics.png
                               the same along the axis, on centerline.csv's
                               stations

The images of those fields come from ../viz-diagnostics.yaml, which
../AllpostCases renders. Every field is defined in docs/rarefaction-
diagnostics.md; plumetools/cai2012ext/diagnostics.py is the adapter.

Needs the reconstructed mesh (constant/polyMesh) and the time-averaged fields;
no OpenFOAM.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012ext import diagnostics, policy
    from plumetools.rarefaction.foam import PostError
    from plumetools.rarefaction.settings import SettingsError
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--time", default=None,
                        help="time directory to read (default: the latest)")
    parser.add_argument("--no-fields", action="store_true",
                        help="do not write the diagnostic fields into <time>/")
    args = parser.parse_args(argv[1:])

    print(f"### {HERE.name}: rarefaction and mesh-resolution diagnostics")
    try:
        fields, written, _ = diagnostics.run(HERE, time=args.time,
                                             write_fields=not args.no_fields)
    except policy.ExtensionError as exc:
        print(f"Alldiagnose: {exc}", file=sys.stderr)
        return 2
    except (SettingsError, PostError, ValueError) as exc:
        print(f"Alldiagnose: {exc}", file=sys.stderr)
        return 1

    if fields:
        print(f"    wrote {len(fields)} fields into {Path(fields[0]).parent}")
    for path in written:
        print(f"    wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
