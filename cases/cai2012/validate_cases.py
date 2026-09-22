#!/usr/bin/env python3
"""Check the generated Cai 2012 tree before committing HPC hours to it.

    ./validate_cases.py               check every case in manifest.yaml
    ./validate_cases.py --verbose     also print what passed
    ./validate_cases.py --help

Exit status is 0 when the tree matches the study and 1 when it does not, so this
drops straight into a job script:

    ./validate_cases.py && sbatch ...

What it checks
--------------
That the tree on disk IS the study study.yaml describes, and in particular that
each physical Cai case's four numerical-particle variants differ in the particle
weight and in nothing else:

    np1x    nEquivalentParticles == the value this case has always run at
    np2x    == np1x / 2
    np5x    == np1x / 5
    np10x   == np1x / 10

    deltaT, blockMeshDict, and every physical key in case.yaml -- the gas, the
    species, the geometry, the boundary conditions, the collision model, the
    density, the transient -- identical across the four

    endTime  == run_time_multiplier x the baseline
    writes   == output_frequency_multiplier x more often, under the SAME
                writeControl the study has always used
    fieldAverage present, averaging from the same timeStart in all four, with
    prime2Mean on for rhoN and dsmcRhoN

It reads the GENERATED OPENFOAM DICTIONARIES, not the Python that wrote them --
re-deriving the values and comparing them with themselves would pass whatever
the directories happened to contain. So the dictionaries have to exist:

    ./AllmeshCases --dict-only     writes them; needs no OpenFOAM
    ./AllmeshCases                 writes them and meshes

There is no fixing here. A failure means regenerate:

    ./generate_cases.py --overwrite && ./AllmeshCases
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012.verify import (
        VerificationError,
        load_manifest,
        verify_tree,
    )
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Check the generated Cai 2012 tree against study.yaml.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="print every check, not only the failures")
    parser.add_argument("--manifest", default=str(HERE / "manifest.yaml"),
                        help="path to manifest.yaml (default: alongside this "
                             "script)")
    return parser.parse_args(argv[1:])


def main(argv) -> int:
    args = parse_args(argv)
    manifest_path = Path(args.manifest)

    try:
        manifest = load_manifest(manifest_path)
    except VerificationError as exc:
        print(f"validate_cases: {exc}", file=sys.stderr)
        return 1

    root = manifest_path.resolve().parent
    cases = manifest.get("cases") or []

    print(f"manifest   {manifest_path}")
    print(f"cases      {len(cases)} "
          f"({manifest.get('n_physical_cases', '?')} physical x "
          f"{manifest.get('n_particle_levels', '?')} particle level(s))")
    print(f"schedule   run time x{manifest.get('run_time_multiplier', 1)}, "
          f"output frequency x{manifest.get('output_frequency_multiplier', 1)}")
    print()

    report = verify_tree(root, manifest)
    for line in report.lines(verbose=args.verbose):
        print(line)

    print()
    if report.passed:
        print(f"validate_cases: {len(cases)} case(s) OK, "
              f"{len(report.passes)} checks passed.")
        print()
        print("Next:")
        print("    ./AllrunCases           run them")
        print("    ./AllrunCases --list    the case list, for a job array")
        return 0

    print(f"validate_cases: {len(report.failures)} FAILURE(S).", file=sys.stderr)
    print("Do not submit this tree. Regenerate it:", file=sys.stderr)
    print("    ./generate_cases.py --overwrite && ./AllmeshCases",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
