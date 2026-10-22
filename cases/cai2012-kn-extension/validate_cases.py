#!/usr/bin/env python3
"""Check the generated Knudsen-extension tree before committing HPC hours to it.

    ./validate_cases.py                 every check, dictionaries included
    ./validate_cases.py --inputs-only   only what needs no generated dictionaries
    ./validate_cases.py --verbose       also print what passed

Exit status is 0 when the tree is exactly the extension and 1 when it is not, so
it drops straight into a job script:

    ./validate_cases.py && sbatch ...

What it checks
--------------
The extension's own promises, each RECOMPUTED rather than read back
(plumetools/cai2012ext/verification.py):

    exactly Kn 1e-4, 1e-3, 1, 10 at the baseline parcel level -- four cases
    no Kn = 0.01, 0.1 or 100, in the manifest or anywhere under Cases/
    n0 reproduces the target Kn (VHS lambda written out again)
    D, S0, T0, the gas, the baseline occupancy, startup and averaging identical
    nEquivalentParticles x 40 = n0 V_exit
    sampling starts at 10^4 Kn = 0.01 collision times in every case
    averaging lasts 4.5 domain transits (to within the snapped last write)
    every case.yaml is still the one the manifest describes

and then cases/cai2012's own tree check (plumetools.cai2012.verify), which
reads the GENERATED DICTIONARIES -- dsmcProperties, controlDict, blockMeshDict --
and checks each case's weight, deltaT, endTime, write schedule and fieldAverage
block against the manifest. That half needs

    ./AllmeshCases --dict-only        writes them; needs no OpenFOAM

Production suitability is NOT a validation failure: it is a recorded verdict.
This script prints the count, and ./AllrunCases enforces it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012.verify import VerificationError, load_manifest
    from plumetools.cai2012ext.verification import verify_extension
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    parser = argparse.ArgumentParser(
        description="Check the Knudsen-extension tree against the study.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="print every check, not only the failures")
    parser.add_argument("--inputs-only", action="store_true",
                        help="skip the checks that read generated dictionaries")
    parser.add_argument("--manifest", default=str(HERE / "manifest.yaml"),
                        help="path to manifest.yaml (default: alongside this "
                             "script)")
    args = parser.parse_args(argv[1:])

    manifest_path = Path(args.manifest)
    try:
        manifest = load_manifest(manifest_path)
    except VerificationError as exc:
        print(f"validate_cases: {exc}", file=sys.stderr)
        return 1

    root = manifest_path.resolve().parent
    cases = manifest.get("cases") or []
    counts = manifest.get("production_status_counts") or {}
    print(f"manifest   {manifest_path}")
    print(f"cases      {len(cases)} ({manifest.get('n_physical_cases', '?')} "
          f"Kn x {manifest.get('n_particle_levels', '?')} particle levels)")
    print(f"excluded   Kn {manifest.get('excluded_knudsen')}")
    print(f"verdicts   {counts}")
    print(f"checks     {'inputs only' if args.inputs_only else 'inputs + generated dictionaries'}")
    print()

    report = verify_extension(root, manifest,
                              with_dictionaries=not args.inputs_only)
    for line in report.lines(verbose=args.verbose):
        print(line)

    print()
    if report.passed:
        print(f"validate_cases: {len(cases)} case(s) OK, "
              f"{len(report.passes)} checks passed.")
        unsuitable = [c["name"] for c in cases
                      if c.get("production_status") != "production-candidate"]
        if unsuitable:
            print(f"                {len(unsuitable)} are UNSUITABLE FOR "
                  f"PRODUCTION and will not be run by ./AllrunCases:")
            print("                " + ", ".join(unsuitable))
        print()
        print("Next:")
        print("    ./AllrunCases --list         what would run, for a job array")
        print("    ./AllrunCases                run the production candidates")
        return 0

    print(f"validate_cases: {len(report.failures)} FAILURE(S).", file=sys.stderr)
    print("Do not submit this tree. Regenerate it:", file=sys.stderr)
    print("    ./generate_cases.py --regenerate && ./AllmeshCases --dict-only",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
