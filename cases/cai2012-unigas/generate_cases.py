#!/usr/bin/env python3
"""Generate the cai2012-unigas cases idempotently, without running anything.

    ./generate_cases.py                 create what is missing, keep what matches
    ./generate_cases.py --dry-run       the plan; write nothing
    ./generate_cases.py --refresh       copy updated runner scripts into cases
                                        whose inputs are unchanged
    ./generate_cases.py --regenerate    replace CONFLICTING cases that hold no
                                        solver output
    ./generate_cases.py --overwrite     replace conflicting cases EVEN IF THEY
                                        HOLD RESULTS -- which destroys them

The cases are the kn-extension's, case for case (study.yaml and
baseCase/case.yaml are checked against cases/cai2012-kn-extension before
anything is staged); unigas.yaml carries the solver settings. Both are inputs:
changing either makes a generated case conflict. See
plumetools/cai2012ugf/generation.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012.config import CaiConfigError
    from plumetools.cai2012.study import StudyError
    from plumetools.cai2012ugf import generation
    from plumetools.cai2012ugf.settings import UniGasConfigError
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')

LABELS = {"create": "create", "keep": "keep (identical)",
          "keep-stale": "keep (STALE scripts/meta; --refresh)",
          "refresh": "refresh scripts/meta", "replace": "REPLACE",
          "refuse": "REFUSED -- conflict"}


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generate the Cai Knudsen sweep on uniGasFoam.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--pilot", metavar="KN_NAME",
                        help="stage only Pilot/<case> for this Kn (e.g. Kn0p001), "
                             "outside the manifest (spec §8)")
    return parser.parse_args(argv[1:])


def main(argv) -> int:
    args = parse_args(argv)
    if args.pilot:
        try:
            path = generation.stage_pilot(HERE, args.pilot)
        except (generation.GenerationError, StudyError, CaiConfigError, UniGasConfigError) as exc:
            print(f"generate_cases: {exc}", file=sys.stderr)
            return 1
        print(path)
        return 0
    try:
        states = generation.generate(HERE, dry_run=args.dry_run, refresh=args.refresh,
                                     regenerate=args.regenerate, overwrite=args.overwrite)
    except (generation.GenerationError, StudyError, CaiConfigError, UniGasConfigError) as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1
    print(f"{'case':<20} {'Kn':>8}  action")
    for state in states:
        print(f"{state.case.name:<20} {state.case.knudsen:>8g}  {LABELS[state.action]}"
              + ("  [holds solver output]" if state.has_results else ""))
        for line in (state.inputs if state.action in ("refuse", "replace") else
                     state.auxiliary if state.action in ("keep-stale", "refresh") else [])[:8]:
            print(f"{'':<22}{line}")
        if getattr(state, "known_limitations", ""):
            print(f"{'':<22}known source limitation: direct-exit inflow surface "
                  f"(KNOWN-LIMITATIONS)")
    if any(s.blocking for s in states):
        print("\ngenerate_cases: conflicting cases were NOT changed and the manifest "
              "was not written. Re-run with --regenerate (no results) after "
              "checking the differences above.", file=sys.stderr)
        return 1
    if args.dry_run:
        print("\n--dry-run: nothing written.")
        return 0
    print(f"\nmanifest  {HERE / 'manifest-unigas.yaml'}")
    print("Every case is post-run-verdict-pending until it has run (spec §7).")
    print("\nNext:  python validate_cases.py   then   qsub runMeshAll.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
