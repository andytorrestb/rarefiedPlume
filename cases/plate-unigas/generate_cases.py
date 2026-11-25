#!/usr/bin/env python3
"""Generate the plate-impingement Kn sweep on uniGasFoam (cases/plate-unigas).

    python generate_cases.py --dry-run         the 21 cases and their estimates; writes nothing
    python generate_cases.py                   Cases/<case>/ and manifest-plate.yaml
    python generate_cases.py --only Kn100_zD25 one case (no manifest)
    python generate_cases.py --smoke           Smoke/Kn100_zD25, the laptop smoke run
    python generate_cases.py --no-estimate     skip the ~30 s-per-case nozzleplume estimate

Nothing here meshes, runs a solver or submits a job.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))

from plumetools.plateugf import generation  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--no-estimate", action="store_true")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    try:
        if args.smoke:
            print(f"wrote {generation.stage_smoke(HERE)}")
            return 0
        entries = generation.generate(HERE, dry_run=args.dry_run, only=args.only,
                                      estimate=not args.no_estimate, workers=args.workers)
    except generation.GenerationError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1
    print("\n".join(generation.table(entries)))
    print("\n(dry run: nothing written)" if args.dry_run else f"\nwrote {len(entries)} case(s) under {HERE / 'Cases'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
