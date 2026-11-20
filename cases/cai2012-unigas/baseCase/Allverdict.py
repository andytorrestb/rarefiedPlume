#!/usr/bin/env python3
"""Write this case's post-run suitability verdict (spec §7).

    ./Allverdict.py               results/unigas_verdict.yaml from the latest time
    ./Allverdict.py --time T      from time T

Judges the DSMC cells (uniGasCollisionModelId = 1) of the reconstructed,
converted fields: effective collision spacing over the mean free path, the time
step over the mean collision time, and particles per cell above the occupancy
floor; then Allconverge.py's stationarity verdict and the known source
limitations. Adds the pre-run gates (from the manifest) and this verdict to
results/provenance-unigas.yaml. See plumetools/cai2012ugf/verdict.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    import yaml

    from plumetools.cai2012ugf import hpc, verdict
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def prerun_gates() -> dict | None:
    """This case's pre-run gates from the family manifest, if it is in it."""
    manifest = HERE.parent.parent / "manifest-unigas.yaml"
    if not manifest.is_file():
        return None
    document = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    for entry in document.get("cases") or []:
        if entry.get("name") == HERE.name:
            return {"prerun_gates": entry.get("prerun_gates"),
                    "dsmcfoam_prerun_status": entry.get("dsmcfoam_prerun_status")}
    return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--time", default=None)
    args = parser.parse_args(argv)
    try:
        result = verdict.assess_run(HERE, time=args.time)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Allverdict: {exc}", file=sys.stderr)
        return 1
    hpc.add_verdicts(HERE, prerun_gates(),
                     {key: result[key] for key in ("status", "production_ready", "unresolved", "time")})
    print(f"  {result['status'].upper()} at t = {result['time']}")
    print(f"  DSMC cells {result['dsmc']['n_cells']:,} ({result['dsmc']['fraction']:.1%}), "
          f"USP cells {result['usp']['n_cells']:,} ({result['usp']['fraction']:.1%})")
    for reason in result["unresolved"]:
        print(f"    - {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
