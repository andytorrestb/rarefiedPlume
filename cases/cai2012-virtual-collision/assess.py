#!/usr/bin/env python3
"""Read-only three-way free-plume assessment; never generates or runs DSMC."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

for _thread_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_var] = str(min(8, max(1, int(os.environ.get(_thread_var, "1")))))

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from plumetools.virtualcollision.assessment import inventory, load_yaml, write_json, write_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("inventory", "smoke", "full"))
    parser.add_argument("--root", type=Path, default=ROOT, help="repository containing source case families")
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("assessment.yaml"))
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results"))
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    cfg = load_yaml(args.config)
    cfg["_root"] = str(args.root)
    document = inventory(args.root, families=tuple(cfg.get("families", ("cai2012", "cai2012-kn-extension"))))
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "inventory.json", document)
    write_rows(args.output / "inventory.csv", document["cases"])
    skips = [{"path": row["path"], "kn": row["nominal_kn"], "quality": row["quality"],
              "reasons": row["reasons"]} for row in document["cases"] if not row["primary_included"]]
    write_json(args.output / "skipped.json", skips)
    if args.mode != "inventory":
        from plumetools.virtualcollision.demonstration import run_model_demonstration
        report = run_model_demonstration(cfg, args.output, smoke=args.mode == "smoke", plots=not args.no_plots)
        from plumetools.virtualcollision.assessment import run_reference_assessment
        assessment = run_reference_assessment(cfg, document, args.output,
                                             smoke=args.mode == "smoke", plots=not args.no_plots)
        report["source_revision"] = document["source_revision"]
        report["reference_assessment"] = assessment
        report["scientific_effectiveness"] = ("exploratory only; assess quality flags before interpreting errors"
            if assessment["processed_cases"] else "not assessed: no compatible reference processed")
        write_json(args.output / "report.json", report)
    print(f"Discovered {len(document['cases'])} records; {len(skips)} explicitly excluded/skipped.")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
