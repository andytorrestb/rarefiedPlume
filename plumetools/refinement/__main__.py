"""Standalone ``python -m plumetools.refinement`` commands for measured cases."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .assess import assess
from .config import load_settings
from .foam import read_snapshot
from .plan import plan_refinement, write_plan, read_plan
from .report import write_assessment


def main(argv=None) -> int:
    """Assess/export without case changes; mesh application is an explicit command."""
    parser = argparse.ArgumentParser(description="Measured DSMC cell sizing and portable refinement plans")
    subs = parser.add_subparsers(dest="command", required=True)
    for name in ("assess", "plan"):
        p = subs.add_parser(name)
        p.add_argument("case", type=Path)
        p.add_argument("--settings", type=Path, required=True)
        p.add_argument("--time")
        p.add_argument("--output", type=Path, required=(name == "plan"))
        p.add_argument("--dry-run", action="store_true")
    apply = subs.add_parser("apply")
    apply.add_argument("case", type=Path)
    apply.add_argument("--plan", type=Path, required=True)
    apply.add_argument("--backend", choices=["hex"], default="hex")
    apply.add_argument("--output", type=Path, required=True)
    apply.add_argument("--max-stages", type=int, default=20)
    apply.add_argument("--max-cells", type=int)
    args = parser.parse_args(argv)
    try:
        if args.command == "apply":
            from .openfoam import HexRefiner
            result = HexRefiner(args.max_stages, args.max_cells).refine(args.case, read_plan(args.plan), args.output)
            print(f"{result.target_status}: {result.n_cells} cells; "
                  + ("new flow run required" if result.requires_flow_rerun else "no mesh correction needed"))
            if result.destination:
                print(f"wrote mesh artifact {result.destination}")
            return 0
        cfg = load_settings(args.settings)
        snap = read_snapshot(args.case, cfg, time=args.time)
        result = assess(snap.geometry, snap.measured, snap.species, convention=snap.convention, targets=cfg.targets)
        print(f"{result.summary['action']}: {result.summary['cells']} cells; "
              f"{result.summary['mesh_failed_cells']} mesh failures; "
              f"{result.summary['unresolved_cells']} unresolved; target edge/MFP <= {cfg.targets.edge_over_mfp:g}")
        if not args.dry_run and args.output is not None:
            if args.command == "plan":
                plan = plan_refinement(snap.geometry, snap.measured, result, particle_weight=snap.particle_weight)
                path = write_plan(plan, args.output)
            else:
                path = write_assessment(snap.geometry, snap.measured, result, args.output)
            print(f"wrote {path}")
        return 0
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"refinement: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
