"""``python -m plumetools.hybridplume <step> <case> [options]`` -- see workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import workflow


def _print(doc) -> None:
    print(json.dumps(doc, indent=2, default=workflow._json_default)[:20000])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m plumetools.hybridplume")
    sub = ap.add_subparsers(dest="step", required=True)

    s = sub.add_parser("setup", help="resolve profile, assess, write the CFD case")
    s.add_argument("case")
    s.add_argument("--config", required=True)
    s.add_argument("--profile")
    s.add_argument("--np", type=int)
    s.add_argument("--background-ratio", type=float,
                   help="write a background-floor sensitivity copy instead of cfd/")

    for name in ("check-mesh", "cfd-report", "breakdown", "interface"):
        p = sub.add_parser(name)
        p.add_argument("case")
        if name in ("check-mesh",):
            p.add_argument("--subdir", default="cfd")
        if name in ("breakdown", "interface", "cfd-report"):
            p.add_argument("--time", help="CFD time directory (default: latest)")
        if name == "interface":
            p.add_argument("--allow-unsuitable", action="store_true",
                           help="build the DSMC region even if the surface fails the "
                                "criteria; the outcome is recorded as a smoke test only")

    args = ap.parse_args(argv)
    try:
        if args.step == "setup":
            if args.background_ratio is not None:
                doc = workflow.setup(Path(args.case), Path(args.config), args.profile,
                                     n_ranks=args.np, background_ratio=args.background_ratio,
                                     subdir=workflow.background_label(args.background_ratio))
            else:
                doc = workflow.setup(Path(args.case), Path(args.config), args.profile,
                                     n_ranks=args.np)
            _print({k: doc[k] for k in ("label", "profile", "mesh", "source_discretisation",
                                        "resources", "n_ranks")})
        elif args.step == "check-mesh":
            _print(workflow.check_mesh(Path(args.case), args.subdir))
        else:
            fn = getattr(workflow, args.step.replace("-", "_"))
            kwargs = {k: v for k, v in vars(args).items() if k not in ("step", "case")}
            _print(fn(Path(args.case), **kwargs))
    except workflow.WorkflowError as exc:
        print(f"hybridplume {args.step}: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
