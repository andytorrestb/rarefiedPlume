#!/usr/bin/env python3
"""Generate the Cai 2012 case hierarchy from study.yaml.

    ./generate_cases.py                    generate every enabled case
    ./generate_cases.py --dry-run          show what would be generated
    ./generate_cases.py --knudsen 100 0.1  only these Knudsen numbers
    ./generate_cases.py --overwrite        replace cases that already have results

No arguments are needed for normal use.

What it does:

  1. copy baseCase into Cases/<name>
  2. rewrite each case-local case.yaml with that case's Knudsen number --
     STRUCTURALLY, by loading and re-emitting YAML, never by string
     substitution
  3. derive every quantity the case will run at, and record all of them in
     manifest.yaml

The OpenFOAM dictionaries are NOT written here. Each case writes its own from
its own case.yaml when ./Allmesh runs, so there is one authoritative source of
physical inputs per case and no chance of a stale dictionary disagreeing with
it. The derivation is still run here, because a case whose mesh cannot be built
or whose checks fail should be found now rather than an hour into ./AllmeshCases.

Existing cases that have results are refused unless --overwrite is given:
regenerating rewrites case.yaml, and doing that under finished output would
leave results whose inputs no longer describe them.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012 import checks, inflow, mesh
    from plumetools.cai2012.config import CaiConfigError, load_case_config
    from plumetools.cai2012.geometry import from_config
    from plumetools.cai2012.study import (
        StudyError,
        apply_case_parameters,
        case_is_complete,
        clone_cases,
        filter_by_knudsen,
        load_study,
        manifest_entry,
        write_manifest,
    )
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generate the Cai 2012 case hierarchy from study.yaml.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="list the cases that would be generated and stop")
    parser.add_argument("--knudsen", nargs="+", type=float, metavar="KN",
                        help="only these Knudsen numbers (default: every "
                             "enabled case)")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace cases that already contain results")
    parser.add_argument("--study", default=str(HERE / "study.yaml"),
                        help="path to study.yaml (default: alongside this script)")
    return parser.parse_args(argv[1:])


def derive(case_dir: Path):
    """Everything the case will run at, from its own rewritten case.yaml.

    Returns ``(cfg, geom, exit_state, plan, run, nozzle, particles, report)``.
    Loading the case-local file rather than the template is the point: it is the
    only way to be sure the manifest records what the case will actually do.
    """
    cfg = load_case_config(case_dir)
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    particles = mesh.estimate_particles(
        plan, geom, exit_state, run.n_equivalent_particles)
    report = checks.run_checks(cfg, geom, exit_state, plan, run,
                               nozzle=nozzle, particles=particles)
    return cfg, geom, exit_state, plan, run, nozzle, particles, report


def main(argv) -> int:
    args = parse_args(argv)

    try:
        study = load_study(Path(args.study))
    except StudyError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    root = Path(args.study).resolve().parent
    base = root / study.base_case

    try:
        selected = filter_by_knudsen(study.enabled_cases(), args.knudsen)
    except StudyError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    disabled = [c for c in study.all_cases() if not c.enabled]

    print(f"study      {args.study}")
    print(f"base case  {base}")
    print(f"reference  {study.meta.get('reference', 'Cai and Wang 2012')}")
    print()

    if not selected:
        print("Nothing to generate: every case is disabled, or the --knudsen "
              "filter matched none.")
        return 1

    print(f"{'case':<10} {'Kn':>8}  destination")
    for case in selected:
        print(f"{case.name:<10} {case.knudsen:>8g}  {study.case_path(case)}")
    for case in disabled:
        print(f"{case.name:<10} {case.knudsen:>8g}  SKIPPED (enabled: false)")
    print()

    if args.dry_run:
        print("--dry-run: nothing written.")
        return 0

    # Refuse the whole run rather than half of it: generating two cases and then
    # stopping on the third leaves a tree nobody can reason about.
    blocked = [c for c in selected
               if case_is_complete(root / study.case_path(c)) and not args.overwrite]
    if blocked:
        print("generate_cases: these cases already contain results:",
              file=sys.stderr)
        for case in blocked:
            print(f"    {study.case_path(case)}", file=sys.stderr)
        print("\nRegenerating rewrites case.yaml, which would leave results "
              "whose\ninputs no longer describe them. Re-run with --overwrite to "
              "replace\nthem, or with --knudsen to select the others.",
              file=sys.stderr)
        return 1

    for case in selected:
        destination = root / study.case_path(case)
        if destination.exists() and args.overwrite:
            shutil.rmtree(destination)

    try:
        clone_cases(root, study, selected)
    except StudyError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    entries = []
    failed = []
    for case in selected:
        destination = root / study.case_path(case)
        try:
            apply_case_parameters(destination, case, study)
        except StudyError as exc:
            print(f"generate_cases: {exc}", file=sys.stderr)
            return 1

        try:
            cfg, geom, exit_state, plan, run, nozzle, particles, report = derive(
                destination)
        except (CaiConfigError, ValueError) as exc:
            print(f"generate_cases: {case.name}: {exc}", file=sys.stderr)
            failed.append(case.name)
            continue

        entries.append(manifest_entry(case, study, cfg, geom, exit_state, plan,
                                      run, nozzle=nozzle, particles=particles))

        print(f"generated {case.name:<10} Kn {case.knudsen:>7g}")
        print(f"          lambda0 {exit_state.mean_free_path_m:.4e} m   "
              f"n0 {exit_state.number_density_per_m3:.4e} 1/m^3   "
              f"U0 {exit_state.velocity_m_per_s:.1f} m/s")
        print(f"          mesh {plan.n_cells:,} cells, core "
              f"{plan.core_cell_size_m:.4e} m, cell/lambda0 "
              f"{plan.cell_over_mfp:.3g}"
              + ("  (COARSENED to fit mesh.max_cells)" if plan.coarsened else ""))
        print(f"          nozzle {nozzle.n_faces} faces, area "
              f"{100 * nozzle.area_error:+.2f}% from pi R0^2")
        print(f"          weight {run.n_equivalent_particles:.4e}, deltaT "
              f"{run.delta_t_s:.4e} s, {run.n_steps:.0f} steps, "
              f"~{particles['total_particles']:.2e} parcels")
        for result in report.results:
            if result.status != "pass":
                print(f"          {result.status.upper()}: {result.name}: "
                      f"{result.message}")
        if report.failures:
            failed.append(case.name)

    manifest = write_manifest(root / study.manifest_name, study, entries)
    print()
    print(f"manifest  {manifest}")

    if failed:
        print()
        print(f"generate_cases: {len(failed)} case(s) have failing DSMC checks: "
              f"{failed}", file=sys.stderr)
        print("They were generated, and the manifest records them, but ./Allmesh "
              "will\nrefuse to write their dictionaries until case.yaml is fixed.",
              file=sys.stderr)
        return 1

    print()
    print("Next:")
    print("    ./AllmeshCases          mesh every generated case")
    print("    ./AllrunCases           run them")
    print("    ./AllpostCases          post-process and build the study table")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
