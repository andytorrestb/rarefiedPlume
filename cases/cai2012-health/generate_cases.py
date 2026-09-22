#!/usr/bin/env python3
"""Generate the statistical-health matrix from study.yaml.

    ./generate_cases.py                     generate the whole matrix
    ./generate_cases.py --dry-run           what it would generate, and the cost
    ./generate_cases.py --weight ppc005     only these weight rows
    ./generate_cases.py --sampling s0p5     only these sampling columns
    ./generate_cases.py --overwrite         replace cases that already have results

No arguments are needed for normal use.

What it does:

  1. clone baseCase into Cases/<weight>/<sampling> with CaseFoam
  2. rewrite each case-local case.yaml with that cell's TWO swept values --
     STRUCTURALLY, by loading and re-emitting YAML, never by string
     substitution -- and nothing else
  3. derive every quantity each case will run at, and record all of them in
     manifest.yaml, including the ones that must be identical across the matrix
  4. cost the matrix, from measurements if any exist and from a prior if not

COST IT BEFORE RUNNING IT. Run time goes as parcels x total transits, so the
top weight row is roughly half the bill however the matrix is shaped, and
finding that out from a progress bar is finding it out too late. --dry-run
prints the estimate per case and the total, and says whether it is a prior or a
fit. cases/cai2012-health/AllrunCases records each case's wall clock as it
finishes, so the estimate sharpens as the matrix runs.

The OpenFOAM dictionaries are NOT written here. Each case writes its own from
its own case.yaml when ./Allmesh runs, so there is one authoritative source of
physical inputs per case and no chance of a stale dictionary disagreeing with
it. The derivation is still run here, because a case whose mesh cannot be built
or whose checks fail should be found now rather than an hour into ./AllmeshCases.

Existing cases that have results are refused unless --overwrite is given.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012 import checks, health, inflow, mesh
    from plumetools.cai2012.config import CaiConfigError, load_case_config
    from plumetools.cai2012.geometry import from_config
    from plumetools.cai2012.study import StudyError, case_is_complete
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test,cases]"')


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generate the Cai 2012 statistical-health matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="list the cases and their predicted cost, then stop")
    parser.add_argument("--weight", nargs="+", metavar="NAME",
                        help="only these weight rows (default: every enabled row)")
    parser.add_argument("--sampling", nargs="+", metavar="NAME",
                        help="only these sampling columns")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace cases that already contain results")
    parser.add_argument("--study", default=str(HERE / "study.yaml"),
                        help="path to study.yaml (default: alongside this script)")
    return parser.parse_args(argv[1:])


def select(study, weights, samplings):
    """The cells matching the command-line filters.

    A filter naming a level that does not exist stops the run rather than
    quietly generating nothing -- the same rule ``--knudsen`` follows in
    ``cases/cai2012``.
    """
    cases = study.enabled_cases()
    for wanted, available, axis in (
            (weights, {w.name for w in study.enabled_weights()}, "weight"),
            (samplings, {s.name for s in study.enabled_samplings()}, "sampling")):
        if not wanted:
            continue
        missing = sorted(set(wanted) - available)
        if missing:
            raise StudyError(
                f"no {axis} level {missing}; the study defines "
                f"{sorted(available)}")
    if weights:
        cases = [c for c in cases if c.weight.name in set(weights)]
    if samplings:
        cases = [c for c in cases if c.sampling.name in set(samplings)]
    return cases


def derive(case_dir: Path):
    """Everything the case will run at, from its own rewritten case.yaml.

    Loading the case-local file rather than the template is the point: it is
    the only way to be sure the manifest records what the case will actually
    do -- including the values that are supposed to be the same everywhere.
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


#: Reconstructed bytes per written time directory, for the disk estimate.
#:
#: MEASURED on cases/cai2012/Cases/Kn100: a sampling-window write is 455 MB, of
#: which ~255 MB is the 1.31e6-cell field set and ~200 MB the parcel cloud at
#: 1.2e6 parcels. The field part does not scale with the particle weight and the
#: cloud does, so they are carried separately -- an estimate that scaled the
#: whole write with the weight would be four times wrong on the cheapest row.
#:
#: AllrunCases prunes the reconstructed cloud from every time but the latest,
#: which is why the default estimate counts fields only.
FIELD_BYTES_PER_CELL = 255.0e6 / 1.31e6
CLOUD_BYTES_PER_PARCEL = 200.0e6 / 1.2e6


def estimate(case, plan, run, particles, model, *, keep_clouds: bool = False):
    """Predicted wall clock and disk for one cell."""
    bytes_per_write = FIELD_BYTES_PER_CELL * plan.n_cells
    if keep_clouds:
        bytes_per_write += CLOUD_BYTES_PER_PARCEL * particles["total_particles"]
    return health.estimate_case(
        case, parcels=particles["total_particles"], cells=plan.n_cells,
        steps=run.n_steps, writes=run.n_writes,
        bytes_per_write=bytes_per_write, model=model)


def main(argv) -> int:
    args = parse_args(argv)

    try:
        study = health.load_study(Path(args.study))
        selected = select(study, args.weight, args.sampling)
    except StudyError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    root = Path(args.study).resolve().parent
    base = root / study.base_case
    model = health.load_cost_model(root / "results" / "cost-model.yaml")

    print(f"study      {args.study}")
    print(f"base case  {base}")
    print(f"matrix     {len(study.enabled_weights())} weight(s) x "
          f"{len(study.enabled_samplings())} sampling(s) = "
          f"{len(study.enabled_cases())} case(s)")
    print(f"Kn         {health.FAMILY_KNUDSEN:g}, fixed for the whole family")
    print()

    if not selected:
        print("Nothing to generate: the filters matched no enabled case.")
        return 1

    # Derive from the TEMPLATE first, so --dry-run can cost the matrix without
    # writing anything. Every case shares the mesh and the time step, so one
    # derivation covers the whole grid except the particle weight.
    try:
        cfg, geom, exit_state, plan, _, _, _, _ = derive(base)
    except (CaiConfigError, ValueError) as exc:
        print(f"generate_cases: {base}: {exc}", file=sys.stderr)
        return 1

    costs = []
    for case in selected:
        run = inflow.derive_run_settings(
            _with_axes(cfg, case), exit_state, geom,
            min_cell_size_m=plan.min_cell_size_m,
            exit_cell_volume_m3=plan.exit_cell_volume_m3)
        particles = mesh.estimate_particles(plan, geom, exit_state,
                                            run.n_equivalent_particles)
        costs.append(estimate(case, plan, run, particles, model))

    print("\n".join(health.cost_report(costs, model)))
    print()

    if args.dry_run:
        print(f"{'case':<16} destination")
        for case in selected:
            print(f"{case.name:<16} {study.case_path(case)}")
        print()
        print("--dry-run: nothing written.")
        return 0

    # Refuse the whole run rather than half of it: generating four cases and
    # then stopping on the fifth leaves a tree nobody can reason about.
    blocked = [c for c in selected
               if case_is_complete(root / study.case_path(c))
               and not args.overwrite]
    if blocked:
        print("generate_cases: these cases already contain results:",
              file=sys.stderr)
        for case in blocked:
            print(f"    {study.case_path(case)}", file=sys.stderr)
        print("\nRegenerating rewrites case.yaml, which would leave results "
              "whose\ninputs no longer describe them. Re-run with --overwrite "
              "to replace\nthem, or filter with --weight / --sampling.",
              file=sys.stderr)
        return 1

    try:
        health.clone_cases(root, study, selected)
    except StudyError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    entries = []
    failed = []
    for case, cost in zip(selected, costs):
        destination = root / study.case_path(case)
        try:
            health.apply_case_parameters(destination, case, study)
        except StudyError as exc:
            print(f"generate_cases: {exc}", file=sys.stderr)
            return 1

        try:
            (case_cfg, case_geom, case_exit, case_plan, run, nozzle, particles,
             report) = derive(destination)
        except (CaiConfigError, ValueError) as exc:
            print(f"generate_cases: {case.name}: {exc}", file=sys.stderr)
            failed.append(case.name)
            continue

        entries.append(health.manifest_entry(
            case, study, case_cfg, case_geom, case_exit, case_plan, run,
            nozzle=nozzle, particles=particles, cost=cost))

        print(f"generated {case.name:<16} "
              f"{case.particles_per_cell:>4g} parcels/cell, "
              f"{case.domain_transits:>4g} transits")
        print(f"          weight {run.n_equivalent_particles:.4e}, "
              f"~{particles['total_particles']:.2e} parcels, "
              f"{run.n_steps:,} steps")
        print(f"          {run.n_writes} writes, "
              f"{health.sampled_writes(run)} of them sampled "
              f"(from t = {run.average_start_s:.4e} s)")
        for result in report.results:
            if result.status != "pass":
                print(f"          {result.status.upper()}: {result.name}: "
                      f"{result.message}")
        if report.failures:
            failed.append(case.name)

    manifest = health.write_manifest(root / study.manifest_name, study, entries,
                                     model)
    print()
    print(f"manifest  {manifest}")

    shared = _held_fixed(manifest)
    if shared:
        print()
        print("generate_cases: THE MATRIX IS NOT COMPARABLE -- these should be "
              "identical in\nevery case and are not:", file=sys.stderr)
        for key in shared:
            print(f"    {key}", file=sys.stderr)
        return 1

    if failed:
        print()
        print(f"generate_cases: {len(failed)} case(s) have failing DSMC "
              f"checks: {failed}", file=sys.stderr)
        print("They were generated, and the manifest records them, but "
              "./Allmesh will\nrefuse to write their dictionaries until "
              "case.yaml is fixed.", file=sys.stderr)
        return 1

    print()
    print("Next:")
    print("    ./AllmeshCases          mesh every generated case")
    print("    ./AllrunCases           run them (costed above)")
    print("    ./AllpostCases          post-process, audit, render, dissect")
    return 0


def _with_axes(cfg, case):
    """A shallow view of *cfg* carrying this cell's two swept values.

    Used only to cost the matrix before anything is written. ``dataclasses``
    are frozen here, so this replaces the two sections rather than mutating
    them -- and it touches only the two, which is the same guarantee
    :func:`plumetools.cai2012.health.apply_case_parameters` gives.
    """
    import dataclasses

    return dataclasses.replace(
        cfg,
        dsmc=dataclasses.replace(
            cfg.dsmc, sampling_domain_transits=case.domain_transits),
        resolution=dataclasses.replace(
            cfg.resolution,
            target_particles_per_cell=case.particles_per_cell),
    )


def _held_fixed(manifest_path: Path) -> list:
    """The values the manifest says are no longer shared, if any."""
    import yaml

    document = yaml.safe_load(Path(manifest_path).read_text(encoding="utf-8"))
    return [f"{key}: {value}"
            for key, value in (document.get("held_fixed") or {}).items()
            if isinstance(value, str) and value.startswith("NOT SHARED")]


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
