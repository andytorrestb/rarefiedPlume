#!/usr/bin/env python3
"""Generate the Cai 2012 Knudsen extension: four cases, idempotently.

    ./generate_cases.py                 create what is missing, keep what matches
    ./generate_cases.py --dry-run       the plan and the feasibility; write nothing
    ./generate_cases.py --refresh       copy updated runner scripts into cases
                                        whose physics is unchanged
    ./generate_cases.py --regenerate    replace CONFLICTING cases that hold no
                                        solver output
    ./generate_cases.py --overwrite     replace conflicting cases EVEN IF THEY
                                        HOLD RESULTS -- which destroys them

No arguments are needed for normal use, and running it twice is free: a case
whose staged copy matches what is on disk is left alone, mesh and results
included. See plumetools/cai2012ext/generation.py.

What it does:

  1. read study.yaml with cases/cai2012's own loader, and refuse anything that
     is not exactly Kn = 1e-4, 1e-3, 1, 10 at the baseline parcel level -- or
     that mentions Cai's Kn = 0.01, 0.1 or 100 at all, even disabled
  2. check baseCase/case.yaml against cases/cai2012-health's template, key by
     key; only the four declared [EXTENSION] keys may differ
  3. stage every case in a scratch directory with the parent family's own
     apply_case_parameters, and compare it with the disk
  4. create the new ones, keep the identical ones, REPORT the conflicts
  5. derive everything each case will run at, assess it against the production
     criteria, and write manifest.yaml

The OpenFOAM dictionaries are not written here; each case writes its own from
its case.yaml when ./Allmesh runs (./AllmeshCases --dict-only needs no
OpenFOAM). Exit status is 0 when the tree matches the study, even if some cases
are unsuitable for production -- that is a recorded verdict, not an error.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

try:
    from plumetools.cai2012.config import CaiConfigError, load_case_config
    from plumetools.cai2012.study import StudyError, load_study
    from plumetools.cai2012ext import EXCLUDED_KNUDSEN, feasibility, generation, policy
    from plumetools.cai2012ext.provenance import git_provenance
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')

HEALTH = REPO / "cases" / "cai2012-health"


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generate the Cai 2012 Knudsen extension.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="stage, compare and assess; write nothing")
    parser.add_argument("--refresh", action="store_true",
                        help="update runner scripts / meta in cases whose "
                             "physics is unchanged")
    parser.add_argument("--regenerate", action="store_true",
                        help="replace conflicting cases with no solver output")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace conflicting cases even with results")
    parser.add_argument("--study", default=str(HERE / "study.yaml"),
                        help="path to study.yaml (default: alongside this script)")
    return parser.parse_args(argv[1:])


def inheritance(template_cfg) -> dict | None:
    """The template against cases/cai2012-health's, or None if that is absent."""
    health_template = HEALTH / "baseCase"
    if not (health_template / "case.yaml").is_file():
        return None
    report = policy.inheritance_report(template_cfg,
                                       load_case_config(health_template))
    problems = policy.check_health_study_rows(HEALTH / "study.yaml")
    return {
        "source": "cases/cai2012-health/baseCase/case.yaml + study.yaml "
                  "(ppc040, s4p5)",
        "n_inherited": len(report["inherited"]),
        "overridden": report["overridden"],
        "undeclared_differences": report["undeclared"],
        "health_study_problems": problems,
    }


def print_states(states, root) -> None:
    labels = {"create": "create", "keep": "keep (identical)",
              "keep-stale": "keep (STALE scripts/meta; --refresh)",
              "refresh": "refresh scripts/meta", "replace": "REPLACE",
              "refuse": "REFUSED -- conflict"}
    print(f"{'case':<20} {'Kn':>8} {'parcels':>8}  action")
    for state in states:
        print(f"{state.case.name:<20} {state.case.knudsen:>8g} "
              f"{state.case.multiplier:>7g}x  {labels[state.action]}"
              + ("  [holds solver output]" if state.has_results else ""))
        if state.action in ("refuse", "replace"):
            for line in state.inputs[:8]:
                print(f"{'':<20}   input  {line}")
        if state.action in ("keep-stale", "refresh"):
            for line in state.auxiliary[:4]:
                print(f"{'':<20}   aux    {line}")
    print()


def print_table(entries) -> None:
    print(f"{'case':<20} {'Kn':>7} {'n0 [1/m3]':>11} {'cells':>10} "
          f"{'cell/lam':>9} {'dt/tau0':>8} {'parcels':>9} {'steps':>7}  verdict")
    for entry in entries:
        f = entry["feasibility"]
        print(f"{entry['name']:<20} {entry['Kn']:>7g} {entry['n0_per_m3']:>11.4e} "
              f"{entry['n_cells']:>10,} "
              f"{f['resolution']['cell_over_mean_free_path']:>9.3g} "
              f"{f['time_step']['dt_over_exit_collision_time']:>8.3g} "
              f"{f['cost']['estimated_parcels']:>9.2e} {entry['n_steps']:>7,}  "
              f"{f['status']}")
    print()


def main(argv) -> int:
    args = parse_args(argv)
    study_path = Path(args.study).resolve()
    root = study_path.parent

    try:
        study = load_study(study_path)
        policy.check_study(study)
        template_cfg = load_case_config(root / study.base_case)
        policy.check_template(template_cfg, f"{study.base_case}/case.yaml")
    except (StudyError, CaiConfigError) as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    inherited = inheritance(template_cfg)
    if inherited is not None and (inherited["undeclared_differences"]
                                  or inherited["health_study_problems"]):
        print("generate_cases: baseCase/case.yaml no longer inherits "
              "cases/cai2012-health's settings:", file=sys.stderr)
        for key, values in inherited["undeclared_differences"].items():
            print(f"    {key}: health {values['health']!r}, here "
                  f"{values['extension']!r}", file=sys.stderr)
        for problem in inherited["health_study_problems"]:
            print(f"    {problem}", file=sys.stderr)
        print("Either restore the inherited value or declare the override in "
              "plumetools/cai2012ext/policy.py::DELIBERATE_OVERRIDES, with its "
              "reason.", file=sys.stderr)
        return 1

    cases = study.expand()
    if len(cases) != policy.expected_case_count():
        print(f"generate_cases: the matrix has {len(cases)} cases, not "
              f"{policy.expected_case_count()}", file=sys.stderr)
        return 1

    model, origin = feasibility.load_cost_model(root, HEALTH)

    print(f"study      {study_path}")
    print(f"base case  {root / study.base_case}")
    print(f"matrix     {len(study.enabled_cases())} Knudsen numbers x "
          f"{len(study.enabled_particle_levels())} particle levels = "
          f"{len(cases)} cases")
    print(f"excluded   Kn {list(EXCLUDED_KNUDSEN)} -- never generated or run")
    if inherited is not None:
        print(f"inherits   {inherited['n_inherited']} settings from "
              f"cases/cai2012-health; overrides "
              f"{sorted(inherited['overridden'])}")
    print(f"cost model {model.source} ({origin})")
    print()

    strays = generation.stray_directories(root, study, cases)
    for stray in strays:
        flag = ("  EXCLUDED Kn -- not part of this study and never scheduled"
                if stray["excluded"] else "  not part of this study; ignored")
        print(f"stray      {stray['path']} (Kn {stray['Kn']}){flag}")
    if strays:
        print()

    with tempfile.TemporaryDirectory(prefix=".staging-", dir=root) as scratch:
        states = generation.decide(
            generation.classify(root, study, Path(scratch), cases),
            refresh=args.refresh, regenerate=args.regenerate,
            overwrite=args.overwrite)
        print_states(states, root)

        if args.dry_run:
            entries = []
            for state in states:
                entry, _ = generation.build_entry(
                    root, study, study_path, state.case, state.staged,
                    cost_model=model, cost_origin=origin)
                entries.append(entry)
            print_table(entries)
            print("--dry-run: nothing written.")
            return 0

        refused = [s for s in states if s.blocking]
        if refused:
            print("generate_cases: these cases conflict with study.yaml / "
                  "baseCase and were NOT changed:", file=sys.stderr)
            for state in refused:
                why = ("it holds solver output; --overwrite would DESTROY it"
                       if state.has_results else "re-run with --regenerate")
                print(f"    {state.destination.relative_to(root)}: {why}",
                      file=sys.stderr)
            print("Nothing was written. manifest.yaml is unchanged.",
                  file=sys.stderr)
            return 1

        for state in states:
            if state.action == "replace" and state.has_results:
                print(f"generate_cases: --overwrite: DESTROYING the results in "
                      f"{state.destination.relative_to(root)}", file=sys.stderr)
        generation.apply(states)

    entries, hard_failures = [], []
    for case in cases:
        destination = root / study.case_path(case)
        try:
            entry, lines = generation.build_entry(
                root, study, study_path, case, destination,
                cost_model=model, cost_origin=origin)
        except (CaiConfigError, ValueError) as exc:
            print(f"generate_cases: {case.name}: {exc}", file=sys.stderr)
            return 1
        entries.append(entry)
        print(f"{case.name:<20} Kn {case.knudsen:g}, {case.multiplier:g}x "
              f"parcels, weight {entry['n_equivalent_particles']:.4e} "
              f"(= {entry['baseline_n_equivalent_particles']:.4e} / "
              f"{case.multiplier:g})")
        print("\n".join(lines))
        if entry["dsmc_checks"]["n_failures"]:
            hard_failures.append(case.name)

    manifest = generation.write_manifest(
        root / study.manifest_name, study, entries, study_path=study_path,
        inheritance=inherited, code=git_provenance(REPO), strays=strays)
    print()
    print_table(entries)
    print(f"manifest  {manifest}")

    _, _, problems = generation.held_fixed(entries)
    if problems:
        print("\ngenerate_cases: values that must be shared are NOT:",
              file=sys.stderr)
        for problem in problems:
            print(f"    {problem}", file=sys.stderr)
        return 1
    if hard_failures:
        print(f"\ngenerate_cases: {hard_failures} fail a hard DSMC check; "
              f"./Allmesh will refuse them.", file=sys.stderr)
        return 1

    unsuitable = [e["name"] for e in entries
                  if e["production_status"] != feasibility.STATUS_CANDIDATE]
    if unsuitable:
        print()
        print(f"{len(unsuitable)} of {len(entries)} cases are UNSUITABLE FOR "
              f"PRODUCTION as configured (see 'unresolved' in the manifest).")
        print("They are generated and can be meshed; ./AllrunCases will not run "
              "them without\n--allow-unsuitable, and their outputs are "
              "labelled non-production.")
    print()
    print("Next:")
    print("    ./AllmeshCases --dict-only   write every case's dictionaries "
          "(no OpenFOAM)")
    print("    ./validate_cases.py          check the tree against the study")
    print("    ./AllmeshCases               mesh (OpenFOAM)")
    print("    ./AllrunCases --list         what would run, for a job array")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
