#!/usr/bin/env python3
"""Generate the interface track of the Cai 2012 Knudsen extension, idempotently.

    ./generate_interfaces.py               create what is missing, keep what matches
    ./generate_interfaces.py --dry-run     the plan and every verdict; write nothing
    ./generate_interfaces.py --report      also write each case's
                                           results/interface-report.yaml, with the
                                           (core angle x threshold) decision table
    ./generate_interfaces.py --refresh     copy updated scripts / provenance into
                                           cases whose inputs are unchanged
    ./generate_interfaces.py --regenerate  replace CONFLICTING cases with no output
    ./generate_interfaces.py --overwrite   replace conflicting cases EVEN WITH
                                           RESULTS -- which destroys them

What it does:

  1. read study.yaml with cases/cai2012's loader and the extension's policy,
     unchanged, and interfaces.yaml beside it
  2. stage every variant: baseCase/ + baseCaseInterface/, case.yaml by the
     parent's apply_case_parameters, interface.yaml from interfaces.yaml
  3. compare with Interface/ (case.yaml AND interface.yaml are inputs): create,
     keep, report stale, refuse conflicts -- results are never touched
  4. derive each case through plumetools.cai2012ext.interface_case -- the same
     path Allmesh, runInflow.py and the PBS jobs take -- and write
     interface-manifest.yaml

Cases/, manifest.yaml and the direct-exit cases are not read or written.
Exit status 0 when the tree matches interfaces.yaml, even when cases are
infeasible or pending validation: those are recorded verdicts.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

try:
    from plumetools.cai2012.study import StudyError, load_study
    from plumetools.cai2012ext import feasibility, generation, interface_case
    from plumetools.cai2012ext import interface_study, policy
    from plumetools.cai2012ext.interface import InterfaceConfigError
    from plumetools.cai2012ext.provenance import git_provenance
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')

HEALTH = REPO / "cases" / "cai2012-health"


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generate the interface track of the Cai 2012 extension.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--report", action="store_true",
                        help="write results/interface-report.yaml in each case")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--study", default=str(HERE / "study.yaml"))
    parser.add_argument("--interfaces", default=None,
                        help="default: interfaces.yaml beside study.yaml")
    return parser.parse_args(argv[1:])


LABELS = {"create": "create", "keep": "keep (identical)",
          "keep-stale": "keep (STALE scripts/provenance; --refresh)",
          "refresh": "refresh scripts/provenance", "replace": "REPLACE",
          "refuse": "REFUSED -- conflict"}


def main(argv) -> int:
    args = parse_args(argv)
    study_path = Path(args.study).resolve()
    root = study_path.parent
    interfaces_path = Path(args.interfaces).resolve() if args.interfaces \
        else root / "interfaces.yaml"
    try:
        study = load_study(study_path)
        policy.check_study(study)
        layout = interface_study.load(interfaces_path, study)
    except (StudyError, InterfaceConfigError) as exc:
        print(f"generate_interfaces: {exc}", file=sys.stderr)
        return 1

    print(f"interfaces {interfaces_path}")
    print(f"tree       {root / layout.tree}  ({len(layout.variants)} variants; "
          f"Cases/ is not touched)")
    print()
    model, origin = feasibility.load_cost_model(root, HEALTH)

    with tempfile.TemporaryDirectory(prefix=".staging-", dir=root) as scratch:
        states = generation.decide(
            interface_study.classify(root, layout, study, Path(scratch)),
            refresh=args.refresh, regenerate=args.regenerate,
            overwrite=args.overwrite)
        print(f"{'case':<26} {'group':<11} {'Kn':>7}  action")
        for state in states:
            variant = state.case
            print(f"{variant.name:<26} {variant.group:<11} {variant.knudsen:>7g}  "
                  f"{LABELS[state.action]}"
                  + ("  [holds solver output]" if state.has_results else ""))
            for line in (state.inputs if state.action in ("refuse", "replace")
                         else state.auxiliary)[:6]:
                print(f"{'':<28}{line}")
        print()

        targets = [(s.case, s.staged if args.dry_run else s.destination)
                   for s in states]
        if not args.dry_run:
            refused = [s for s in states if s.blocking]
            if refused:
                for state in refused:
                    why = ("it holds solver output; --overwrite would DESTROY it"
                           if state.has_results else "re-run with --regenerate")
                    print(f"generate_interfaces: {state.case.name} conflicts with "
                          f"interfaces.yaml: {why}", file=sys.stderr)
                print("Nothing was written.", file=sys.stderr)
                return 1
            generation.apply(states)
        entries = []
        for variant, case_dir in targets:
            try:
                derived = interface_case.derive(case_dir)
            except (InterfaceConfigError, policy.ExtensionError, ValueError) as exc:
                print(f"generate_interfaces: {variant.name}: {exc}", file=sys.stderr)
                return 1
            assessment = interface_case.assess(derived, cost_model=model,
                                               cost_origin=origin)
            entry = interface_case.manifest_entry(
                root if not args.dry_run else Path(scratch), variant.name, derived,
                assessment, group=variant.group)
            if args.dry_run:
                entry["path"] = f"{layout.tree}/{variant.name}"
            entries.append(entry)
            print(f"{variant.name}  ({variant.group}, Kn {variant.knudsen:g}, "
                  f"B = {derived.settings.bird_threshold:g})")
            print("\n".join(interface_case.describe(derived, assessment)))
            if args.report and not args.dry_run:
                path = interface_case.write_report(
                    case_dir, interface_case.report(derived, assessment,
                                                    decision_table=True))
                print(f"          wrote {path}")
            print()

    if args.dry_run:
        print("--dry-run: nothing written.")
        return 0
    manifest = interface_study.write_manifest(
        root / layout.manifest_name, layout, entries, code=git_provenance(REPO))
    summary = interface_study.validation_summary(entries, layout)
    print(f"manifest  {manifest}")
    if summary and not summary["overlap_study_possible"]:
        floor = summary["B_at_smallest_admissible_radius"]
        print()
        print(f"VALIDATION: no admissible sphere at Kn = {summary['knudsen']:g} for "
              f"thresholds {summary['thresholds']}.")
        print(f"  B on the axis at the smallest admissible radius is "
              f"{floor['axis']:.4g}" if floor["axis"] is not None else
              "  (B floor not evaluated)")
        print("  The overlap study cannot run as declared. The thresholds are a "
              "study decision\n  (docs/cai2012-kn-interface.md, section 5); nothing was "
              "changed to make one fit.")
    print()
    print("Next (per case, docs/cai2012-kn-interface.md, section 9):")
    print("    cd Interface/<case> && bash ./Allmesh --dict-only     # no OpenFOAM")
    print("    cd Interface/<case> && qsub runMesh.sh && qsub runCase.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
