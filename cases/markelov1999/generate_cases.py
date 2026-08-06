#!/usr/bin/env python3
"""Generate the AIAA 99-3455 case hierarchy from study.yaml.

    ./generate_cases.py                  generate every enabled case
    ./generate_cases.py --dry-run        show what would be generated
    ./generate_cases.py --pressures 5 25 only these reservoir pressures
    ./generate_cases.py --overwrite      replace cases that already have results

No arguments are needed for normal use.

What it does, per case:

  1. clone baseCase into Cases/<gap>/<case>   (CaseFoam if installed, else the
     equivalent built-in copy -- same tree either way)
  2. rewrite the case-local case.yaml with this case's reservoir pressure and
     particle weight -- STRUCTURALLY, by loading and re-emitting YAML, never by
     string substitution
  3. record everything in manifest.yaml

The OpenFOAM dictionaries are NOT written here. Each case writes its own from
its own case.yaml when ./Allmesh runs, so there is one authoritative source of
physical inputs per case and no chance of a stale dictionary disagreeing with it.

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
    from plumetools.config import load_case_config
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')

from plumetools.markelov1999 import mesh as mesh_module
from plumetools.markelov1999 import resolution as resolution_module
from plumetools.markelov1999.geometry import from_config
from plumetools.markelov1999.study import (
    StudyError,
    apply_case_parameters,
    case_is_complete,
    clone_base_case,
    filter_by_pressure,
    generation_backend,
    load_study,
    write_manifest,
)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generate the AIAA 99-3455 case hierarchy from study.yaml.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument(
        "--dry-run", action="store_true",
        help="list the cases that would be generated and stop")
    parser.add_argument(
        "--pressures", nargs="+", type=float, metavar="PSI",
        help="only these reservoir pressures (default: every enabled case)")
    parser.add_argument(
        "--overwrite", action="store_true",
        help="replace cases that already contain results")
    parser.add_argument(
        "--study", default=str(HERE / "study.yaml"),
        help="path to study.yaml (default: alongside this script)")
    return parser.parse_args(argv[1:])


def derive_particle_weight(case_dir: Path, pressure_pa: float):
    """Particle weight and resolution estimate for one generated case.

    Loads the case's own ``case.yaml`` -- after the pressure has been written --
    so the weight is derived from the same configuration the case will run with,
    not from the template's.

    Returns ``(weight, estimate)``, or ``(None, None)`` when the case pins the
    weight explicitly or ``resolution.report_only`` is set.
    """
    cfg = load_case_config(case_dir)
    if cfg.resolution.report_only and cfg.dsmc.n_equivalent_particles is not None:
        return None, None

    geom = from_config(cfg)
    cell_sizes = mesh_module.surface_cell_sizes(cfg)
    estimate = resolution_module.estimate(
        cfg, geom, cell_sizes,
        n_equivalent_particles=cfg.dsmc.n_equivalent_particles)
    return estimate.n_equivalent_particles, estimate


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
        selected = filter_by_pressure(study.enabled_cases(), args.pressures)
    except StudyError as exc:
        print(f"generate_cases: {exc}", file=sys.stderr)
        return 1

    disabled = [c for c in study.all_cases() if not c.enabled]
    backend = generation_backend()

    print(f"study      {args.study}")
    print(f"base case  {base}")
    print(f"gap        {study.gap_in:g} in -> {study.gap_dir}/")
    print(f"backend    {backend}"
          + ("  (CaseFoam is not installed; the built-in clone produces the "
             "same tree)" if backend == "builtin" else ""))
    print()

    if not selected:
        print("Nothing to generate: every case is disabled, or the --pressures "
              "filter matched none.")
        return 1

    print(f"{'case':<10} {'psi':>6} {'Pa':>16}  destination")
    for case in selected:
        print(f"{case.name:<10} {case.pressure_psi:>6g} {case.pressure_pa:>16.6f}  "
              f"{study.case_path(case)}")
    for case in disabled:
        print(f"{case.name:<10} {case.pressure_psi:>6g} {'':>16}  "
              f"SKIPPED (enabled: false)")
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
        print("\nRegenerating rewrites case.yaml, which would leave results whose\n"
              "inputs no longer describe them. Re-run with --overwrite to replace\n"
              "them, or with --pressures to select the others.", file=sys.stderr)
        return 1

    entries = []
    for case in selected:
        destination = root / study.case_path(case)
        if destination.exists() and args.overwrite:
            shutil.rmtree(destination)

        clone_base_case(base, destination)
        applied = apply_case_parameters(destination, case, study)

        weight, estimate = derive_particle_weight(destination, case.pressure_pa)
        if weight is not None:
            applied.update(apply_case_parameters(
                destination, case, study, n_equivalent_particles=weight))

        entry = {
            "name": case.name,
            "path": str(study.case_path(case)).replace("\\", "/"),
            "pressure_psi": case.pressure_psi,
            "pressure_pa": case.pressure_pa,
            "gap_in": study.gap_in,
            "n_equivalent_particles": applied.get("n_equivalent_particles"),
        }
        if estimate is not None:
            entry["estimated_particles_per_cell"] = {
                r.region.name: round(r.particles_per_cell, 3)
                for r in estimate.regions
            }
            entry["estimated_total_particles"] = float(f"{estimate.total_particles:.6g}")
            entry["sizing_region"] = estimate.sizing_region
        if case.note:
            entry["note"] = case.note
        entries.append(entry)

        print(f"generated {case.name:<10} {case.pressure_psi:>6g} psi"
              + (f"   weight {weight:.4e}" if weight else ""))
        if estimate is not None:
            occupancy = "  ".join(
                f"{r.region.name} {r.particles_per_cell:.1f}" for r in estimate.regions)
            print(f"          particles/cell: {occupancy}   "
                  f"total {estimate.total_particles:.3e}")

    manifest = write_manifest(root / study.manifest_name, study, entries, backend)
    print()
    print(f"manifest  {manifest}")
    print()
    print("Next:")
    print(f"    ./AllmeshCases          mesh every generated case")
    print(f"    ./AllrunCases           run them")
    print(f"    ./AllpostCases          post-process and build the study table")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
