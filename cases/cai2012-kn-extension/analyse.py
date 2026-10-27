#!/usr/bin/env python3
"""Assemble the extension: the study table and the Knudsen comparison.

    ./analyse.py

Runs after every case has been post-processed (./AllpostCases does it last).
Reads each case's results/ and the manifest, and writes under results/:

    study-table.csv            one row per case -- all four, including any not
                               run yet (state says so) -- with the production
                               verdict, the pre-run resolution, the measured
                               resolution from Alldiagnose.py (valid volume,
                               exit KnD and Delta/lambda, the fraction coarser
                               than lambda/3, the cells a resolved mesh needs),
                               the centreline / r/D = 10 / centre-plane errors,
                               measured occupancy and the stationarity verdict
    study-table.md             a readable subset of it
    knudsen-comparison.png     the four Kn on one figure: centreline n/n0, U,
                               T and the r/D = 10 density and mass flux against
                               Cai's collisionless curve, and the centreline
                               Delta_max/lambda beside them

EVERY comparison is against Cai's collisionless ANALYTICAL solution. There is
no published DSMC result at any Kn in this study, and the table has no "Cai
reported" column for that reason.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    import yaml

    from plumetools.cai2012ext import policy, summary
    from plumetools.cai2012.study import load_study
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    study = load_study(HERE / "study.yaml")
    path = HERE / study.manifest_name
    if not path.is_file():
        sys.exit(f"analyse: no {study.manifest_name}. Generate the cases first:\n    ./generate_cases.py")
    manifest = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    problems = policy.check_manifest(manifest)
    if problems:
        print("analyse: manifest.yaml breaks the extension policy:",
              file=sys.stderr)
        for problem in problems:
            print(f"    {problem}", file=sys.stderr)
        return 2

    entries = manifest.get("cases") or []
    rows = [summary.study_row(HERE, entry) for entry in entries]
    results = HERE / "results" / ("full-suite" if policy.suite_name(study) == "full_decades" else "")
    written = [summary.write_table(results / "study-table.csv", rows),
               summary.write_markdown(results / "study-table.md", rows)]
    written += summary.plot_knudsen(HERE, entries, results)

    print(f"{'case':<20} {'status':<26} {'state':<36} {'n mean%':>8} "
          f"{'F mean%':>8}  stationarity")
    for row in rows:
        def show(value):
            return f"{value:8.2f}" if isinstance(value, float) else f"{'':>8}"
        print(f"{row['case']:<20} {str(row['production_status']):<26} "
              f"{row['state']:<36} {show(row['centerline_density_mean_rel_pct'])} "
              f"{show(row['angular_mass_flux_mean_rel_pct'])}  "
              f"{row['stationarity'] or ''}")
    print()
    print("Reference for every row: Cai's collisionless analytical solution. "
          "No external-study result fields are imported.")
    print()
    for item in written:
        print(f"    wrote {item}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
