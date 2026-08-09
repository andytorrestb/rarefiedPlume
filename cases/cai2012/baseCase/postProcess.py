#!/usr/bin/env python3
"""Extract this case's four validation outputs and score them against Cai.

    ./postProcess.py                 the latest time
    ./postProcess.py --time 0.05     a specific time directory

Produces, under ``results/``:

    centerline.csv        n/n0, U sqrt(beta0), T/T0 vs X/D, DSMC and analytical
    density_plane.csv     n/n0 at the cell centres of a centre plane
    metrics.yaml          max / mean / RMS relative error against Cai
    *.png                 the four figures (needs matplotlib)

Reads the TIME-AVERAGED fields only. A single-timestep rhoN is one step's worth
of parcels -- DSMCCloud::resetFields zeroes it every step -- and comparing that
with an analytical curve measures shot noise.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012 import CAI_REPORTED_MAX_DENSITY_ERROR_PERCENT
    from plumetools.cai2012 import inflow, mesh, post
    from plumetools.cai2012.config import load_case_config
    from plumetools.cai2012.geometry import from_config
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--time", default=None,
                        help="time directory to read (default: the latest)")
    args = parser.parse_args(argv[1:])

    cfg = load_case_config(HERE)
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)

    try:
        sampled = post.read_case(HERE, cfg.gas.mass_kg, time=args.time)
    except post.PostError as exc:
        print(f"postProcess: {exc}", file=sys.stderr)
        return 1

    label = cfg.meta.get("case_name") or f"Kn = {exit_state.knudsen:g}"
    print(f"### {label}")
    print(f"    time                 {sampled.time}")
    print(f"    cells                {sampled.n_cells:,}")

    profile = post.centerline(sampled, cfg, geom, exit_state)
    plane = post.density_plane(sampled, cfg, geom, exit_state)
    metrics = post.centerline_metrics(profile)

    print(f"    centreline samples   {profile['x_over_D'].size} bins "
          f"(of {cfg.post.centerline_points})")
    print(f"    plane cells          {plane['n_over_n0'].size}")
    print()
    print("Centreline against the Cai collisionless solution")
    print(f"  {'quantity':<12} {'max':>10} {'mean':>10} {'RMS':>10}")
    for name in ("density", "velocity", "temperature"):
        print(f"  {name:<12} "
              f"{100 * metrics[f'{name}_max_rel_error']:>9.2f}% "
              f"{100 * metrics[f'{name}_mean_rel_error']:>9.2f}% "
              f"{100 * metrics[f'{name}_rms_rel_error']:>9.2f}%")

    reported = CAI_REPORTED_MAX_DENSITY_ERROR_PERCENT.get(
        float(exit_state.knudsen))
    if reported is not None:
        print()
        print(f"  Cai reports a maximum centreline density error of "
              f"{reported:g}% at Kn = {exit_state.knudsen:g}.")
        print(f"  That is a REFERENCE VALUE, not a tolerance: this is a "
              f"different DSMC code")
        print(f"  on a Cartesian mesh. The trend across the three cases is what "
              f"is being")
        print(f"  validated -- see docs/cai2012-case.md.")

    document = {
        "case": {
            "name": cfg.meta.get("case_name", HERE.name),
            "Kn": float(exit_state.knudsen),
            "S0": float(exit_state.speed_ratio),
            "D_m": float(geom.diameter_m),
            "T0_K": float(exit_state.T0_K),
            "n0_per_m3": float(exit_state.number_density_per_m3),
            "lambda0_m": float(exit_state.mean_free_path_m),
            "cell_over_mean_free_path": float(plan.cell_over_mfp),
            "n_cells": int(sampled.n_cells),
            "time": sampled.time,
            "cai_reported_max_percent": reported,
        },
        "centerline": {k: v for k, v in metrics.items()},
    }

    results = HERE / "results"
    for path in post.write_results(results, profile, plane, document):
        print(f"    wrote {path}")

    exact_plane = post.analytical_plane(cfg, geom, exit_state)
    figures = post.plot_case(results, profile, plane, exact_plane, cfg, label)
    if figures:
        for figure in figures:
            print(f"    wrote {figure}")
    else:
        print('    no plots: matplotlib is not installed '
              '(pip install -e ".[plots]")')
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
