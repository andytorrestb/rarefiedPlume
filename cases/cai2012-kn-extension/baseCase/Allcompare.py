#!/usr/bin/env python3
"""Compare this case with Cai's collisionless solution, off the axis too.

    ./Allcompare.py                    the latest sampled frame
    ./Allcompare.py --time 0.06        a specific time directory
    ./Allcompare.py --radius 10        the angular sphere, in D (default 10)

Run after ./Allpost, which writes the cell centres and volumes and the four
outputs every Cai family produces (centreline n/n0, U sqrt(beta0), T/T0 and the
n/n0 contour). This adds, under results/:

    centerline_extended.csv   the centreline with p/p0 and absolute errors
    angular_r10D.csv          n/n0, u_R, u_theta, T/T0, p/p0 and the mass flux
                              F on the r/D = 10 sphere, DSMC / same-cells
                              analytical / exact analytical
    meridional_plane.csv      every centre-plane cell: DSMC and analytical
                              n/n0, u_x, u_r, T/T0, p/p0
    comparison.yaml           every metric, every exclusion count, the
                              definitions, and this case's production verdict
    *.png                     the figures, on fixed normalised colour scales

Definitions, normalisations, the projection of the 3-D velocity onto the
meridional plane and the error rules are in plumetools/cai2012ext/comparison.py
and docs/cai2012-kn-extension.md.

THE REFERENCE IS CAI'S ANALYTICAL COLLISIONLESS SOLUTION, nothing else. There
are no published DSMC results at this Knudsen number, and Cai's reported errors
at Kn = 100, 0.1 and 0.01 are not used as tolerances here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    import yaml

    from plumetools.cai2012 import post
    from plumetools.cai2012ext import comparison, convergence
    from plumetools.cai2012ext import figures, interface_case
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--time", default=None,
                        help="time directory to read (default: the latest "
                             "sampled one)")
    parser.add_argument("--mesh-time", default=None,
                        help="time directory to read C and V from")
    parser.add_argument("--radius", type=float,
                        default=comparison.ANGULAR_RADIUS_OVER_D,
                        help="angular sphere radius in D (default 10)")
    parser.add_argument("--shell", type=float,
                        default=comparison.SHELL_HALF_WIDTH_OVER_D,
                        help="shell half-width in D (default 0.5)")
    parser.add_argument("--bin", type=float, default=comparison.ANGULAR_BIN_DEG,
                        help="polar-angle bin in degrees (default 2.5)")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv[1:])

    # The direct-exit derivation, or -- with an enabled interface.yaml -- the
    # interface one; the comparison is with Cai's exit solution either way.
    context = interface_case.post_context(HERE)
    if context.refusal:
        print(f"Allcompare: {context.refusal}", file=sys.stderr)
        return 2
    cfg, geom, exit_state, run = (context.cfg, context.geom, context.exit_state,
                                  context.run)
    verdict = context.verdict
    ready = verdict["production_ready"]
    banner = figures.status_banner(ready)
    label = cfg.meta.get("case_name") or HERE.name

    times = convergence.sampled_frames(HERE, run.average_start_s)
    time = args.time or (times[-1] if times else None)
    if time is None:
        print(f"Allcompare: no sampled frame at or after "
              f"t = {run.average_start_s:.6e} s. Has the case reached its "
              f"averaging window?", file=sys.stderr)
        return 1
    mesh_time = args.mesh_time or convergence.geometry_time(HERE, [time])
    if mesh_time is None:
        print("Allcompare: no cell centres in this case. Run ./Allpost first "
              "(postProcess -func writeCellCentres).", file=sys.stderr)
        return 1

    try:
        sampled = post.read_case(HERE, cfg.gas.mass_kg, time=time,
                                 mesh_time=mesh_time)
    except post.PostError as exc:
        print(f"Allcompare: {exc}", file=sys.stderr)
        return 1
    parcel_path = HERE / time / "dsmcRhoNMean"
    parcels = (post.read_internal_field(parcel_path)
               if parcel_path.is_file() else None)

    print(f"### {label}")
    print(f"    verdict              {verdict['status'].upper()}")
    print(f"    time                 {sampled.time}  (averaging from "
          f"{run.average_start_s:.6e} s, "
          f"{(float(sampled.time) - run.average_start_s) / run.domain_transit_s:.3f}"
          f" transits)")
    print(f"    cells                {sampled.n_cells:,}")

    s0 = float(exit_state.speed_ratio)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    profile, centerline_metrics = comparison.centerline_extended(profile, s0)

    angular = comparison.angular_profile(
        sampled, geom, exit_state, radius_over_D=args.radius,
        half_width_over_D=args.shell, bin_deg=args.bin, parcel_count=parcels)
    angular_metrics = comparison.angular_metrics(angular, s0)

    plane = comparison.meridional_plane(sampled, cfg, geom, exit_state)
    plane_metrics = comparison.plane_metrics(plane, s0)

    print()
    print("Against Cai's collisionless ANALYTICAL solution (matched estimators)")
    print(f"  {'quantity':<28} {'valid':>6} {'max rel':>9} {'mean rel':>9} "
          f"{'max abs':>9}")
    for name, metrics in (("density", centerline_metrics),
                          ("velocity", centerline_metrics),
                          ("temperature", centerline_metrics),
                          ("pressure", centerline_metrics),
                          ("angular_density", angular_metrics),
                          ("angular_mass_flux", angular_metrics),
                          ("plane_density", plane_metrics),
                          ("plane_axial_velocity", plane_metrics),
                          ("plane_radial_velocity", plane_metrics),
                          ("plane_temperature", plane_metrics),
                          ("plane_pressure", plane_metrics)):
        rel = metrics.get(f"{name}_max_rel_error")
        mean = metrics.get(f"{name}_mean_rel_error")
        shown = "centreline " + name if metrics is centerline_metrics else name
        print(f"  {shown:<28} {metrics[f'{name}_n_valid']:>6} "
              + (f"{100 * rel:>8.2f}% {100 * mean:>8.2f}% " if rel is not None
                 else f"{'--':>9} {'--':>9} ")
              + f"{metrics[f'{name}_max_abs_error']:>9.4f}")
    print(f"  mass flux through r/D = {args.radius:g}: DSMC "
          f"{angular['integrals']['F']:.4f}, same cells "
          f"{angular['integrals']['F_analytical']:.4f}, collisionless "
          f"{angular['integrals']['F_exact']:.4f} (of the exit flow)")
    if context.nozzle_area_error is not None:
        print(f"  meshed nozzle area / pi R0^2 = "
              f"{1.0 + context.nozzle_area_error:.4f}")
    else:
        print(f"  INTERFACE MODE: R < {context.interface['R_interface_over_D']:.4f} D "
              f"is not simulated; nothing there is compared")

    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    written = [
        post.write_csv(results / "centerline_extended.csv", {
            key: profile[key] for key in (
                "x_over_D", "n_cells",
                "n_over_n0", "n_over_n0_analytical", "n_over_n0_analytical_axis",
                "U_sqrt_beta0", "U_sqrt_beta0_analytical",
                "U_sqrt_beta0_analytical_axis",
                "T_over_T0", "T_over_T0_analytical", "T_over_T0_analytical_axis",
                "p_over_p0", "p_over_p0_analytical", "p_over_p0_analytical_axis")}),
        post.write_csv(results / f"angular_r{args.radius:g}D.csv", {
            "theta_deg": angular["theta_deg"], "n_cells": angular["n_cells"],
            "solid_angle_sr": angular["solid_angle"],
            **{f"{key}{suffix}": angular[f"{key}{suffix}"]
               for key in ("n_over_n0", "F", "uR", "utheta", "T_over_T0",
                           "p_over_p0", "uphi")
               for suffix in ("", "_analytical", "_exact")},
            **({"parcels": angular["parcels"]} if "parcels" in angular else {}),
        }),
        post.write_csv(results / "meridional_plane.csv", {
            "x_over_D": plane["x_over_D"],
            f"{plane['lateral_axis'].lower()}_over_D": plane["lateral_over_D"],
            "normal_m": plane["measured"]["centres"][
                :, 1 if plane["lateral_axis"] == "Z" else 2],
            "nozzle_lip_zone": plane["lip"].astype(float),
            "occupied": plane["measured"]["occupied"].astype(float),
            **{f"{key}{suffix}": source[key]
               for key in ("n_over_n0", "ux", "ur", "uphi", "uR", "T_over_T0",
                           "p_over_p0")
               for suffix, source in (("", plane["measured"]),
                                      ("_analytical", plane["exact"]))
               if key in source},
        }),
    ]

    document = {
        "case": {
            "name": label,
            "Kn": float(exit_state.knudsen),
            "cai_case": cfg.meta.get("cai_case"),
            "particle_level": cfg.meta.get("particle_level"),
            "numerical_particle_multiplier": float(
                cfg.dsmc.numerical_particle_multiplier),
            "time": sampled.time,
            "mesh_time": mesh_time,
            "average_start_s": float(run.average_start_s),
            "averaged_domain_transits": (float(sampled.time) - run.average_start_s)
            / run.domain_transit_s,
            "n_cells": int(sampled.n_cells),
            "nozzle_area_ratio": (None if context.nozzle_area_error is None
                                  else 1.0 + float(context.nozzle_area_error)),
            "mode": context.mode,
        },
        "interface": context.interface,
        "production": {
            "status": verdict["status"],
            "production_ready": bool(ready),
            "unresolved": verdict["unresolved"],
            "non_production_marker": (HERE / "NON-PRODUCTION").is_file(),
            "statement": ("production candidate: still requires the convergence "
                          "and occupancy evidence of Allconverge.py and "
                          "Allaudit.py" if ready else
                          "UNSUITABLE FOR PRODUCTION: nothing here is a "
                          "validated result"),
        },
        "reference": {
            "compared_with": "Cai's collisionless analytical solution "
                             "(plumetools.cai2012.analytical.moments)",
            "published_dsmc": None,
            "note": "No DSMC data are published at this Kn. Cai's reported "
                    "centreline errors at Kn = 100, 0.1 and 0.01 are not "
                    "tolerances for this case and are not applied.",
        },
        "definitions": {
            "normalisation": "n/n0, u sqrt(beta0) with beta0 = 1/(2 R T0), "
                             "T/T0, p/p0",
            "pressure": "p = n k T, T the local translational temperature "
                        "from 3 R T = <v^2> - |u|^2; p/p0 = (n/n0)(T/T0)",
            "u_r": "cylindrical radial (y u_y + z u_z)/r, away from the axis",
            "u_R": "spherical radial from the exit centre, "
                   "cos(theta) u_x + sin(theta) u_r",
            "F": "n u_R R^2 / (Gamma0 pi R0^2) per steradian; "
                 "integral F 2 pi sin(theta) dtheta = 1 when collisionless",
            "estimators": "analytical values at the same cell centres with "
                          "the same weights (volume for n, p, F; molecules "
                          "for u, T)",
            "density_floor_over_n0": comparison.DENSITY_FLOOR_OVER_N0,
            "angular_floor_fraction_of_peak": comparison.ANGULAR_FLOOR_FRACTION,
            "lip_exclusion_over_D": comparison.LIP_EXCLUSION_OVER_D,
            "flux_floor_fraction_of_peak": comparison.FLUX_FLOOR_FRACTION,
            "relative_error": "|DSMC - analytical| / |analytical| where valid",
            "absolute_error": "|DSMC - analytical| / scale (1 for n, T, p; "
                              "S0 for velocities; peak F for the flux)",
        },
        "centerline": centerline_metrics,
        "angular": {"shell": angular["shell"], **angular_metrics},
        "plane": plane_metrics,
    }
    path = results / "comparison.yaml"
    path.write_text(
        "# GENERATED by Allcompare.py -- do not edit. Errors are fractions.\n"
        f"# {verdict['status'].upper()}\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    written.append(path)

    if not args.no_plots:
        grid = comparison.analytical_meridional_grid(cfg, geom, exit_state)
        written += figures.plot_centerline_pressure(results, profile, label, banner)
        written += figures.plot_angular(results, angular, label, banner)
        written += figures.plot_plane(results, plane, grid, label, banner)

    print()
    for item in written:
        print(f"    wrote {item}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
