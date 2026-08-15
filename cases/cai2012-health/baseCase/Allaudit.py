#!/usr/bin/env python3
"""Audit this case's statistical health, frame by frame.

    ./Allaudit.py                  every sampled frame, plus the occupancy audit
    ./Allaudit.py --latest-only    the occupancy audit alone (much quicker)
    ./Allaudit.py --stride 2       every other frame

Produces, under ``results/``:

    occupancy.yaml     measured parcels per cell, per region, at the last frame
    convergence.yaml   centreline error against Cai at EVERY sampled frame

This is the half of ``cases/cai2012-health`` that needs no ParaView and no
solver -- only fields the run already wrote. ``./Allpost`` produces the four
validation outputs at the latest time and is unchanged from ``cases/cai2012``;
this adds the two things that family does not measure.

**Occupancy.** ``plumetools.cai2012.checks`` prints an *estimate* of the exit
cell's parcel count and says only a post-run audit of the sampled ``dsmcRhoN``
field can say what it was. This is that audit: ``dsmcRhoNMean * V`` is the mean
parcel count per cell over the sampling window, exactly, and it is reported per
region because the target is met at the exit cell and nowhere else by
construction.

**Convergence, physically.** The running mean is written many times across the
sampling window, so every frame is a shorter average of the same run. Scoring
each one against Cai's collisionless solution gives the error as a function of
how long the average has run -- the physical curve that
``plumetools.viz.dissect``'s image metric has to be read beside.

Cell centres and volumes are read **once**, from whichever time
``postProcess -func writeCellCentres`` wrote them into, and every frame is read
against them. They describe the mesh, the mesh does not move, and writing 60 MB
of them into each of 37 frames would cost more disk than the fields do.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012 import audit, inflow, mesh, post
    from plumetools.cai2012.config import load_case_config
    from plumetools.cai2012.geometry import from_config
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def sampled_times(case_dir: Path, average_start_s: float) -> list[str]:
    """Time directories at or after ``average_start_s``, ascending.

    The only frames that carry a ``*Mean`` field. Everything earlier was written
    during the transient, when ``fieldAverage`` had not started -- reading one
    would silently score an instantaneous field against an analytical curve and
    measure shot noise.
    """
    times = []
    for entry in Path(case_dir).iterdir():
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        # A tolerance of one part in 1e9 rather than an exact >=: the write
        # schedule is counted in steps and the boundary frame can land a
        # rounding below the start.
        if value >= average_start_s * (1.0 - 1.0e-9) and value > 0.0:
            if (entry / "rhoNMean").is_file():
                times.append((value, entry.name))
    return [name for _, name in sorted(times)]


def geometry_time(case_dir: Path, times: list[str]) -> str | None:
    """The time directory holding ``C``, or ``None`` if none does."""
    for name in reversed(times):
        if (Path(case_dir) / name / "C").is_file():
            return name
    for entry in sorted(Path(case_dir).iterdir()):
        if entry.is_dir() and (entry / "C").is_file():
            return entry.name
    return None


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--latest-only", action="store_true",
                        help="occupancy audit only; skip the per-frame curve")
    parser.add_argument("--stride", type=int, default=1,
                        help="audit every Nth frame (default: all of them)")
    args = parser.parse_args(argv[1:])

    cfg = load_case_config(HERE)
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom, min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)

    label = cfg.meta.get("case_name") or HERE.name
    print(f"### {label}")

    times = sampled_times(HERE, run.average_start_s)
    if not times:
        print(f"    no sampled frame at or after t = {run.average_start_s:.4e} s. "
              f"This case has\n    not reached its sampling window yet, or was "
              f"run without fieldAverage.")
        return 0

    mesh_time = geometry_time(HERE, times)
    if mesh_time is None:
        print("    no cell centres anywhere in this case. They come from\n"
              "        postProcess -func writeCellCentres -latestTime\n"
              "    which ./Allpost runs. Run it first.", file=sys.stderr)
        return 1

    print(f"    sampled frames       {len(times)} "
          f"({times[0]} .. {times[-1]})")
    print(f"    geometry from        {mesh_time}")
    print(f"    averaging starts at  {run.average_start_s:.6e} s")
    print()

    results = HERE / "results"

    # ----------------------------------------------------------------- #
    # the occupancy audit, at the final frame
    # ----------------------------------------------------------------- #
    final = times[-1]
    sampled = post.read_case(HERE, cfg.gas.mass_kg, time=final,
                             mesh_time=mesh_time)
    parcel_path = HERE / final / "dsmcRhoNMean"
    if not parcel_path.is_file():
        print(f"    no dsmcRhoNMean at {parcel_path}; the occupancy audit "
              f"needs the\n    PARCEL number density, not just rhoNMean.",
              file=sys.stderr)
        return 1

    averaged_steps = (float(final) - run.average_start_s) / run.delta_t_s
    occupancy = audit.occupancy_audit(
        sampled, cfg, geom, exit_state,
        n_equivalent_particles=run.n_equivalent_particles,
        parcel_density=post.read_internal_field(parcel_path),
        averaged_steps=averaged_steps,
        core_cell_size_m=plan.core_cell_size_m)
    occupancy["case"] = label
    occupancy["time"] = final
    occupancy["estimated_exit_particles_per_cell"] = float(
        run.exit_particles_per_cell)

    print("\n".join(audit.occupancy_report(occupancy)))
    print()
    measured = [r for r in occupancy["regions"] if r["name"] == "exit"][0]
    print(f"  checks.py ESTIMATED the exit cell at "
          f"{run.exit_particles_per_cell:.1f} parcels/cell;")
    print(f"  the sampled field MEASURES {measured['median']:.1f} "
          f"(median over {measured['n_cells']} exit cells).")
    print()
    audit.write_yaml(results / "occupancy.yaml", occupancy,
                     f"Measured parcel occupancy for {label} at t = {final}.")
    print(f"    wrote {results / 'occupancy.yaml'}")

    if args.latest_only:
        return 0

    # ----------------------------------------------------------------- #
    # the physical convergence curve, frame by frame
    # ----------------------------------------------------------------- #
    print()
    print("Centreline error against Cai, per frame (the running mean settling)")
    print(f"  {'frame':>5} {'time':>12} {'transits':>9} {'density max':>12} "
          f"{'mean':>9} {'RMS':>9}")

    frames = []
    for index, time in enumerate(times[::max(1, args.stride)]):
        try:
            snapshot = post.read_case(HERE, cfg.gas.mass_kg, time=time,
                                      mesh_time=mesh_time)
        except post.PostError as exc:
            print(f"  {index:>5} {time:>12} skipped: {exc}")
            continue
        profile = post.centerline(snapshot, cfg, geom, exit_state)
        metrics = post.centerline_metrics(profile)
        elapsed = (float(time) - run.average_start_s) / run.domain_transit_s
        frames.append({
            "frame": index,
            "time": float(time),
            "averaged_transits": float(elapsed),
            "averaged_steps": float((float(time) - run.average_start_s)
                                    / run.delta_t_s),
            **{k: float(v) for k, v in metrics.items()},
        })
        print(f"  {index:>5} {float(time):>12.6g} {elapsed:>9.3f} "
              f"{100 * metrics['density_max_rel_error']:>11.2f}% "
              f"{100 * metrics['density_mean_rel_error']:>8.2f}% "
              f"{100 * metrics['density_rms_rel_error']:>8.2f}%")

    document = {
        "case": label,
        "target_particles_per_cell": float(cfg.resolution.target_particles_per_cell),
        "sampling_domain_transits": float(cfg.dsmc.sampling_domain_transits),
        "average_start_s": float(run.average_start_s),
        "domain_transit_s": float(run.domain_transit_s),
        "delta_t_s": float(run.delta_t_s),
        "n_equivalent_particles": float(run.n_equivalent_particles),
        "note": "One row per written frame of the running average. THE "
                "PHYSICAL curve -- read plumetools.viz.dissect's image metric "
                "beside this, never instead of it.",
        "frames": frames,
    }
    audit.write_yaml(results / "convergence.yaml", document,
                     f"Centreline error against Cai per frame for {label}.")
    print()
    print(f"    wrote {results / 'convergence.yaml'}  ({len(frames)} frame(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
