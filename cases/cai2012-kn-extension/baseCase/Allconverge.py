#!/usr/bin/env python3
"""Was the sampling window steady where it matters? Frame by frame.

    ./Allconverge.py                   every sampled frame
    ./Allconverge.py --stride 2        every other frame
    ./Allconverge.py --tolerance 0.01  a tighter stationarity criterion

Produces, under results/:

    stationarity.yaml   per region and per frame: the running-mean n/n0 and
                        the instantaneous n/n0, and the verdict
    stationarity.png    the same, drawn (needs matplotlib)

A flat total parcel count is NOT evidence that the far field has converged: it
is dominated by the near field and plateaus while the far centreline is still
filling (docs/cai2012-case.md 5.1). So the regions here are the centreline tube
in bands of X/D and the r/D = 9-10 shell in polar sectors; each must show a
running mean that has stopped moving AND an instantaneous field with no trend
across the window. The domain total is printed as context and is not part of
the verdict. See plumetools/cai2012ext/convergence.py.

Allaudit.py (from cases/cai2012-health, unchanged) is the companion: measured
parcels per cell, and the centreline error against Cai at every frame.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    import yaml

    from plumetools.cai2012 import post
    from plumetools.cai2012ext import convergence, figures, interface_case
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stride", type=int, default=1,
                        help="read every Nth sampled frame (default: all)")
    parser.add_argument("--tolerance", type=float,
                        default=convergence.DEFAULT_TOLERANCE,
                        help="relative change across the window allowed "
                             f"(default {convergence.DEFAULT_TOLERANCE:g})")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv[1:])

    # Direct exit or interface mode: the same run settings Allmesh wrote.
    context = interface_case.post_context(HERE)
    cfg, geom, exit_state, run = (context.cfg, context.geom, context.exit_state,
                                  context.run)
    verdict = context.verdict
    banner = figures.status_banner(verdict["production_ready"])
    label = cfg.meta.get("case_name") or HERE.name

    times = convergence.sampled_frames(HERE, run.average_start_s)
    print(f"### {label}")
    if not times:
        print(f"    no sampled frame at or after t = {run.average_start_s:.4e} s")
        return 0
    mesh_time = convergence.geometry_time(HERE, times)
    if mesh_time is None:
        print("    no cell centres in this case; run ./Allpost first",
              file=sys.stderr)
        return 1
    centres = post.read_internal_field(HERE / mesh_time / "C")
    volume_path = HERE / mesh_time / "V"
    if not volume_path.is_file():
        print(f"    no cell volumes at {volume_path}; the regional means are "
              f"volume-weighted. Run ./Allpost first.", file=sys.stderr)
        return 1
    volumes = post.read_internal_field(volume_path)
    masks = convergence.region_masks(centres, cfg, geom)
    n0 = float(exit_state.number_density_per_m3)

    selected = times[::max(1, args.stride)]
    if selected[-1] != times[-1]:
        selected.append(times[-1])
    print(f"    sampled frames       {len(times)} ({times[0]} .. {times[-1]}), "
          f"reading {len(selected)}")
    print(f"    averaging from       {run.average_start_s:.6e} s")
    print()

    frames = []
    for time in selected:
        mean = post.read_internal_field(HERE / time / "rhoNMean")
        if mean.size != volumes.size:
            print(f"    {time}: rhoNMean has {mean.size:,} values against "
                  f"{volumes.size:,} cells; skipped", file=sys.stderr)
            continue
        frame = {"time": float(time),
                 "averaged_transits": (float(time) - run.average_start_s)
                 / run.domain_transit_s,
                 "mean": convergence.regional_densities(mean, volumes, masks, n0)}
        instant_path = HERE / time / "rhoN"
        if instant_path.is_file():
            instant = post.read_internal_field(instant_path)
            if instant.size == volumes.size:
                frame["instant"] = convergence.regional_densities(
                    instant, volumes, masks, n0)
        frames.append(frame)
        print(f"    {time:>12}  {frame['averaged_transits']:6.3f} transits")

    assessment = convergence.assess(frames, tolerance=args.tolerance)
    print()
    print(f"  {'region':<36} {'drift':>8} {'trend':>8} {'+-2se':>8}  status")
    for name, region in assessment.get("regions", {}).items():
        print(f"  {name:<36} {100 * region['running_mean_drift']:>7.2f}% "
              f"{100 * region['instant_relative_change']:>7.2f}% "
              f"{200 * region['instant_trend_standard_error']:>7.2f}%  "
              + region["status"].upper()
              + ("" if region["in_verdict"] else "  (context only)"))
    print()
    print(f"  verdict: {assessment['verdict']}")

    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    document = {
        "case": label,
        "production_status": verdict["status"],
        "average_start_s": float(run.average_start_s),
        "domain_transit_s": float(run.domain_transit_s),
        "assessment": assessment,
        "frames": frames,
    }
    path = results / "stationarity.yaml"
    path.write_text(
        "# GENERATED by Allconverge.py -- do not edit.\n"
        f"# {verdict['status'].upper()}: {assessment['verdict']}\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    print(f"    wrote {path}")
    if not args.no_plots:
        for item in figures.plot_stationarity(
                results, frames, assessment, run.domain_transit_s,
                run.average_start_s, label, banner):
            print(f"    wrote {item}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
