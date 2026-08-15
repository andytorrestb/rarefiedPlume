#!/usr/bin/env python3
"""Assemble the matrix: overlap, sweep table, contact sheets, curves.

    ./analyse.py                     everything
    ./analyse.py --no-sheets         skip the ffmpeg steps
    ./analyse.py --field rhoN        dissect a different field's series

Runs after every case has been post-processed and rendered. Produces, under
``results/``:

    overlap.yaml            long runs against short ones, where they overlap
    sweep-table.csv         one row per case: occupancy, budget, error
    contact-final.png       the matrix at each case's own final frame
    contact-common.png      the matrix at one time every case reached
    convergence.yaml        image-difference curves, one per case

The four artefacts answer different questions and the two sheets in particular
are easy to confuse:

**contact-final** is the headline. Each cell is the answer that case actually
delivers -- its own longest average -- so the columns differ because the runs
differ in how long they sampled.

**contact-common** is the control. At a time every case reached, the three
cases of a weight row are *the same run truncated*, so their cells should be
pixel-identical. A row that varies in that sheet is a determinism failure you
can see, and ``overlap.yaml`` is the same statement measured on the fields.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

try:
    from plumetools.cai2012 import audit, health
    from plumetools.cai2012.post import PostError
    from plumetools.viz import dissect
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test,cases]"')


def load_manifest() -> dict:
    path = HERE / "manifest.yaml"
    if not path.is_file():
        sys.exit("analyse: no manifest.yaml. Generate the cases first:\n"
                 "    ./generate_cases.py")
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def case_rows(manifest: dict, study) -> list:
    """Manifest entries, with their grid position attached."""
    weights = [w.name for w in study.enabled_weights()]
    samplings = [s.name for s in study.enabled_samplings()]
    rows = []
    for entry in manifest.get("cases", []):
        try:
            rows.append({
                **entry,
                "row": weights.index(entry["weight_level"]),
                "column": samplings.index(entry["sampling_level"]),
            })
        except ValueError:
            # A case whose level is no longer in study.yaml. Reported rather
            # than dropped: a matrix quietly missing a row still draws.
            print(f"    {entry['name']}: not in the current study axes; "
                  f"skipped")
    return rows


# --------------------------------------------------------------------------- #
# 1. the overlap
# --------------------------------------------------------------------------- #

def check_overlap(rows: list) -> dict:
    """Compare every longer case in a weight row against the shorter ones.

    Within a row the cases differ **only** in ``endTime``: same mesh, same
    particle weight, same time step, same decomposition, same seeding. The
    short run is therefore a prefix of the long one, and their running means
    must be the same numbers at any time both reached.
    """
    by_row = {}
    for entry in rows:
        by_row.setdefault(entry["weight_level"], []).append(entry)

    results = []
    for weight, entries in sorted(by_row.items()):
        ordered = sorted(entries, key=lambda e: e["sampling_domain_transits"])
        longest = ordered[-1]
        for shorter in ordered[:-1]:
            long_dir, short_dir = HERE / longest["path"], HERE / shorter["path"]
            if not (long_dir.is_dir() and short_dir.is_dir()):
                continue
            try:
                result = audit.compare_overlap(long_dir, short_dir)
            except PostError as exc:
                print(f"    {weight}: {exc}")
                continue
            result["weight_level"] = weight
            result["long_case"] = longest["name"]
            result["short_case"] = shorter["name"]
            results.append(result)
            print("\n".join("    " + line
                            for line in audit.overlap_report(result)))
            print()

    disagreed = [r for r in results if not r["all_identical"]]
    untested = [r for r in results if r["n_compared"] == 0]
    return {
        "n_pairs": len(results),
        "n_disagreed": len(disagreed),
        "n_untested": len(untested),
        "deterministic": bool(results) and not disagreed and not untested,
        "pairs": results,
    }


# --------------------------------------------------------------------------- #
# 2. the sweep table
# --------------------------------------------------------------------------- #

def build_table(rows: list, manifest: dict) -> list:
    """One row per case, from its metrics and its occupancy audit."""
    shared = manifest.get("held_fixed") or {}
    delta_t = shared.get("deltaT_s")
    table = []
    for entry in rows:
        case_dir = HERE / entry["path"]
        metrics = _read(case_dir / "results" / "metrics.yaml")
        occupancy = _read(case_dir / "results" / "occupancy.yaml")
        transit = None
        convergence = _read(case_dir / "results" / "convergence.yaml")
        if convergence:
            transit = convergence.get("domain_transit_s")
        steps_per_transit = (float(transit) / float(delta_t)
                             if transit and delta_t else 0.0)
        table.append(audit.sweep_row(entry, metrics, occupancy,
                                     steps_per_transit))
    return table


def _read(path: Path):
    if not Path(path).is_file():
        return None
    return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}


def print_table(table: list) -> None:
    columns = ("case", "particles_per_cell", "sampling_transits",
               "exit_occupancy_estimated", "exit_occupancy_measured",
               "plume_median_occupancy", "statistical_budget",
               "density_mean_rel_error_percent")
    heads = {"case": "case", "particles_per_cell": "ppc",
             "sampling_transits": "transits",
             "exit_occupancy_estimated": "exit est",
             "exit_occupancy_measured": "exit meas",
             "plume_median_occupancy": "plume med",
             "statistical_budget": "budget",
             "density_mean_rel_error_percent": "density mean"}
    widths = {c: max(len(heads[c]), 12) for c in columns}
    widths["case"] = 16
    print("    " + "  ".join(f"{heads[c]:>{widths[c]}}" for c in columns))
    for row in table:
        cells = []
        for column in columns:
            value = row.get(column)
            if value is None:
                cells.append(f"{'':>{widths[column]}}")
            elif isinstance(value, float):
                cells.append(f"{value:>{widths[column]}.4g}")
            else:
                cells.append(f"{str(value):>{widths[column]}}")
        print("    " + "  ".join(cells))


# --------------------------------------------------------------------------- #
# 3. the contact sheets
# --------------------------------------------------------------------------- #

def frame_at(viz_dir: Path, field: str, time: float | None):
    """``(filename, time)`` for one frame of a series.

    ``time = None`` takes the last frame. An explicit time takes the newest
    frame at or before it, so a common-time sheet lines up even when a case's
    write schedule put its frames a rounding apart.
    """
    frames, times = dissect.series_frames(viz_dir, field)
    if not frames:
        return None, None
    if time is None:
        return frames[-1], times[-1]
    candidates = [(t, f) for f, t in zip(frames, times) if t <= time * (1 + 1e-9)]
    if not candidates:
        return None, None
    chosen = max(candidates)
    return chosen[1], chosen[0]


def build_sheet(rows: list, field: str, output: Path, *, common_time=None,
                title: str, notes: list, ffmpeg: str = "ffmpeg") -> bool:
    """Assemble one contact sheet. Returns whether it was written."""
    cells = []
    for entry in rows:
        viz_dir = HERE / entry["path"] / "results" / "viz"
        if not (viz_dir / "manifest.yaml").is_file():
            continue
        try:
            name, time = frame_at(viz_dir, field, common_time)
        except dissect.DissectError:
            continue
        if name is None:
            continue
        directory = dissect.series_directory(viz_dir, field)
        cells.append(dissect.Cell(
            row=entry["row"], column=entry["column"],
            path=directory / name,
            label=f"{entry['name']}  t={time:.4g}"))

    if not cells:
        print(f"    no rendered {field} frames yet; no sheet")
        return False

    font = dissect.find_font()
    try:
        dissect.build_contact_sheet(
            cells, output, font=font, cell_width=300, cell_height=700,
            trim_bottom=110, ffmpeg=ffmpeg)
    except dissect.DissectError as exc:
        print(f"    contact sheet FAILED: {exc}")
        return False
    dissect.write_legend(
        output.with_suffix(".yaml"), cells, title=title,
        rows=sorted({e["weight_level"] for e in rows}),
        columns=sorted({e["sampling_level"] for e in rows}), notes=notes)
    print(f"    wrote {output}  ({len(cells)} cell(s))")
    return True


def common_frame_time(rows: list, field: str) -> float | None:
    """The latest time every rendered case has a frame at.

    The shortest case's last frame, in practice. It is found rather than
    assumed so the sheet still assembles while the matrix is part way through.
    """
    latest = []
    for entry in rows:
        viz_dir = HERE / entry["path"] / "results" / "viz"
        if not (viz_dir / "manifest.yaml").is_file():
            continue
        try:
            _, times = dissect.series_frames(viz_dir, field)
        except dissect.DissectError:
            continue
        if times:
            latest.append(max(times))
    return min(latest) if latest else None


# --------------------------------------------------------------------------- #
# 4. the convergence curves
# --------------------------------------------------------------------------- #

def build_curves(rows: list, field: str, ffmpeg: str = "ffmpeg") -> list:
    """An image-difference curve per case, beside its physical error."""
    curves = []
    for entry in rows:
        viz_dir = HERE / entry["path"] / "results" / "viz"
        if not (viz_dir / "manifest.yaml").is_file():
            continue
        try:
            frames, times = dissect.series_frames(viz_dir, field)
            if len(frames) < 2:
                continue
            curve = dissect.convergence_curve(
                dissect.series_directory(viz_dir, field), frames, times=times,
                ffmpeg=ffmpeg)
        except dissect.DissectError as exc:
            print(f"    {entry['name']}: {exc}")
            continue
        curve["case"] = entry["name"]
        curve["target_particles_per_cell"] = entry["target_particles_per_cell"]
        curve["sampling_domain_transits"] = entry["sampling_domain_transits"]

        physical = _read(HERE / entry["path"] / "results" / "metrics.yaml")
        print("\n".join("    " + line for line in dissect.convergence_report(
            curve, (physical or {}).get("centerline"))))
        print()
        curves.append(curve)
    return curves


# --------------------------------------------------------------------------- #

def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--field", default="dsmcRhoN",
                        help="which rendered series to dissect "
                             "(default: dsmcRhoN, the parcel density)")
    parser.add_argument("--no-sheets", action="store_true",
                        help="skip everything that needs ffmpeg")
    parser.add_argument("--no-overlap", action="store_true")
    parser.add_argument("--ffmpeg", default="ffmpeg")
    args = parser.parse_args(argv[1:])

    manifest = load_manifest()
    study = health.load_study(HERE / "study.yaml")
    rows = case_rows(manifest, study)
    if not rows:
        print("analyse: the manifest lists no case in the current axes")
        return 1

    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)

    print("### overlap: is the pipeline deterministic?")
    overlap = {"deterministic": None}
    if not args.no_overlap:
        overlap = check_overlap(rows)
        audit.write_yaml(results / "overlap.yaml", overlap,
                         "Long runs against short ones where they overlap.")
        print(f"    wrote {results / 'overlap.yaml'}")
    print()

    print("### sweep table")
    table = build_table(rows, manifest)
    print_table(table)
    written = audit.write_sweep_table(results / "sweep-table.csv", table)
    print(f"    wrote {written}")
    print()

    if not args.no_sheets:
        if not dissect.have_ffmpeg(args.ffmpeg):
            print(f"### contact sheets skipped: {args.ffmpeg} is not on PATH.")
            print("    Everything above is already on disk.")
            return 0

        print("### contact sheets")
        pinned = ("Colour ranges are pinned in viz.yaml, so every cell is on "
                  "the same scale. Without that each case autoscales to its "
                  "own data and no two cells are comparable.")
        build_sheet(rows, args.field, results / "contact-final.png",
                    title=f"{args.field}: each case at its own final frame",
                    notes=[pinned,
                           "Rows are parcels per cell, columns are averaging "
                           "duration. Each cell is the answer that case "
                           "delivers."],
                    ffmpeg=args.ffmpeg)

        common = common_frame_time(rows, args.field)
        if common is not None:
            build_sheet(rows, args.field, results / "contact-common.png",
                        common_time=common,
                        title=f"{args.field}: every case at t = {common:.6g} s",
                        notes=[pinned,
                               "At a time every case reached, the three cases "
                               "of a weight row are THE SAME RUN TRUNCATED and "
                               "their cells should be pixel-identical. A row "
                               "that varies here is a determinism failure. "
                               "See results/overlap.yaml for the same "
                               "statement measured on the fields.",
                               f"common time {common:.6g} s"],
                        ffmpeg=args.ffmpeg)
        print()

        print("### convergence curves")
        curves = build_curves(rows, args.field, ffmpeg=args.ffmpeg)
        if curves:
            path = dissect.write_curve(results / "convergence.yaml", curves)
            print(f"    wrote {path}  ({len(curves)} series)")
        print()

    if overlap.get("deterministic") is False:
        print("analyse: THE OVERLAP CHECK FAILED. Cases in a weight row differ "
              "only in\n         endTime, so agreeing was not optional. The "
              "sweep attributes its\n         differences to statistics, and "
              "this says part of them are not.\n"
              "         See results/overlap.yaml.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
