r"""Dissect a rendered sweep: a contact sheet, and a convergence curve per case.

Runs after the rendering, over the PNGs already on disk. **No ParaView**, the
same contract as :mod:`plumetools.viz.video`: importable, testable without
ffmpeg installed, runnable with plain ``python``::

    python plumetools/viz/dissect.py sheet   --matrix cases/cai2012-health/matrix.yaml
    python plumetools/viz/dissect.py converge cases/cai2012-health/Cases/*/*/results/viz

Two things, both built with ffmpeg filters and no new dependency.

**1. The contact sheet.** Rows are particle weight, columns are averaging
duration, one cell per case -- the whole sweep as a single image.
``scale``/``pad`` normalise the cells (a slice's pixel width follows the
geometry, so it is not the same in every case), ``drawtext`` labels them, and
``xstack`` assembles the grid at explicit pixel offsets.

**2. The convergence curve.** Successive frames of a running average stop
differing once the average has converged, and ``ffmpeg -lavfi psnr`` measures
by how much, straight from the image data.

What an image metric is, and is not
-----------------------------------
PSNR between two log-scaled, colour-mapped, palette-quantised PNGs is **a
perceptual proxy for "has this stopped changing"**. It is not a physical error
and it is not a convergence criterion:

* the colour map is **logarithmic**, so a fixed PSNR means a different relative
  change at the plume edge than on the axis;
* it is **clipped** to the pinned range, so anything outside the range moves
  without registering at all;
* the mapping to 8-bit RGB **quantises**, so changes below a colour step are
  invisible, and it saturates -- two converged frames give ``psnr = inf``,
  which is a statement about the palette, not about the flow;
* it weights every pixel equally, so the vacuum -- most of the picture --
  counts as much as the plume core.

So it is always reported **beside** the physical error from
:func:`plumetools.cai2012.post.centerline_metrics`, never on its own, and
:func:`convergence_report` says so on every run. Where the two disagree, the
disagreement is the finding: an image that has stopped changing while the
centreline error is still falling means the change has moved below a colour
step, and the picture has stopped being able to show it.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import yaml

from plumetools.viz.video import CONCAT_NAME, VideoError, concat_list, have_ffmpeg

#: Fonts ``drawtext`` is tried against, in order.
#:
#: A missing font is not fatal -- the sheet is drawn unlabelled and the legend
#: file still says which cell is which -- but an unlabelled 3 x 3 of nearly
#: identical plumes is nearly useless, so it is worth looking in more than one
#: place.
FONT_CANDIDATES = (
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    "/Library/Fonts/Arial.ttf",
    "C:/Windows/Fonts/arial.ttf",
)

#: One line of ``psnr``/``ssim``'s ``stats_file``: ``n:1 mse_avg:12.3 ...``.
#:
#: Parsed generically rather than against a fixed key list: the keys depend on
#: the pixel format ffmpeg chose (``mse_y/u/v`` for YUV input, ``mse_r/g/b`` for
#: RGB), and a parser that knew only one of them would return nothing at all for
#: the other -- an empty curve, which looks like a converged one.
_STAT = re.compile(r"(\w+):([-+0-9.eE]+|inf|nan)")


class DissectError(RuntimeError):
    """ffmpeg failed, or the sweep cannot be assembled."""


# --------------------------------------------------------------------------- #
# 1. the contact sheet
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Cell:
    """One image in the grid.

    Attributes:
        row: 0-based row index -- the particle-weight axis.
        column: 0-based column index -- the averaging-duration axis.
        path: the PNG.
        label: drawn in the cell's header strip.
    """

    row: int
    column: int
    path: Path
    label: str = ""


def grid_layout(cells: Sequence[Cell], cell_width: int,
                cell_height: int) -> str:
    """``xstack``'s ``layout``, in explicit pixels.

    Args:
        cells: the grid, in the order they will be given to ffmpeg.
        cell_width, cell_height: the size every cell is scaled to.

    Returns:
        ``0_0|300_0|600_0|0_740|...``

    Explicit pixels rather than ``w0``/``h0`` arithmetic: every cell is padded
    to the same size first, so the offsets are known exactly, and the symbolic
    form silently produces overlapping cells when an input turns out to be a
    different size than assumed.
    """
    return "|".join(f"{cell.column * cell_width}_{cell.row * cell_height}"
                    for cell in cells)


def find_font(candidates: Iterable[str] = FONT_CANDIDATES) -> str | None:
    """The first usable ``drawtext`` font, or ``None``."""
    for candidate in candidates:
        if Path(candidate).is_file():
            return candidate
    return None


def escape_drawtext(text: str) -> str:
    """Escape a label for ``drawtext``'s ``text=`` option.

    ``:`` separates filter options and ``'`` quotes them, so an unescaped label
    like ``ppc020/s1p5 t:0.0099`` truncates the filter graph or fails it -- and
    a label is exactly the sort of string that grows a colon later.
    """
    out = text
    for character in ("\\", ":", "'", "%", "[", "]", ","):
        out = out.replace(character, "\\" + character)
    return out


def cell_filter(index: int, cell: Cell, *, cell_width: int, cell_height: int,
                header: int, font: str | None, font_size: int,
                trim_bottom: int, background: str) -> str:
    """The filter chain that normalises and labels one cell.

    Args:
        index: the input index in the ffmpeg command.
        cell: the cell.
        cell_width, cell_height: the padded size, header included.
        header: pixels reserved at the top for the label.
        font: a ``drawtext`` font file, or ``None`` to skip labelling.
        font_size: label size in pixels.
        trim_bottom: pixels to crop off the bottom of the source image before
            scaling. The renderer puts a colour bar there, and with the range
            pinned across the matrix all nine bars are identical -- nine copies
            of one legend, taking a fifth of the sheet.
        background: pad colour.

    Returns:
        One ``[in]...[out]`` chain, ready to join with ``;``.
    """
    steps = []
    if trim_bottom > 0:
        # max(...,1) rather than a bare subtraction: a cell shorter than the
        # trim would otherwise ask for a zero-height crop, which ffmpeg rejects
        # with a message about the filter rather than about the trim.
        steps.append(f"crop=iw:max(ih-{int(trim_bottom)}\\,1):0:0")
    inner = max(1, cell_height - header)
    steps.append(f"scale={cell_width}:{inner}:force_original_aspect_ratio=decrease")
    steps.append(f"pad={cell_width}:{cell_height}:(ow-iw)/2:{header}:"
                 f"color={background}")
    if font and cell.label:
        steps.append(
            f"drawtext=fontfile={font}:text={escape_drawtext(cell.label)}"
            f":x=8:y={max(0, (header - font_size) // 2)}"
            f":fontsize={font_size}:fontcolor=black")
    return f"[{index}:v]" + ",".join(steps) + f"[c{index}]"


def contact_sheet_command(cells: Sequence[Cell], output: str, *,
                          cell_width: int = 300, cell_height: int = 780,
                          header: int = 40, font: str | None = None,
                          font_size: int = 22, trim_bottom: int = 0,
                          background: str = "white",
                          ffmpeg: str = "ffmpeg") -> list[str]:
    """The argv that assembles the contact sheet.

    Pure, so the layout arithmetic can be checked without ffmpeg installed --
    and a grid whose offsets are wrong produces an image, just not the one
    anybody meant.

    Args:
        cells: the grid.
        output: the PNG to write.
        cell_width, cell_height: size of one cell, header included.
        header: label strip height.
        font: ``drawtext`` font file, or ``None`` for an unlabelled sheet.
        font_size: label size.
        trim_bottom: pixels cropped off each source image before scaling.
        background: pad and fill colour.
        ffmpeg: the binary.

    Returns:
        The argument vector.

    Raises:
        DissectError: if there are no cells.
    """
    if not cells:
        raise DissectError("a contact sheet of no cases is not an image")

    command = [ffmpeg, "-y", "-hide_banner", "-loglevel", "error"]
    for cell in cells:
        command += ["-i", str(cell.path)]

    chains = [cell_filter(i, cell, cell_width=cell_width,
                          cell_height=cell_height, header=header, font=font,
                          font_size=font_size, trim_bottom=trim_bottom,
                          background=background)
              for i, cell in enumerate(cells)]
    stacked = "".join(f"[c{i}]" for i in range(len(cells)))
    if len(cells) == 1:
        # xstack refuses inputs=1, and a one-case sheet is a legitimate thing
        # to ask for while a matrix is part way through.
        graph = ";".join(chains) + ";[c0]null[out]"
    else:
        graph = ";".join(chains) + (
            f";{stacked}xstack=inputs={len(cells)}"
            f":layout={grid_layout(cells, cell_width, cell_height)}"
            f":fill={background}[out]")

    command += ["-filter_complex", graph, "-map", "[out]", "-frames:v", "1",
                output]
    return command


def build_contact_sheet(cells: Sequence[Cell], output: Path, **options) -> Path:
    """Assemble and write the contact sheet.

    Args:
        cells: the grid.
        output: the PNG to write.
        **options: as :func:`contact_sheet_command`.

    Returns:
        The path written.

    Raises:
        DissectError: if ffmpeg is missing, fails, or writes nothing.
    """
    ffmpeg = options.get("ffmpeg", "ffmpeg")
    if not have_ffmpeg(ffmpeg):
        raise DissectError(f"{ffmpeg} is not on PATH")
    missing = [str(cell.path) for cell in cells if not Path(cell.path).is_file()]
    if missing:
        raise DissectError(
            f"{len(missing)} cell image(s) do not exist, so the grid would be "
            f"missing cases without saying which: {missing[:3]}")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        contact_sheet_command(cells, str(output), **options),
        capture_output=True, text=True)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        raise DissectError(f"ffmpeg failed assembling the contact sheet "
                           f"(exit {result.returncode}): "
                           + (detail[0] if detail else "no output"))
    if not output.is_file() or output.stat().st_size == 0:
        raise DissectError("ffmpeg wrote no contact sheet")
    return output


def write_legend(path: Path, cells: Sequence[Cell], *, title: str,
                 rows: Sequence[str], columns: Sequence[str],
                 notes: Sequence[str] = ()) -> Path:
    """Record which grid position is which case.

    A contact sheet is an image, and an image records nothing about where its
    cells came from. It is also the one artefact of a study most likely to be
    pasted somewhere on its own, so what it shows has to be written down
    beside it -- including when ``drawtext`` was unavailable and the sheet came
    out unlabelled.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "title": title,
        "rows": list(rows),
        "columns": list(columns),
        "notes": list(notes),
        "cells": [
            {"row": cell.row, "column": cell.column, "label": cell.label,
             "image": str(cell.path).replace("\\", "/")}
            for cell in cells
        ],
    }
    path.write_text(
        "# GENERATED by plumetools.viz.dissect -- do not edit.\n"
        "# What each cell of the contact sheet beside this file is.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


# --------------------------------------------------------------------------- #
# 2. the convergence curve
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class FrameMetric:
    """One frame compared against another.

    Attributes:
        frame: 1-based index into the compared sequence.
        reference: what it was compared with -- ``"previous"`` or ``"final"``.
        value: PSNR in dB, or the SSIM ``All`` figure.
        mse: mean squared error, when the metric reports one.
        saturated: the two frames were identical, so the metric is infinite.
            **Not** a convergence result: it means the change fell below one
            colour step of an 8-bit palette.
    """

    frame: int
    reference: str
    value: float
    mse: float = float("nan")
    saturated: bool = False

    def as_dict(self) -> dict:
        return {"frame": self.frame, "reference": self.reference,
                "value": None if self.saturated else self.value,
                "mse": None if self.mse != self.mse else self.mse,
                "saturated": self.saturated}


def parse_stats(text: str, metric: str = "psnr") -> list[dict]:
    """Parse ffmpeg's ``stats_file`` output into one dict per frame.

    Args:
        text: the stats stream.
        metric: ``"psnr"`` or ``"ssim"``, which decides the summary key.

    Returns:
        One dict per frame, each with at least ``n`` and ``value``.

    The key names depend on the pixel format ffmpeg picked -- ``mse_y``/
    ``psnr_y`` for YUV, ``mse_r``/``psnr_r`` for RGB -- so every ``key:value``
    pair is read and the summary is looked up among the names either format
    could produce. A parser hard-coded to one of them returns an empty curve for
    the other, and an empty curve is indistinguishable from a converged one.
    """
    summary_keys = ("psnr_avg", "psnr_y", "psnr_r") if metric == "psnr" \
        else ("All", "all")
    frames = []
    for line in text.splitlines():
        pairs = dict(_STAT.findall(line))
        if "n" not in pairs:
            continue
        record = {}
        for key, raw in pairs.items():
            try:
                record[key] = float(raw)
            except ValueError:
                continue
        value = next((record[k] for k in summary_keys if k in record), None)
        if value is None:
            continue
        record["value"] = value
        record["mse"] = next(
            (record[k] for k in ("mse_avg", "mse_y", "mse_r") if k in record),
            float("nan"))
        frames.append(record)
    return frames


def frame_pairs(frames: Sequence[str], reference: str = "previous"
                ) -> tuple[list[str], list[str]]:
    """The two sequences to compare, frame for frame.

    Args:
        frames: file names in playing order.
        reference: ``"previous"`` compares each frame with the one before it --
            "has it stopped changing". ``"final"`` compares each with the last
            one -- "how far is it from where it ends up".

    Returns:
        ``(left, right)``, equal length.

    Raises:
        DissectError: for fewer than two frames, or an unknown reference.

    Both are needed and they answer different questions. Successive-frame
    differences go to zero for a *stalled* series as readily as a converged
    one; distance-to-final says how much of the answer the frame already had.
    Neither alone is a convergence statement.
    """
    frames = list(frames)
    if len(frames) < 2:
        raise DissectError(
            f"a convergence curve needs at least two frames, got {len(frames)}")
    if reference == "previous":
        return frames[:-1], frames[1:]
    if reference == "final":
        return frames[:-1], [frames[-1]] * (len(frames) - 1)
    raise DissectError(
        f"unknown reference {reference!r}; expected 'previous' or 'final'")


def metric_command(left_list: str, right_list: str, *, metric: str = "psnr",
                   ffmpeg: str = "ffmpeg") -> list[str]:
    """The argv that compares two concat lists frame by frame.

    Args:
        left_list, right_list: concat lists, relative to the working directory.
        metric: ``"psnr"`` or ``"ssim"``.
        ffmpeg: the binary.

    Returns:
        The argument vector.

    ``-r 1`` on both inputs is not cosmetic. The concat demuxer gives a list of
    stills no timestamps of their own, and without a rate ffmpeg emits
    ``non monotonically increasing dts`` for every frame and interleaves the
    warnings with the statistics being parsed off the same stream.

    ``stats_file=-`` sends the per-frame numbers to **stdout**, which is why
    the output goes to ``-f null -`` rather than a file: the measurement is the
    text, not a video.
    """
    if metric not in ("psnr", "ssim"):
        raise DissectError(f"unknown metric {metric!r}; expected psnr or ssim")
    return [
        ffmpeg, "-hide_banner", "-loglevel", "error",
        "-r", "1", "-f", "concat", "-safe", "0", "-i", left_list,
        "-r", "1", "-f", "concat", "-safe", "0", "-i", right_list,
        "-lavfi", f"{metric}=stats_file=-",
        "-f", "null", "-",
    ]


def measure_series(directory: Path, frames: Sequence[str], *,
                   reference: str = "previous", metric: str = "psnr",
                   ffmpeg: str = "ffmpeg") -> list[FrameMetric]:
    """Compare a folder of frames against itself, offset by one or against the last.

    Args:
        directory: where the frames are.
        frames: file names in playing order.
        reference: ``"previous"`` or ``"final"``.
        metric: ``"psnr"`` or ``"ssim"``.
        ffmpeg: the binary.

    Returns:
        One :class:`FrameMetric` per comparison.

    Raises:
        DissectError: if ffmpeg is missing or fails.
    """
    if not have_ffmpeg(ffmpeg):
        raise DissectError(f"{ffmpeg} is not on PATH")

    directory = Path(directory)
    left, right = frame_pairs(frames, reference)
    left_path = directory / f"{CONCAT_NAME}.a"
    right_path = directory / f"{CONCAT_NAME}.b"
    left_path.write_text(concat_list(left), encoding="utf-8", newline="\n")
    right_path.write_text(concat_list(right), encoding="utf-8", newline="\n")
    try:
        result = subprocess.run(
            metric_command(left_path.name, right_path.name, metric=metric,
                           ffmpeg=ffmpeg),
            cwd=directory, capture_output=True, text=True)
    finally:
        left_path.unlink(missing_ok=True)
        right_path.unlink(missing_ok=True)

    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        raise DissectError(f"ffmpeg {metric} failed (exit {result.returncode}): "
                           + (detail[0] if detail else "no output"))

    measured = []
    for record in parse_stats(result.stdout, metric):
        value = float(record["value"])
        measured.append(FrameMetric(
            frame=int(record["n"]),
            reference=reference,
            value=value,
            mse=float(record.get("mse", float("nan"))),
            # inf comes back from two byte-identical frames. Recording it as a
            # number would put a spike at the end of every curve and make a
            # palette limit look like perfect convergence.
            saturated=not (value < float("inf")),
        ))
    return measured


def convergence_curve(directory: Path, frames: Sequence[str], *,
                      times: Sequence[float] | None = None,
                      metric: str = "psnr",
                      ffmpeg: str = "ffmpeg") -> dict:
    """Both curves for one folder of frames, with their caveat attached.

    Args:
        directory: where the frames are.
        frames: file names in playing order.
        times: the solver time of each frame, from the render manifest.
        metric: ``"psnr"`` or ``"ssim"``.
        ffmpeg: the binary.

    Returns:
        A dict with ``previous``, ``final``, ``times`` and ``caveat``.

    The caveat travels **with the data**, not only in the report: this dict is
    what gets written to disk and read back by whoever plots it, and a column
    of PSNR values with no note attached is exactly how an image metric ends up
    being quoted as a physical one.
    """
    return {
        "directory": str(directory).replace("\\", "/"),
        "metric": metric,
        "n_frames": len(frames),
        "times": [float(t) for t in times] if times else None,
        "previous": [m.as_dict() for m in measure_series(
            directory, frames, reference="previous", metric=metric,
            ffmpeg=ffmpeg)],
        "final": [m.as_dict() for m in measure_series(
            directory, frames, reference="final", metric=metric,
            ffmpeg=ffmpeg)],
        "caveat":
            f"{metric.upper()} between log-scaled, colour-mapped, "
            f"palette-quantised PNGs. A PERCEPTUAL PROXY for 'has this stopped "
            f"changing', not a physical error: the scale is logarithmic and "
            f"clipped to the pinned range, the 8-bit palette hides changes "
            f"below one colour step, and the vacuum counts as many pixels as "
            f"the plume. Read it beside the centreline error from "
            f"plumetools.cai2012.post, never instead of it.",
    }


def convergence_report(curve: dict, physical: dict | None = None) -> list[str]:
    """The curve as printed lines, with the physical error beside it.

    Args:
        curve: what :func:`convergence_curve` returned.
        physical: the case's centreline metrics, or ``None``.

    Returns:
        Report lines.

    When ``physical`` is absent the report says so in place of the column,
    rather than printing the image metric alone and leaving the reader to
    assume it means convergence. That is the whole discipline this function
    exists to enforce.
    """
    lines = [
        f"Convergence of the running average -- {curve['directory']}",
        f"  {curve['metric'].upper()} against the previous frame, and against "
        f"the final one",
        "",
        f"  {'frame':>6} {'time':>12} {'vs previous':>13} {'vs final':>13}",
    ]
    times = curve.get("times") or []
    previous = {m["frame"]: m for m in curve.get("previous", [])}
    final = {m["frame"]: m for m in curve.get("final", [])}
    for index in sorted(set(previous) | set(final)):
        time = times[index] if index < len(times) else float("nan")
        lines.append(
            f"  {index:>6} {time:>12.6g} "
            f"{_format_metric(previous.get(index)):>13} "
            f"{_format_metric(final.get(index)):>13}")

    lines.append("")
    if physical:
        for name in ("density", "velocity", "temperature"):
            key = f"{name}_mean_rel_error"
            if key in physical:
                lines.append(f"  centreline {name:<12} "
                             f"{100.0 * float(physical[key]):.3f}% mean "
                             f"relative error against Cai")
    else:
        lines.append(
            "  NO PHYSICAL METRIC ALONGSIDE. The numbers above say the picture "
            "stopped\n"
            "  changing, which is not the same as the answer being right -- run "
            "./Allpost\n"
            "  so results/metrics.yaml exists and read them together.")
    lines += ["", "  " + curve["caveat"]]
    return lines


def _format_metric(record) -> str:
    """A metric cell: ``saturated`` rather than ``inf``."""
    if record is None:
        return "-"
    if record.get("saturated"):
        return "saturated"
    value = record.get("value")
    return "-" if value is None else f"{float(value):.2f}"


def write_curve(path: Path, curves: Sequence[dict]) -> Path:
    """Write the convergence curves for one or more series."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by plumetools.viz.dissect -- do not edit.\n"
        "# Image-difference curves. Read the 'caveat' before quoting these.\n\n"
        + yaml.safe_dump({"curves": list(curves)}, sort_keys=False,
                         default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


# --------------------------------------------------------------------------- #
# reading a render manifest
# --------------------------------------------------------------------------- #

def series_frames(viz_dir: Path, field: str) -> tuple[list[str], list[float]]:
    """One field's frames and their solver times, from ``manifest.yaml``.

    Args:
        viz_dir: a ``results/viz`` directory.
        field: the ``requested`` field name, e.g. ``dsmcRhoN``.

    Returns:
        ``(filenames, times)``, ordered by frame number and relative to the
        field's own folder.

    Raises:
        DissectError: if there is no manifest.

    Order comes from the manifest's frame numbers, not from the filenames --
    the same reason :mod:`plumetools.viz.video` does: ``%g``-formatted times
    sort wrongly, and a curve computed on out-of-order frames measures
    something that never happened.
    """
    manifest = Path(viz_dir) / "manifest.yaml"
    if not manifest.is_file():
        raise DissectError(
            f"no manifest.yaml in {viz_dir}; it records which frame is which, "
            f"and the filenames are not a reliable order")

    document = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    records = []
    for entry in document.get("images") or []:
        requested = str(entry.get("requested") or entry.get("field") or "")
        if requested != field or not entry.get("file"):
            continue
        records.append((int(entry.get("frame", 0)),
                        Path(entry["file"]),
                        float(entry.get("time", float("nan")))))
    records.sort()
    return ([str(path.name) for _, path, _ in records],
            [time for _, _, time in records])


def series_directory(viz_dir: Path, field: str) -> Path:
    """Where a field's frames live, from the manifest's own relative paths."""
    manifest = Path(viz_dir) / "manifest.yaml"
    document = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    for entry in document.get("images") or []:
        requested = str(entry.get("requested") or entry.get("field") or "")
        if requested == field and entry.get("file"):
            return (Path(viz_dir) / entry["file"]).parent
    raise DissectError(f"{viz_dir} has no frames for {field!r}")


# --------------------------------------------------------------------------- #
# command line
# --------------------------------------------------------------------------- #

def _sheet(args) -> int:
    """``dissect.py sheet --matrix matrix.yaml``."""
    document = yaml.safe_load(Path(args.matrix).read_text(encoding="utf-8")) or {}
    rows = [str(r) for r in document.get("rows") or []]
    columns = [str(c) for c in document.get("columns") or []]
    cells = [
        Cell(row=int(entry["row"]), column=int(entry["column"]),
             path=Path(entry["image"]), label=str(entry.get("label", "")))
        for entry in document.get("cells") or []
    ]
    if not cells:
        print(f"dissect: {args.matrix} lists no cells")
        return 1

    font = args.font or find_font()
    if font is None:
        print("dissect: no drawtext font found; the sheet will be unlabelled "
              "and the legend beside it says which cell is which.")

    output = Path(args.output)
    try:
        build_contact_sheet(cells, output, cell_width=args.cell_width,
                            cell_height=args.cell_height, font=font,
                            trim_bottom=args.trim_bottom, ffmpeg=args.ffmpeg)
    except DissectError as exc:
        print(f"dissect: {exc}")
        return 1
    write_legend(output.with_suffix(".yaml"), cells,
                 title=str(document.get("title", "contact sheet")),
                 rows=rows, columns=columns,
                 notes=document.get("notes") or [])
    print(f"    wrote {output}  ({len(cells)} cell(s), "
          f"{len(rows)} x {len(columns)})")
    return 0


def _converge(args) -> int:
    """``dissect.py converge <viz-dir> ... --field dsmcRhoN``."""
    curves = []
    for directory in args.viz_dir:
        viz_dir = Path(directory)
        try:
            frames, times = series_frames(viz_dir, args.field)
            if len(frames) < 2:
                print(f"### {viz_dir}: {len(frames)} frame(s) of "
                      f"{args.field}; a curve needs two")
                continue
            curve = convergence_curve(series_directory(viz_dir, args.field),
                                      frames, times=times, metric=args.metric,
                                      ffmpeg=args.ffmpeg)
        except DissectError as exc:
            print(f"### {viz_dir}: {exc}")
            continue

        physical = None
        metrics = viz_dir.parent / "metrics.yaml"
        if metrics.is_file():
            physical = (yaml.safe_load(metrics.read_text(encoding="utf-8"))
                        or {}).get("centerline")
        print("\n".join(convergence_report(curve, physical)))
        print()
        curves.append(curve)

    if not curves:
        print("dissect: nothing measured")
        return 1
    written = write_curve(Path(args.output), curves)
    print(f"    wrote {written}  ({len(curves)} series)")
    return 0


def main(argv: Sequence[str]) -> int:
    """``python plumetools/viz/dissect.py {sheet,converge} ...``."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="dissect.py",
        description="Contact sheet and convergence curves from rendered frames.")
    parser.add_argument("--ffmpeg", default="ffmpeg")
    sub = parser.add_subparsers(dest="command", required=True)

    sheet = sub.add_parser("sheet", help="assemble a matrix contact sheet")
    sheet.add_argument("--matrix", required=True,
                       help="a YAML file listing rows, columns and cells")
    sheet.add_argument("--output", default="contact-sheet.png")
    sheet.add_argument("--cell-width", type=int, default=300)
    sheet.add_argument("--cell-height", type=int, default=780)
    sheet.add_argument("--trim-bottom", type=int, default=0,
                       help="pixels of colour bar to crop off each cell; with "
                            "the range pinned across the matrix every bar is "
                            "the same one")
    sheet.add_argument("--font", default=None)
    sheet.set_defaults(handler=_sheet)

    converge = sub.add_parser("converge",
                              help="image-difference curves per series")
    converge.add_argument("viz_dir", nargs="+")
    converge.add_argument("--field", default="dsmcRhoN")
    converge.add_argument("--metric", default="psnr", choices=("psnr", "ssim"))
    converge.add_argument("--output", default="convergence.yaml")
    converge.set_defaults(handler=_converge)

    args = parser.parse_args(list(argv[1:]))
    if not have_ffmpeg(args.ffmpeg):
        print(f"dissect: {args.ffmpeg} is not on PATH; nothing to do.")
        return 0
    return args.handler(args)


if __name__ == "__main__":
    import sys

    raise SystemExit(main(sys.argv))
