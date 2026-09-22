"""Tier 0 for :mod:`plumetools.viz.dissect` -- everything except running ffmpeg.

The contact sheet and the convergence curves are built from filter graphs and
concat lists, and every way they go wrong produces output rather than an error:
a grid whose offsets overlap is still an image, a metric parser that matches no
key returns an empty curve that looks converged, and frames compared in
filename order measure a sequence that never happened.

Nothing here needs ffmpeg installed. The commands are checked as argument
vectors, which is the same split :mod:`plumetools.viz.video` draws.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from plumetools.viz.dissect import (
    Cell,
    DissectError,
    build_contact_sheet,
    cell_filter,
    contact_sheet_command,
    convergence_curve,
    convergence_report,
    escape_drawtext,
    find_font,
    frame_pairs,
    grid_layout,
    metric_command,
    parse_stats,
    series_directory,
    series_frames,
    write_curve,
    write_legend,
)

FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"


def matrix(rows=3, columns=3, root="img"):
    return [Cell(row=r, column=c, path=Path(f"{root}/{r}{c}.png"),
                 label=f"r{r}c{c}")
            for r in range(rows) for c in range(columns)]


# --------------------------------------------------------------------------- #
# the contact sheet
# --------------------------------------------------------------------------- #

def test_the_grid_offsets_tile_without_overlapping():
    """Row-major cells at multiples of the cell size. An arithmetic slip here
    stacks two cases on top of each other, and the sheet still renders."""
    cells = matrix(2, 3)
    layout = grid_layout(cells, cell_width=300, cell_height=780)
    assert layout == ("0_0|300_0|600_0|"
                      "0_780|300_780|600_780")


def test_every_cell_of_a_three_by_three_gets_a_distinct_position():
    cells = matrix(3, 3)
    positions = grid_layout(cells, 300, 780).split("|")
    assert len(positions) == 9
    assert len(set(positions)) == 9


def test_the_command_stacks_exactly_the_cells_it_was_given():
    cells = matrix(3, 3)
    command = contact_sheet_command(cells, "sheet.png", font=FONT)
    joined = " ".join(command)
    assert command.count("-i") == 9
    assert "xstack=inputs=9" in joined
    assert "-frames:v 1" in joined
    assert command[-1] == "sheet.png"
    assert "[out]" in command[command.index("-map") + 1]


def test_one_cell_does_not_go_through_xstack():
    """xstack refuses inputs=1, and a one-case sheet is a legitimate thing to
    ask for while the matrix is part way through."""
    command = contact_sheet_command(matrix(1, 1), "sheet.png")
    joined = " ".join(command)
    assert "xstack" not in joined
    assert "null[out]" in joined


def test_a_sheet_of_no_cases_is_refused():
    with pytest.raises(DissectError, match="no cases"):
        contact_sheet_command([], "sheet.png")


def test_a_label_with_a_colon_does_not_truncate_the_filter_graph():
    """`:` separates filter options. An unescaped label -- and a label is
    exactly the sort of string that grows a `t:0.0099` later -- either fails
    the graph or silently drops everything after the colon."""
    assert escape_drawtext("ppc020/s1p5 t:0.0099") == \
        "ppc020/s1p5 t\\:0.0099"
    command = contact_sheet_command(
        [Cell(0, 0, Path("a.png"), "t:1, [x]")], "s.png", font=FONT)
    graph = command[command.index("-filter_complex") + 1]
    assert "t\\:1\\, \\[x\\]" in graph


def test_the_colour_bar_is_only_cropped_when_asked_for():
    """With the range pinned across the matrix, all nine bars are the same
    legend -- but cropping by default would silently remove the scale from a
    sheet whose ranges were NOT pinned."""
    plain = cell_filter(0, Cell(0, 0, Path("a.png")), cell_width=300,
                        cell_height=780, header=40, font=None, font_size=22,
                        trim_bottom=0, background="white")
    trimmed = cell_filter(0, Cell(0, 0, Path("a.png")), cell_width=300,
                          cell_height=780, header=40, font=None, font_size=22,
                          trim_bottom=110, background="white")
    assert "crop=" not in plain
    assert "crop=iw:max(ih-110\\,1):0:0" in trimmed


def test_a_missing_font_gives_an_unlabelled_cell_not_a_broken_graph():
    graph = cell_filter(0, Cell(0, 0, Path("a.png"), "label"), cell_width=300,
                        cell_height=780, header=40, font=None, font_size=22,
                        trim_bottom=0, background="white")
    assert "drawtext" not in graph
    assert graph.startswith("[0:v]") and graph.endswith("[c0]")


def test_a_missing_cell_image_is_refused_before_ffmpeg_runs(tmp_path):
    """A grid quietly missing a case is worse than no grid: the sheet is the
    headline artefact, and nothing in the image says a cell is absent."""
    (tmp_path / "there.png").write_bytes(b"x")
    cells = [Cell(0, 0, tmp_path / "there.png"),
             Cell(0, 1, tmp_path / "missing.png")]
    with pytest.raises(DissectError, match="do not exist"):
        build_contact_sheet(cells, tmp_path / "sheet.png")


def test_the_legend_records_which_cell_is_which(tmp_path):
    """A contact sheet is an image and records nothing about its own
    provenance, and it is the artefact most likely to be pasted somewhere
    alone."""
    cells = matrix(2, 2)
    path = write_legend(tmp_path / "sheet.yaml", cells, title="sweep",
                        rows=["ppc005", "ppc020"], columns=["s0p5", "s1p5"],
                        notes=["range pinned from ppc020/s1p5"])
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert document["rows"] == ["ppc005", "ppc020"]
    assert len(document["cells"]) == 4
    assert document["cells"][0]["image"].endswith("00.png")


def test_find_font_returns_none_rather_than_a_path_that_does_not_exist():
    assert find_font(["/no/such/font.ttf"]) is None


# --------------------------------------------------------------------------- #
# parsing ffmpeg's statistics
# --------------------------------------------------------------------------- #

RGB_STATS = (
    "n:1 mse_avg:1213.13 mse_r:176.85 mse_g:3272.10 mse_b:190.43 "
    "psnr_avg:17.29 psnr_r:25.65 psnr_g:12.98 psnr_b:25.33 \n"
    "n:2 mse_avg:1027.60 mse_r:173.20 mse_g:2733.98 mse_b:175.63 "
    "psnr_avg:18.01 psnr_r:25.75 psnr_g:13.76 psnr_b:25.68 \n"
)

YUV_STATS = (
    "n:1 mse_avg:0.50 mse_y:0.40 mse_u:0.10 mse_v:0.10 "
    "psnr_avg:51.14 psnr_y:52.11 psnr_u:58.13 psnr_v:58.13\n"
)


def test_both_pixel_format_key_sets_parse():
    """ffmpeg emits mse_r/g/b for RGB input and mse_y/u/v for YUV, and which
    one you get depends on the frames. A parser that knew only one would return
    an empty curve for the other -- indistinguishable from a converged one."""
    rgb = parse_stats(RGB_STATS)
    yuv = parse_stats(YUV_STATS)
    assert [record["value"] for record in rgb] == [17.29, 18.01]
    assert [record["value"] for record in yuv] == [51.14]
    assert rgb[0]["mse"] == 1213.13


def test_identical_frames_come_back_as_infinity_not_as_a_number():
    """Two byte-identical frames give psnr inf. Recorded as a number it would
    put a spike at the end of every curve and make an 8-bit palette limit look
    like perfect convergence."""
    records = parse_stats("n:1 mse_avg:0.00 psnr_avg:inf\n")
    assert records[0]["value"] == float("inf")


def test_lines_without_a_frame_number_are_ignored():
    """ffmpeg interleaves warnings on the same stream the statistics come
    down."""
    noisy = ("[null @ 0x1] Application provided invalid, non monotonically "
             "increasing dts to muxer in stream 0: 1 >= 0\n") + RGB_STATS
    assert len(parse_stats(noisy)) == 2


def test_ssim_reads_its_own_summary_key():
    records = parse_stats("n:1 R:0.99 G:0.98 B:0.99 All:0.986 (18.7)\n",
                          metric="ssim")
    assert records[0]["value"] == pytest.approx(0.986)


# --------------------------------------------------------------------------- #
# pairing frames
# --------------------------------------------------------------------------- #

FRAMES = ["f_0000.png", "f_0001.png", "f_0002.png", "f_0003.png"]


def test_successive_frames_are_offset_by_exactly_one():
    left, right = frame_pairs(FRAMES, "previous")
    assert left == FRAMES[:-1]
    assert right == FRAMES[1:]
    assert len(left) == len(right)


def test_the_final_frame_reference_compares_everything_with_the_last():
    left, right = frame_pairs(FRAMES, "final")
    assert right == [FRAMES[-1]] * 3
    assert left == FRAMES[:-1]


def test_a_single_frame_is_not_a_curve():
    with pytest.raises(DissectError, match="at least two frames"):
        frame_pairs(["only.png"], "previous")


def test_an_unknown_reference_is_refused_rather_than_defaulted():
    with pytest.raises(DissectError, match="previous"):
        frame_pairs(FRAMES, "nearest")


def test_the_metric_command_gives_both_inputs_a_frame_rate():
    """The concat demuxer gives stills no timestamps. Without -r, ffmpeg emits
    a dts warning per frame onto the same stream the statistics are parsed
    from."""
    command = metric_command("a.ffconcat", "b.ffconcat")
    assert command.count("-r") == 2
    assert command[command.index("-lavfi") + 1] == "psnr=stats_file=-"
    assert command[-3:] == ["-f", "null", "-"]


def test_an_unknown_metric_is_refused():
    with pytest.raises(DissectError, match="psnr or ssim"):
        metric_command("a", "b", metric="rmse")


# --------------------------------------------------------------------------- #
# the report
# --------------------------------------------------------------------------- #

def _curve():
    return {
        "directory": "case/results/viz/dsmcRhoN",
        "metric": "psnr",
        "n_frames": 3,
        "times": [0.006, 0.008, 0.010],
        "previous": [{"frame": 1, "reference": "previous", "value": 41.2,
                      "mse": 5.0, "saturated": False},
                     {"frame": 2, "reference": "previous", "value": None,
                      "mse": 0.0, "saturated": True}],
        "final": [{"frame": 1, "reference": "final", "value": 40.0,
                   "mse": 6.0, "saturated": False},
                  {"frame": 2, "reference": "final", "value": None,
                   "mse": 0.0, "saturated": True}],
        "caveat": "PSNR between log-scaled, colour-mapped PNGs.",
    }


def test_the_caveat_travels_with_the_data_not_only_the_report(tmp_path):
    """A column of PSNR values with no note attached is how an image metric
    ends up quoted as a physical one."""
    path = write_curve(tmp_path / "convergence.yaml", [_curve()])
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert "caveat" in document["curves"][0]


def test_the_report_refuses_to_stand_alone_without_a_physical_metric():
    text = "\n".join(convergence_report(_curve(), physical=None))
    assert "NO PHYSICAL METRIC ALONGSIDE" in text
    assert _curve()["caveat"] in text


def test_the_physical_error_is_printed_beside_the_image_metric():
    text = "\n".join(convergence_report(
        _curve(), physical={"density_mean_rel_error": 0.0107}))
    assert "1.070%" in text
    assert "NO PHYSICAL METRIC" not in text


def test_saturation_is_labelled_rather_than_printed_as_a_number():
    """`inf` in a column of dB reads as the best possible result. It means the
    change fell below one colour step."""
    text = "\n".join(convergence_report(_curve()))
    assert "saturated" in text
    assert "inf" not in text


# --------------------------------------------------------------------------- #
# reading a render manifest
# --------------------------------------------------------------------------- #

def write_manifest(tmp_path, images):
    viz = tmp_path / "viz"
    viz.mkdir(parents=True, exist_ok=True)
    (viz / "manifest.yaml").write_text(
        yaml.safe_dump({"images": images}), encoding="utf-8")
    return viz


def test_frames_come_back_in_frame_number_order_not_manifest_order(tmp_path):
    """Order comes from the manifest's frame numbers, not from the filenames
    and not from the order the entries happen to appear in. %g-formatted times
    sort wrongly in every glob, and a curve computed on out-of-order frames
    measures a sequence that never happened."""
    viz = write_manifest(tmp_path, [
        {"file": "dsmcRhoN/dsmcRhoN_0002.png", "requested": "dsmcRhoN",
         "frame": 2, "time": 0.00990533},
        {"file": "dsmcRhoN/dsmcRhoN_0000.png", "requested": "dsmcRhoN",
         "frame": 0, "time": 0.0021211},
        {"file": "dsmcRhoN/dsmcRhoN_0001.png", "requested": "dsmcRhoN",
         "frame": 1, "time": 0.00424357},
    ])
    frames, times = series_frames(viz, "dsmcRhoN")
    assert frames == ["dsmcRhoN_0000.png", "dsmcRhoN_0001.png",
                      "dsmcRhoN_0002.png"]
    assert times == [0.0021211, 0.00424357, 0.00990533]


def test_only_the_requested_field_is_taken(tmp_path):
    viz = write_manifest(tmp_path, [
        {"file": "rhoN/rhoN_0000.png", "requested": "rhoN", "frame": 0,
         "time": 0.002},
        {"file": "dsmcRhoN/dsmcRhoN_0000.png", "requested": "dsmcRhoN",
         "frame": 0, "time": 0.002},
    ])
    assert series_frames(viz, "dsmcRhoN")[0] == ["dsmcRhoN_0000.png"]
    assert series_directory(viz, "dsmcRhoN").name == "dsmcRhoN"


def test_a_missing_manifest_says_why_the_filenames_are_not_enough(tmp_path):
    (tmp_path / "viz").mkdir()
    with pytest.raises(DissectError, match="not a reliable order"):
        series_frames(tmp_path / "viz", "dsmcRhoN")
