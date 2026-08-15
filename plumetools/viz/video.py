r"""Encode each folder of frames into a video, with ffmpeg.

Runs after the rendering, over the PNGs already on disk. **No ParaView**, which
is the point: re-encoding at a different frame rate is seconds of work and does
not mean re-rendering a study for twenty minutes. It is importable, testable and
runnable with plain ``python``::

    python plumetools/viz/video.py cases/cai2012/Cases/Kn100/results/viz
    python plumetools/viz/video.py <viz-dir> --framerate 10
    python plumetools/viz/video.py <viz-dir> --spec cases/cai2012/viz.yaml

One video per folder of frames, written into that same folder beside them:

.. code-block:: text

    results/viz/rhoN/rhoN_0000.png ... rhoN_0004.png
                     rhoN.mp4

Three things ffmpeg will not do for you
---------------------------------------
1. **Odd pixel dimensions are fatal.** ``libx264`` with ``yuv420p`` needs both
   sides even, and a slice's width comes from the geometry's aspect ratio, so it
   is odd about half the time -- the plume images here are 345 x 800:

   .. code-block:: text

       [libx264] width not divisible by 2 (345x800)
       Error while opening encoder

   Every command therefore carries ``pad=ceil(iw/2)*2:ceil(ih/2)*2``, which adds
   at most one pixel on each side.

2. **Frame order is taken from the manifest, not from the filenames.** The
   obvious ``-pattern_type glob`` sorts lexicographically, which is right only
   while the names happen to be zero-padded. ``docs/viz-slices.md`` warns that a
   ``{time}``-named series sorts wrongly -- ``%g`` gives ``0.0021211`` and
   ``0.00990533`` different widths -- and an encoder that inherited that trap
   would silently play the frames out of order. A concat list built from the
   recorded frame numbers cannot.

3. **``-pix_fmt yuv420p`` is not the default** for a PNG input, and without it
   the result plays in ffplay and VLC but not in a browser or QuickTime.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

import yaml

#: Name of the concat list written beside the frames while encoding, then
#: removed. Leading dot so it sorts out of the way if a run is interrupted.
CONCAT_NAME = ".frames.ffconcat"

#: Codecs that take ``-crf``. Anything else gets the encoder's own default,
#: because passing ``-crf`` to a codec without it is a hard error, not a
#: warning.
CRF_CODECS = frozenset({"libx264", "libx264rgb", "libx265", "libvpx", "libvpx-vp9"})

#: ffmpeg's filter for rounding a frame up to even dimensions. See the module
#: docstring: without it, half of all slice images fail to encode at all.
EVEN_DIMENSIONS = "pad=ceil(iw/2)*2:ceil(ih/2)*2"


class VideoError(RuntimeError):
    """ffmpeg failed, or could not be run."""


@dataclass(frozen=True)
class Series:
    """One folder's worth of frames, in the order they should be played.

    Attributes:
        directory: where the frames are, and where the video goes.
        name: the field the frames are of; also the video's basename.
        frames: file names, ordered by frame number.
    """

    directory: Path
    name: str
    frames: tuple[str, ...]

    def output(self, container: str) -> Path:
        """Where this series' video goes."""
        return self.directory / f"{self.name}.{container}"


@dataclass
class EncodedVideo:
    """One video that was written."""

    path: Path
    name: str
    frames: int
    framerate: float

    @property
    def seconds(self) -> float:
        """How long it plays."""
        return self.frames / self.framerate if self.framerate else 0.0


# --------------------------------------------------------------------------- #
# building the command
# --------------------------------------------------------------------------- #

def concat_list(frames: Sequence[str]) -> str:
    """The ffconcat document listing *frames* in playing order.

    No ``duration`` directives: every frame is shown for the same time, so the
    input frame rate sets it, and a trailing ``duration`` would need the last
    file repeated -- which silently adds one extra frame to every video.

    Args:
        frames: file names, relative to the directory ffmpeg will run in.

    Returns:
        The file's contents.
    """
    lines = ["ffconcat version 1.0"]
    for name in frames:
        # Single quotes are the concat demuxer's escape; a quote inside a name
        # has to be closed, escaped and reopened.
        escaped = name.replace("'", r"'\''")
        lines.append(f"file '{escaped}'")
    return "\n".join(lines) + "\n"


def ffmpeg_command(list_file: str, output: str, framerate: float,
                   codec: str = "libx264", quality: int = 18,
                   ffmpeg: str = "ffmpeg") -> list[str]:
    """The argv for encoding one series.

    Pure, so the flags that matter can be checked without ffmpeg installed --
    and three of them do matter; see the module docstring.

    Args:
        list_file: the concat list, relative to the working directory.
        output: the video, relative to the working directory.
        framerate: frames per second, set on both input and output. On the
            input it decides how long each PNG is held; on the output it is the
            stream's rate.
        codec: an ffmpeg video encoder.
        quality: ``-crf``, lower being better. Ignored by codecs without it.
        ffmpeg: the binary.

    Returns:
        The argument vector.
    """
    command = [
        ffmpeg, "-y", "-hide_banner", "-loglevel", "error",
        "-r", f"{framerate:g}",
        "-f", "concat", "-safe", "0",
        "-i", list_file,
        "-vf", EVEN_DIMENSIONS,
        "-c:v", codec,
        "-pix_fmt", "yuv420p",
    ]
    if codec in CRF_CODECS:
        command += ["-crf", str(int(quality))]
    command += ["-r", f"{framerate:g}", output]
    return command


# --------------------------------------------------------------------------- #
# finding the frames
# --------------------------------------------------------------------------- #

def series_from_records(records: Iterable[tuple[Path, str, int]],
                        min_frames: int = 2) -> list[Series]:
    """Group rendered frames into series, one per folder and field.

    Grouping is by *both*, not by folder alone: a flat layout puts every field's
    images in one directory, and grouping by folder there would splice unrelated
    fields into a single video.

    Args:
        records: ``(path, field, frame_number)`` for each image written.
        min_frames: series shorter than this are not videos. One frame is a
            still, and encoding it produces a file nobody wants.

    Returns:
        Series in a stable order, each with its frames sorted by frame number.
    """
    grouped: dict[tuple[Path, str], list[tuple[int, str]]] = {}
    for path, name, frame in records:
        grouped.setdefault((path.parent, name), []).append((frame, path.name))

    series = []
    for (directory, name), frames in sorted(grouped.items(),
                                            key=lambda item: str(item[0])):
        if len(frames) < min_frames:
            continue
        ordered = tuple(filename for _, filename in sorted(frames))
        series.append(Series(directory=directory, name=name, frames=ordered))
    return series


def series_from_manifest(viz_dir: Path, min_frames: int = 2) -> list[Series]:
    """Read ``manifest.yaml`` and group what it records into series.

    The manifest is the authority on frame order -- see the module docstring on
    why the filenames are not.

    Args:
        viz_dir: a ``results/viz`` directory.
        min_frames: as :func:`series_from_records`.

    Returns:
        The series found.

    Raises:
        VideoError: if there is no manifest to read.
    """
    manifest = Path(viz_dir) / "manifest.yaml"
    if not manifest.is_file():
        raise VideoError(
            f"no manifest.yaml in {viz_dir}. It records which frame is which, "
            f"and is written unless output.manifest was turned off.")

    document = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
    records = []
    for entry in document.get("images") or []:
        relative = entry.get("file")
        if not relative:
            continue
        records.append((
            Path(viz_dir) / relative,
            str(entry.get("requested") or entry.get("field") or "video"),
            int(entry.get("frame", 0)),
        ))
    return series_from_records(records, min_frames=min_frames)


# --------------------------------------------------------------------------- #
# encoding
# --------------------------------------------------------------------------- #

def have_ffmpeg(ffmpeg: str = "ffmpeg") -> bool:
    """Whether *ffmpeg* can be found on PATH."""
    return shutil.which(ffmpeg) is not None


def encode_series(series: Series, framerate: float = 4.0,
                  codec: str = "libx264", container: str = "mp4",
                  quality: int = 18, ffmpeg: str = "ffmpeg") -> EncodedVideo:
    """Encode one series into a video beside its frames.

    ffmpeg runs *in* the frame directory, so the concat list holds bare file
    names. That keeps ``-safe 0`` honest and means a path with spaces or an
    absolute prefix cannot confuse the demuxer.

    Args:
        series: the frames, in order.
        framerate: frames per second.
        codec: an ffmpeg video encoder.
        container: the file extension, and so the muxer.
        quality: ``-crf`` for codecs that take one.
        ffmpeg: the binary.

    Returns:
        A record of the video written.

    Raises:
        VideoError: if ffmpeg is missing, fails, or writes nothing.
    """
    if not have_ffmpeg(ffmpeg):
        raise VideoError(f"{ffmpeg} is not on PATH")

    output = series.output(container)
    listing = series.directory / CONCAT_NAME
    listing.write_text(concat_list(series.frames), encoding="utf-8",
                       newline="\n")
    try:
        result = subprocess.run(
            ffmpeg_command(CONCAT_NAME, output.name, framerate, codec,
                           quality, ffmpeg),
            cwd=series.directory, capture_output=True, text=True)
    finally:
        listing.unlink(missing_ok=True)

    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        raise VideoError(
            f"ffmpeg failed for {series.name} (exit {result.returncode}): "
            + (detail[0] if detail else "no output"))
    if not output.is_file() or output.stat().st_size == 0:
        # ffmpeg exits 0 having written nothing when every frame was rejected
        # by the encoder, so the exit code alone does not mean success.
        raise VideoError(f"ffmpeg wrote no video for {series.name}")

    return EncodedVideo(path=output, name=series.name,
                        frames=len(series.frames), framerate=framerate)


def encode_all(series: Sequence[Series], framerate: float = 4.0,
               codec: str = "libx264", container: str = "mp4",
               quality: int = 18, ffmpeg: str = "ffmpeg",
               report=print) -> list[EncodedVideo]:
    """Encode every series, reporting as it goes.

    One series failing does not stop the others: the videos are a convenience
    over frames that are already on disk, and losing the rest of them because
    one field would not encode is a poor trade.

    Args:
        series: what to encode.
        framerate, codec, container, quality, ffmpeg: as
            :func:`encode_series`.
        report: where progress lines go.

    Returns:
        The videos successfully written.
    """
    written: list[EncodedVideo] = []
    for entry in series:
        try:
            video = encode_series(entry, framerate, codec, container, quality,
                                  ffmpeg)
        except VideoError as exc:
            report(f"    video FAILED {entry.name}: {exc}")
            continue
        written.append(video)
        report(f"    video {video.path.name:<28} "
               f"{video.frames} frames @ {framerate:g} fps "
               f"= {video.seconds:.2g}s")
    return written


def main(argv: Sequence[str]) -> int:
    """``python plumetools/viz/video.py <viz-dir> [...]``."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="video.py",
        description="Encode each folder of rendered frames into a video.")
    parser.add_argument("viz_dir", nargs="+",
                        help="a results/viz directory holding manifest.yaml")
    parser.add_argument("--spec", default=None, metavar="FILE",
                        help="take the settings from a sampling spec's video: "
                             "section, so a re-encode matches what the study "
                             "asked for. Explicit flags still win.")
    parser.add_argument("--framerate", type=float, default=None)
    parser.add_argument("--codec", default=None)
    parser.add_argument("--container", default=None)
    parser.add_argument("--quality", type=int, default=None,
                        help="-crf; lower is better (default 18)")
    parser.add_argument("--ffmpeg", default=None)
    args = parser.parse_args(list(argv[1:]))

    from plumetools.viz.spec import VideoSpec, VizSpecError, load_spec

    settings = VideoSpec()
    if args.spec is not None:
        try:
            settings = load_spec(args.spec).video
        except VizSpecError as exc:
            print(f"video: {exc}")
            return 2
    for key in ("framerate", "codec", "container", "quality", "ffmpeg"):
        value = getattr(args, key)
        if value is not None:
            settings = replace(settings, **{key: value})
    args.framerate, args.codec = settings.framerate, settings.codec
    args.container, args.quality = settings.container, settings.quality
    args.ffmpeg = settings.ffmpeg

    if not have_ffmpeg(args.ffmpeg):
        print(f"video: {args.ffmpeg} is not on PATH; nothing to do.")
        return 0

    total = 0
    for directory in args.viz_dir:
        try:
            series = series_from_manifest(Path(directory))
        except VideoError as exc:
            print(f"video: {directory}: {exc}")
            return 1
        if not series:
            print(f"### {directory}: no series of two or more frames")
            continue
        print(f"### {directory}")
        total += len(encode_all(series, args.framerate, args.codec,
                                args.container, args.quality, args.ffmpeg))
    print(f"### {total} video(s)")
    return 0


if __name__ == "__main__":
    import sys

    raise SystemExit(main(sys.argv))
