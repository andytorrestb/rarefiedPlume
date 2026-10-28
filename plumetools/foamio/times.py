"""Numeric time directories: the one rule every reader in plumetools applies.

A time directory is a subdirectory whose name ``float()`` accepts. Time ``0``
holds the initial fields, which no solver wrote, so every query here starts
strictly after a given time -- ``0`` by default. Names are returned exactly as
OpenFOAM wrote them: ``0.00990534`` and ``0.0099053`` are different
directories, and comparing them as floats would invent a tolerance where there
is an exact answer.

Twelve copies of this loop had grown up across the families; they differed only
in what they returned, which is what the three functions below cover.
"""

from __future__ import annotations

from pathlib import Path

from plumetools.foamio.errors import FoamDataError


def time_value(name: str) -> float | None:
    """The time a directory name stands for, or ``None`` if it is not a time."""
    try:
        return float(name)
    except ValueError:
        return None


def time_directories(directory: Path, *, after: float = 0.0,
                     missing_ok: bool = False) -> list[tuple[float, str]]:
    """``(time, name)`` of every time directory later than ``after``, ascending.

    Args:
        directory: a case directory, or one ``processor<N>`` directory of it.
        after: only times strictly greater than this [s].
        missing_ok: return ``[]`` for a directory that does not exist, rather
            than raising ``FileNotFoundError``.

    Returns:
        Sorted by time, then by name.
    """
    directory = Path(directory)
    if missing_ok and not directory.is_dir():
        return []
    found = []
    for entry in directory.iterdir():
        if not entry.is_dir():
            continue
        value = time_value(entry.name)
        if value is not None and value > after:
            found.append((value, entry.name))
    return sorted(found)


def time_names(directory: Path, *, after: float = 0.0,
               missing_ok: bool = False) -> list[str]:
    """The names of :func:`time_directories`, ascending."""
    return [name for _, name in time_directories(directory, after=after,
                                                 missing_ok=missing_ok)]


def latest_time(case_dir: Path) -> str:
    """The name of the largest time directory past ``0``.

    Raises:
        FoamDataError: if there is none -- the solver wrote nothing.
    """
    found = time_directories(case_dir)
    if not found:
        raise FoamDataError(
            f"{case_dir} has no time directory past 0; the solver wrote nothing. "
            f"Run ./Allrun first.")
    return found[-1][1]
