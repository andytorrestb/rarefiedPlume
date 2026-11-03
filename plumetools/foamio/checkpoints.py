r"""Restart checkpoints of a case, decomposed or not, and what a resume must do.

A checkpoint is a time directory a solver can continue from. In a decomposed
case it exists once per rank, and it is usable only when **every** rank wrote
**all** of it. ``startFrom latestTime`` makes each rank take its own newest
directory, so a partial one is never harmless: a rank that lacks it aborts the
run, and a rank that has it without its parcel cloud resumes that subdomain
with no parcels at all, which nothing reports.

What "complete" means here
--------------------------
A time ``T`` is complete on a rank when all of these hold:

* ``T/uniform/time`` exists -- ``Time`` writes it on every rank at every write;
* every field file in ``T/`` and in ``T/lagrangian/<cloud>/`` ends with
  OpenFOAM's end divider (``// ***...*** //``), which is written last, so a
  file cut off by a killed job does not;
* ``T/`` holds every field file any other rank holds at ``T``, and its cloud
  holds every file any other rank's cloud holds;
* a rank whose previous time carried a cloud carries one at ``T`` too. (A rank
  whose subdomain has never held a parcel legitimately writes none.)

A time is a complete **decomposed** checkpoint when it is complete on every one
of ``processor0 .. processor<N-1>``. Times in the case directory itself -- a
serial run's, or ``reconstructPar``'s -- are checked the same way, as one rank.
Only the newest times are inspected file by file: the walk stops at the first
complete one, since nothing older is resumed from.

Nothing here deletes, moves or edits anything.

.. code-block:: text

    python -m plumetools.foamio.checkpoints report   <case>
    python -m plumetools.foamio.checkpoints finished <case>
    python -m plumetools.foamio.checkpoints plan     <case> --np N

``plan`` is what a case's ``Allrun`` reads: shell assignments on stdout, and
exit status 3 with the reason on stderr when the case must not be continued
as it is.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from plumetools.foamio.times import time_directories

#: What OpenFOAM writes last in every field file.
END_DIVIDER = b"// " + b"*" * 73 + b" //"

#: Written by ``Time`` on every rank at every write.
TIME_STATE = Path("uniform") / "time"

_PROCESSOR = re.compile(r"processor(\d+)$")


# --------------------------------------------------------------------------- #
# what one rank wrote at one time
# --------------------------------------------------------------------------- #

def _ends_with_divider(path: Path) -> bool:
    if path.suffix == ".gz":                 # writeCompression on: no plain footer
        return path.stat().st_size > 0
    with open(path, "rb") as handle:
        handle.seek(0, 2)
        size = handle.tell()
        handle.seek(max(0, size - 4 * len(END_DIVIDER)))
        return END_DIVIDER in handle.read()


def _files(directory: Path) -> set:
    if not directory.is_dir():
        return set()
    return {entry.name for entry in directory.iterdir() if entry.is_file()}


def _clouds(time_dir: Path) -> dict:
    lagrangian = time_dir / "lagrangian"
    if not lagrangian.is_dir():
        return {}
    return {cloud.name: _files(cloud) for cloud in lagrangian.iterdir()
            if cloud.is_dir()}


@dataclass(frozen=True)
class _Write:
    fields: set
    clouds: dict
    problems: tuple


def _inspect(time_dir: Path) -> _Write:
    problems = []
    if not (time_dir / TIME_STATE).is_file():
        problems.append(f"no {TIME_STATE.as_posix()}")
    fields = _files(time_dir)
    clouds = _clouds(time_dir)
    paths = [time_dir / name for name in fields]
    for cloud, names in clouds.items():
        paths += [time_dir / "lagrangian" / cloud / name for name in names]
    truncated = sorted(p.relative_to(time_dir).as_posix() for p in paths
                       if not _ends_with_divider(p))
    if truncated:
        problems.append("cut off before the end: " + ", ".join(truncated))
    return _Write(fields=fields, clouds=clouds, problems=tuple(problems))


@dataclass(frozen=True)
class _Rank:
    """One directory that holds time directories: ``processor<N>``, or the case."""

    label: str
    base: Path
    times: tuple          # ((value, name), ...) ascending

    def has(self, name: str) -> bool:
        return any(n == name for _, n in self.times)

    def previous(self, name: str) -> str | None:
        names = [n for _, n in self.times]
        index = names.index(name)
        return names[index - 1] if index else None


def _newest_complete(ranks: list) -> tuple:
    """Walk the union of the ranks' times, newest first, to the first complete one.

    Returns:
        ``(complete, incomplete)``: ``(time, name)`` or ``None``; and
        ``[(time, name, problems), ...]`` for every newer time.
    """
    union = sorted({t for rank in ranks for t in rank.times}, reverse=True)
    incomplete = []
    for value, name in union:
        problems = [f"{rank.label} has no {name}/" for rank in ranks
                    if not rank.has(name)]
        writes = {rank.label: (rank, _inspect(rank.base / name))
                  for rank in ranks if rank.has(name)}
        every_field = set().union(*(w.fields for _, w in writes.values()))
        every_cloud = {}
        for _, write in writes.values():
            for cloud, names in write.clouds.items():
                every_cloud.setdefault(cloud, set()).update(names)
        for label, (rank, write) in writes.items():
            problems += [f"{label}: {p}" for p in write.problems]
            lacking = sorted(every_field - write.fields)
            if lacking:
                problems.append(f"{label}: lacks {', '.join(lacking)}")
            before = rank.previous(name)
            had = _clouds(rank.base / before) if before else {}
            for cloud in sorted(set(every_cloud) | set(had)):
                if cloud in write.clouds:
                    short = sorted(every_cloud[cloud] - write.clouds[cloud])
                    if short:
                        problems.append(f"{label}: cloud '{cloud}' lacks "
                                        f"{', '.join(short)}")
                elif cloud in had:
                    problems.append(f"{label}: wrote a '{cloud}' cloud at "
                                    f"{before} but none at {name}")
        if not problems:
            return (value, name), incomplete
        incomplete.append((value, name, tuple(problems)))
    return None, incomplete


# --------------------------------------------------------------------------- #
# the whole case
# --------------------------------------------------------------------------- #

def rank_directories(case_dir: Path) -> list[Path]:
    """``processor<N>`` directories, by rank number."""
    ranks = []
    for entry in Path(case_dir).iterdir():
        match = _PROCESSOR.match(entry.name)
        if match and entry.is_dir():
            ranks.append((int(match.group(1)), entry))
    return [path for _, path in sorted(ranks)]


@dataclass(frozen=True)
class Checkpoints:
    """The restartable state of one case.

    Attributes:
        case_dir: the case.
        n_ranks: the ``processor<N>`` directories found.
        missing_ranks: rank numbers below the highest one with no directory.
        reconstructed: ``(time, name)`` of the newest complete time in the case
            directory itself, or ``None``.
        decomposed: ``(time, name)`` of the newest time complete on every rank,
            or ``None``.
        incomplete: ``{"case" | "processors": [(time, name, problems), ...]}``
            for every time newer than the newest complete one there.
    """

    case_dir: Path
    n_ranks: int
    missing_ranks: tuple
    reconstructed: tuple | None
    decomposed: tuple | None
    incomplete: dict = field(default_factory=dict)

    @property
    def latest(self) -> tuple | None:
        """``(time, name, where)`` of the newest usable checkpoint, or ``None``.

        ``where`` is ``"case"`` or ``"processors"``. On a tie the processors
        win: they are what a parallel resume reads without decomposing again.
        """
        candidates = []
        if self.decomposed:
            candidates.append((*self.decomposed, "processors"))
        if self.reconstructed:
            candidates.append((*self.reconstructed, "case"))
        return max(candidates, key=lambda c: (c[0], c[2] == "processors"),
                   default=None)

    def blocking(self, newer_than: float | None,
                 where=("processors", "case")) -> list[str]:
        """Partial times newer than ``newer_than`` in ``where``, as messages.

        ``startFrom latestTime`` -- and ``decomposePar -latestTime``, and
        ``reconstructPar -newTimes`` -- would pick these up, each rank its own
        newest, so a run must not be continued while they are in place.
        """
        newest = float("-inf") if newer_than is None else newer_than
        lines = []
        for place in where:
            for value, name, problems in self.incomplete.get(place, ()):
                if value > newest:
                    lines += [f"{place} t = {name}: {p}" for p in problems]
        if "processors" in where and self.missing_ranks:
            lines.append(f"processor directories {list(self.missing_ranks)} are "
                         f"missing among 0..{self.n_ranks + len(self.missing_ranks) - 1}")
        return lines


def scan(case_dir: Path) -> Checkpoints:
    """Find the newest complete checkpoint of a case, reconstructed and decomposed."""
    case_dir = Path(case_dir)
    case_rank = _Rank("case", case_dir, tuple(time_directories(case_dir)))
    reconstructed, case_incomplete = _newest_complete([case_rank])

    directories = rank_directories(case_dir)
    numbers = [int(_PROCESSOR.match(d.name).group(1)) for d in directories]
    missing = tuple(sorted(set(range(max(numbers) + 1)) - set(numbers))) if numbers else ()
    ranks = [_Rank(d.name, d, tuple(time_directories(d))) for d in directories]
    decomposed, rank_incomplete = (_newest_complete(ranks) if ranks and not missing
                                   else (None, []))
    incomplete = {}
    if case_incomplete:
        incomplete["case"] = case_incomplete
    if rank_incomplete:
        incomplete["processors"] = rank_incomplete
    return Checkpoints(case_dir=case_dir, n_ranks=len(directories),
                       missing_ranks=missing, reconstructed=reconstructed,
                       decomposed=decomposed, incomplete=incomplete)


# --------------------------------------------------------------------------- #
# controlDict, and the questions a runner asks
# --------------------------------------------------------------------------- #

def control_entries(case_dir: Path) -> dict:
    """``endTime`` and ``deltaT`` from ``system/controlDict``, as floats."""
    control = Path(case_dir) / "system" / "controlDict"
    entries = {}
    if control.is_file():
        for line in control.read_text(encoding="utf-8").splitlines():
            match = re.match(r"\s*(endTime|deltaT)\s+([^;]+);", line)
            if match:
                entries[match.group(1)] = float(match.group(2))
    return entries


def at_end(value: float | None, end_time: float | None, delta_t: float) -> bool:
    """``Time::run()``'s own test, with its half-step slack.

    ``endTime`` and a time directory's name are the same number printed by two
    formatters -- 0.00990534 against 0.00990533 -- so an exact comparison would
    call a finished case unfinished.
    """
    return value is not None and end_time is not None and value >= end_time - 0.5 * delta_t


def reached_end(case_dir: Path) -> bool:
    """Whether the newest COMPLETE checkpoint, on any side, is at ``endTime``.

    A time some rank wrote only partly does not count: a job killed while
    writing its last time has not finished.
    """
    entries = control_entries(case_dir)
    latest = scan(case_dir).latest
    return at_end(latest[0] if latest else None, entries.get("endTime"),
                  entries.get("deltaT", 0.0))


@dataclass(frozen=True)
class RestartPlan:
    """What a parallel ``Allrun`` must do with the state it finds.

    Attributes:
        resume: the time to continue from, or ``None`` for a fresh start.
        case_time: newest complete time in the case directory, or ``None``.
        proc_time: newest time complete on every rank, or ``None``.
        proc_count: ``processor*`` directories present.
        reached_end: the resume time is at ``endTime``.
        decompose: ``"none"`` (the processors hold the resume time and the
            count matches), ``"all"`` (a fresh start), or ``"latest"``
            (decompose only the resume time).
        reconstruct_first: before re-decomposing, reconstruct the processors'
            newer times, which ``decomposePar -force`` would otherwise delete.
        refusal: why the case must not be continued as it is, or ``None``.
    """

    resume: str | None
    case_time: str | None
    proc_time: str | None
    proc_count: int
    reached_end: bool
    decompose: str
    reconstruct_first: bool
    refusal: str | None

    def shell(self) -> str:
        """``KEY=value`` lines a runner script can ``eval``."""
        values = {
            "RESUME": self.resume or "",
            "CASE_TIME": self.case_time or "",
            "PROC_TIME": self.proc_time or "",
            "PROC_COUNT": str(self.proc_count),
            "REACHED_END": "1" if self.reached_end else "0",
            "DECOMPOSE": self.decompose,
            "RECONSTRUCT_FIRST": "1" if self.reconstruct_first else "0",
        }
        for value in values.values():
            if not re.fullmatch(r"[A-Za-z0-9.+\-]*", value):
                raise ValueError(f"refusing to print {value!r} for a shell")
        return "".join(f"{key}={value}\n" for key, value in values.items())


def plan_restart(case_dir: Path, n_subdomains: int) -> RestartPlan:
    """What a run of ``n_subdomains`` must do with this case's state.

    The decisions ``Allrun`` used to make in shell from ``processor0`` alone,
    made from every rank:

    * resume from the newest checkpoint complete everywhere it is needed --
      never from a partial one, which is refused instead (the run would abort,
      or lose a subdomain's parcels);
    * reuse the processor directories when they hold that checkpoint and their
      count is the run's;
    * otherwise decompose only the resume time -- after reconstructing the
      processors' newer times first, if the decomposition is changing and they
      hold state the case directory does not.
    """
    state = scan(case_dir)
    entries = control_entries(case_dir)
    proc_time = state.decomposed[1] if state.decomposed else None
    case_time = state.reconstructed[1] if state.reconstructed else None
    reconstruct_first = False

    if n_subdomains <= 1:
        # A serial dsmcFoam reads the case directory and nothing else.
        latest = (*state.reconstructed, "case") if state.reconstructed else None
        blocking = state.blocking(latest[0] if latest else None, where=("case",))
        decompose = "none"
        if state.decomposed and (latest is None or state.decomposed[0] > latest[0]):
            blocking.append(
                f"the processor directories hold t = {proc_time}, newer than the "
                f"case directory's {case_time or 'nothing'}; a serial run would "
                f"resume from the case directory and ignore them. Reconstruct "
                f"first (reconstructPar -newTimes), or run with -np "
                f"{state.n_ranks}")
    else:
        latest = state.latest
        blocking = state.blocking(latest[0] if latest else None)
        if latest is None:
            decompose = "all"
        elif latest[2] == "processors" and state.n_ranks == n_subdomains:
            decompose = "none"
        else:
            decompose = "latest"
            # The processors hold a newer complete state than the case
            # directory, and decomposePar -force deletes them: reconstruct it.
            # (Not when the case directory already holds that time too.)
            reconstruct_first = latest[2] == "processors" and (
                state.reconstructed is None
                or state.decomposed[0] > state.reconstructed[0])

    refusal = None
    if blocking:
        refusal = ("this case cannot be continued as it is:\n"
                   + "\n".join(f"    {line}" for line in blocking)
                   + "\nAn incomplete time is one a killed job was writing. Move "
                     "it aside on every rank -- 'mv processorN/<t> "
                     "processorN/<t>.incomplete', and the same in the case "
                     "directory -- then run this again. Nothing was changed.")
    return RestartPlan(
        resume=latest[1] if latest else None, case_time=case_time,
        proc_time=proc_time, proc_count=state.n_ranks,
        reached_end=at_end(latest[0] if latest else None, entries.get("endTime"),
                           entries.get("deltaT", 0.0)),
        decompose=decompose, reconstruct_first=reconstruct_first, refusal=refusal)


# --------------------------------------------------------------------------- #
# the command line
# --------------------------------------------------------------------------- #

def report_lines(state: Checkpoints) -> list[str]:
    def show(entry):
        return entry[1] if entry else "none"

    lines = [f"case directory       newest complete t = {show(state.reconstructed)}",
             f"processor0..{max(state.n_ranks - 1, 0):<8} newest complete t = "
             f"{show(state.decomposed)} ({state.n_ranks} rank(s))"]
    for where, entries in state.incomplete.items():
        for _, name, problems in entries:
            lines.append(f"INCOMPLETE {where} t = {name}")
            lines += [f"    {p}" for p in problems]
    return lines


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.foamio.checkpoints",
        description="Restart checkpoints of a case, checked on every rank.")
    steps = parser.add_subparsers(dest="step", required=True)
    for name, text in (("report", "list the newest complete checkpoints"),
                       ("finished", "exit 0 if the newest complete checkpoint "
                                    "is at endTime, else 1"),
                       ("plan", "what Allrun must do; shell assignments")):
        step = steps.add_parser(name, help=text)
        step.add_argument("case_dir", type=Path)
        if name == "plan":
            step.add_argument("--np", type=int, required=True, dest="n_subdomains")
    args = parser.parse_args(argv)
    if not args.case_dir.is_dir():
        print(f"checkpoints: no case directory {args.case_dir}", file=sys.stderr)
        return 2
    if args.step == "report":
        print("\n".join(report_lines(scan(args.case_dir))))
        return 0
    if args.step == "finished":
        return 0 if reached_end(args.case_dir) else 1
    plan = plan_restart(args.case_dir, args.n_subdomains)
    if plan.refusal:
        print(plan.refusal, file=sys.stderr)
        return 3
    print(plan.shell(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
