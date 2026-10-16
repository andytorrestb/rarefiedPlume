"""Start the next mesh iteration of a case: its inputs, and nothing from its run.

The workflow this serves::

    run      ./Allrun                         samples accumulate in <time>/
    inspect  ./Allpost, ParaView              plumetools.rarefaction's fields
    remesh   python -m plumetools.iteration   a NEW directory with the inputs
             edit its case.yaml               the refinement, by hand
             ./Allmesh                        dictionaries and mesh from it
    fresh    ./Allrun                         new parcels from t = 0

The previous case is only ever read. Its mesh, results and configuration stay
exactly where they are, so every iteration can be compared with the one before.

What is copied
--------------
A case **is** its ``case.yaml``, its other YAML (``viz.yaml``), its scripts
(``All*``) and its ``*.foam`` marker -- :data:`CASE_INPUTS`, the files at the top
of the directory. Everything else there was made from them by ``./Allmesh``,
``./Allrun`` or ``./Allpost``, and is made again from the new ``case.yaml``:
the dictionaries, the mesh, ``0/``, the time step and particle weight derived
from the smallest cell.

What a solver could resume from -- a time directory, ``processor*``, the parcels
in ``lagrangian/``, ``fieldAverage``'s accumulators in ``uniform/`` -- is never
copied, by pattern or by ``--include``, and the new directory is checked for all
of it before this returns. ``controlDict``'s ``startFrom latestTime`` therefore
finds nothing but ``0`` to start from.

A case-specific asset the patterns miss (an STL under ``constant/triSurface``)
is named with ``--include``.

Provenance
----------
``iteration.yaml`` in the new directory records the parent, the hash of the
``case.yaml`` it was copied from, its latest time, the repository commit and a
free-text note; ``provenance/`` keeps the parent's ``case.yaml`` and its
diagnostics summary. ``diff provenance/parent-case.yaml case.yaml`` is then
exactly the edit that made this iteration.
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import yaml

#: The files that ARE a case, at the top of its directory.
CASE_INPUTS = ("case.yaml", "*.yaml", "All*", "*.foam")

#: Matched by :data:`CASE_INPUTS` but written by a tool, so never carried over.
GENERATED_FILES = ("iteration.yaml", "case-summary.json")

#: Path components that mark solver state, a mesh, or results. Nothing under
#: one is copied, whatever asked for it.
STATE_PARTS = ("processor*", "polyMesh", "lagrangian", "uniform", "postProcessing",
               "results", "provenance", "0", "0.orig")

PROVENANCE_DIR = "provenance"


class IterationError(ValueError):
    """The iteration cannot be prepared as asked."""


@dataclass(frozen=True)
class Iteration:
    """What :func:`prepare` made."""

    source: Path
    destination: Path
    copied: tuple
    record: dict


def _is_time(name: str) -> bool:
    try:
        float(name)
    except ValueError:
        return False
    return True


def _sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def times(case_dir: Path) -> list[str]:
    """Numeric time directories past 0, in a case and its ``processor0``."""
    found = set()
    for base in (Path(case_dir), Path(case_dir) / "processor0"):
        if base.is_dir():
            found |= {entry.name for entry in base.iterdir()
                      if entry.is_dir() and _is_time(entry.name)
                      and float(entry.name) > 0.0}
    return sorted(found, key=float)


def resumable_state(case_dir: Path) -> list[str]:
    """Everything in a case a solver or mesher would pick up rather than make.

    Empty for a case that ``./Allmesh`` and ``./Allrun`` would build from
    nothing -- which is what a fresh iteration must be.
    """
    case_dir = Path(case_dir)
    problems = []
    for entry in sorted(case_dir.iterdir()) if case_dir.is_dir() else []:
        if entry.is_dir() and (_is_time(entry.name)
                               or fnmatch.fnmatch(entry.name, "processor*")):
            problems.append(f"{entry.name}/ (solver state)")
    for relative in ("constant/polyMesh", "constant/dsmcProperties",
                     "system/controlDict", "postProcessing", "results"):
        if (case_dir / relative).exists():
            problems.append(f"{relative} (generated)")
    return problems


def _forbidden(relative: Path) -> str | None:
    for part in relative.parts:
        if _is_time(part) or any(fnmatch.fnmatch(part, p) for p in STATE_PARTS):
            return part
    return None


def _git_commit(path: Path) -> str | None:
    try:
        result = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                                capture_output=True, text=True, timeout=10,
                                check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def prepare(source: Path, destination: Path, *, note: str | None = None,
            include=()) -> Iteration:
    """Make ``destination`` a fresh case with ``source``'s inputs.

    Args:
        source: a case directory with a ``case.yaml``. Only read.
        destination: a directory that does not exist yet, or is empty.
        note: why this iteration exists -- recorded in ``iteration.yaml``.
        include: further paths, relative to ``source``, to copy as they are.

    Raises:
        IterationError: for a source that is not a case, a destination that is
            not empty or is nested with the source, or an ``include`` that
            names solver state, a mesh or results.
    """
    source = Path(source).resolve()
    destination = Path(destination).resolve()
    if not (source / "case.yaml").is_file():
        raise IterationError(f"{source} has no case.yaml; it is not a case")
    if destination == source or source in destination.parents:
        raise IterationError(
            f"{destination} is inside {source}. A case's scripts treat their "
            f"directory as theirs -- ./Allclean removes by pattern -- so the next "
            f"iteration goes beside it, not in it.")
    if destination in source.parents:
        raise IterationError(f"{destination} contains {source}")
    if destination.exists() and (not destination.is_dir()
                                 or any(destination.iterdir())):
        raise IterationError(
            f"{destination} already exists and is not empty. An iteration is "
            f"prepared once, into a new directory, so it can never inherit "
            f"anything a previous run left there.")

    files = []
    for pattern in CASE_INPUTS:
        for path in sorted(source.glob(pattern)):
            if (path.is_file() and path.name not in GENERATED_FILES
                    and not path.name.startswith("log.") and path not in files):
                files.append(path)
    extras = []
    for item in include:
        relative = Path(item)
        if relative.is_absolute() or ".." in relative.parts:
            raise IterationError(f"--include {item}: give a path relative to "
                                 f"{source}, inside it")
        bad = _forbidden(relative)
        if bad:
            raise IterationError(
                f"--include {item}: '{bad}' is solver state, a mesh or results, "
                f"which a fresh iteration must generate rather than inherit")
        if not (source / relative).exists():
            raise IterationError(f"--include {item}: no such path in {source}")
        extras.append(relative)

    destination.mkdir(parents=True, exist_ok=True)
    copied = []
    for path in files:
        shutil.copy2(path, destination / path.name)      # keeps the exec bit
        copied.append(path.name)
    for relative in extras:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if (source / relative).is_dir():
            shutil.copytree(source / relative, target)
        else:
            shutil.copy2(source / relative, target)
        copied.append(relative.as_posix())

    provenance = destination / PROVENANCE_DIR
    provenance.mkdir()
    shutil.copy2(source / "case.yaml", provenance / "parent-case.yaml")
    parent_summary = source / "results" / "diagnostics_summary.yaml"
    if parent_summary.is_file():
        shutil.copy2(parent_summary, provenance / "parent-diagnostics_summary.yaml")

    parent_record = {}
    if (source / "iteration.yaml").is_file():
        parent_record = yaml.safe_load(
            (source / "iteration.yaml").read_text(encoding="utf-8")) or {}
    parent_times = times(source)
    record = {
        "name": destination.name,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": note,
        "parent": {
            "case_dir": str(source),
            "name": source.name,
            "case_yaml_sha256": _sha256(source / "case.yaml"),
            "latest_time": parent_times[-1] if parent_times else None,
            "diagnostics_summary_sha256": _sha256(parent_summary),
        },
        "lineage": [*parent_record.get("lineage", []), source.name],
        "repository_commit": _git_commit(Path(__file__).resolve().parent),
        "copied": copied,
        "not_copied": "the mesh, every generated dictionary, 0/, time and "
                      "processor directories, logs and results -- ./Allmesh and "
                      "./Allrun make them again from case.yaml",
        "start": "fresh: new parcels from t = 0, no mapped fields, no averaging "
                 "carried over",
    }
    (destination / "iteration.yaml").write_text(
        "# GENERATED by plumetools.iteration -- the record of where this case\n"
        "# came from. diff provenance/parent-case.yaml case.yaml shows what was\n"
        "# changed to make it.\n\n"
        + yaml.safe_dump(record, sort_keys=False), encoding="utf-8", newline="\n")

    left = resumable_state(destination)
    if left:                                  # cannot happen; checked anyway
        raise IterationError(f"{destination} holds {left}; refusing to call it fresh")
    return Iteration(source=source, destination=destination, copied=tuple(copied),
                     record=record)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.iteration",
        description="Prepare the next mesh iteration of a case in a new "
                    "directory: its inputs only. The source is not modified.")
    parser.add_argument("source", help="the case to start from")
    parser.add_argument("destination", help="a new, empty directory")
    parser.add_argument("--note", default=None,
                        help="why this iteration exists; kept in iteration.yaml")
    parser.add_argument("--include", action="append", default=[], metavar="PATH",
                        help="a further input to copy, relative to the source "
                             "(repeatable)")
    args = parser.parse_args(argv)
    try:
        made = prepare(args.source, args.destination, note=args.note,
                       include=args.include)
    except IterationError as exc:
        print(f"iteration: {exc}", file=sys.stderr)
        return 1

    destination = made.destination
    print(f"Prepared {destination}")
    print(f"  from     {made.source}"
          + (f"   (latest time {made.record['parent']['latest_time']})"
             if made.record["parent"]["latest_time"] else ""))
    print(f"  copied   {', '.join(made.copied)}")
    print("  not      the mesh, generated dictionaries, 0/, time and processor")
    print("           directories, logs, results. The source is unchanged.")
    print()
    print("Next:")
    print(f"  1. edit {destination / 'case.yaml'}  (mesh.refinement_regions, ...)")
    print(f"     diff {PROVENANCE_DIR}/parent-case.yaml case.yaml   shows the edit")
    print(f"  2. cd {destination} && ./Allmesh --describe && ./Allmesh")
    print("  3. ./Allrun     -- starts from t = 0: new parcels, new averaging")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
