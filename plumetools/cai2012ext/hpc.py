r"""What one PBS job does around ``./Allrun``, for one case of the extension.

Each case is its own job: ``cd Cases/<case> && qsub runCase.sh``. The script
calls ``./Allrun`` itself, and around it the four steps here, through

.. code-block:: text

    python -m plumetools.cai2012ext.hpc guard    <case_dir>
    python -m plumetools.cai2012ext.hpc label    <case_dir>
    python -m plumetools.cai2012ext.hpc finished <case_dir>
    python -m plumetools.cai2012ext.hpc record   <case_dir> <seconds>

``guard``
    enforce the matrix profile recorded in case metadata. Legacy cases refuse
    Cai's own Kn values; the explicit full_decades profile permits all seven.
    Nothing on the
    run path reads the manifest any more, so this is where the exclusion is
    enforced for a run -- on the case's own ``case.yaml``.
``label``
    re-assess the case from its own ``case.yaml`` (:func:`.generation.derive`
    and :func:`.feasibility.assess`, not the manifest, which may be stale) and,
    if it is unsuitable for production, write a ``NON-PRODUCTION`` marker into
    it listing the unresolved criteria. The case still runs; ``Allcompare.py``
    and the study table carry the label.
``finished``
    exit 0 if the case already reached ``endTime`` (``processor*/`` included,
    with ``Time::run()``'s half-step slack), 1 if not.
``record``
    add the job's wall clock to the family's ``results/cost-model.yaml``
    through :func:`.feasibility.record_measurement`, under an exclusive
    ``fcntl.flock`` on a lock file beside it: jobs for different cases can
    finish at the same moment, and an unlocked read-modify-write loses one.

``finished`` and the measurement are ``cases/cai2012-health/AllrunCases``'
completed-run test and wall-clock bookkeeping, per case.

An **interface-mode** case (an enabled ``interface.yaml``) is handed to
:mod:`.interface_case` at each step: ``guard`` also refuses an infeasible
interface and results built on another mesh, and ``label`` exits 3 -- which
stops ``runCase.sh`` before ``./Allrun`` -- for a production case that is not
a candidate unless its ``run_mode`` is ``exploratory``. A direct-exit case is
handled exactly as before.
"""

from __future__ import annotations

import argparse
import contextlib
import re
import sys
from pathlib import Path

from plumetools.cai2012.config import load_case_config
from plumetools.cai2012ext import feasibility, generation, interface, policy
from plumetools.cai2012ext import interface_case

MARKER = "NON-PRODUCTION"


def _interface_mode(case_dir: Path) -> bool:
    """Whether the case injects at a downstream sphere (malformed raises)."""
    return interface.load_interface(case_dir) is not None

#: Where a case's wall clock is recorded, relative to the family directory.
COST_MODEL = Path("results") / "cost-model.yaml"


# --------------------------------------------------------------------------- #
# guard
# --------------------------------------------------------------------------- #

def guard(case_dir: Path) -> str | None:
    """Why the case must not run, or ``None`` if it may.

    A case without an explicit suite declaration retains the legacy paper-Kn
    exclusions. Full-decade cases must carry the recognized profile themselves;
    editing only a manifest cannot bypass the case guard.
    """
    case_dir = Path(case_dir)
    try:
        if _interface_mode(case_dir):
            return interface_case.guard(case_dir)
    except (interface.InterfaceConfigError, policy.ExtensionError,
            ValueError) as exc:
        return f"interface.yaml: {exc}"
    cfg = load_case_config(case_dir)
    return policy.direct_exit_problem(cfg, f"{case_dir.resolve().name}/case.yaml")


# --------------------------------------------------------------------------- #
# label
# --------------------------------------------------------------------------- #

def assess(case_dir: Path, cost_model=None, cost_origin=None) -> tuple:
    """``(derived, assessment)`` for the case, from its own ``case.yaml``."""
    if _interface_mode(Path(case_dir)):
        derived = interface_case.derive(Path(case_dir))
        return derived, interface_case.assess(derived, cost_model=cost_model,
                                              cost_origin=cost_origin)
    derived = generation.derive(Path(case_dir))
    return derived, feasibility.assess(*derived, cost_model=cost_model,
                                       cost_origin=cost_origin)


def marker_text(assessment: dict) -> str:
    """The ``NON-PRODUCTION`` marker's content: the verdict and why."""
    return ("This case was run by its runCase.sh job.\n"
            "It is UNSUITABLE FOR PRODUCTION; nothing from it is a validated "
            "result.\n\nUnresolved:\n"
            + "".join(f"  - {line}\n" for line in assessment["unresolved"]))


def label(case_dir: Path) -> str:
    """Mark the case ``NON-PRODUCTION`` if it is unsuitable; return what to print.

    A candidate is left as it is -- including a marker an earlier, unsuitable
    configuration left behind, since results from that run may still be in the
    case.
    """
    case_dir = Path(case_dir)
    if _interface_mode(case_dir):
        text, runnable = interface_case.label(case_dir, MARKER)
        if not runnable:
            raise policy.ExtensionError(text)
        return text
    derived, assessment = assess(case_dir)
    limitations = policy.source_limitations(derived[0])
    source_text = ("Known source limitations (not resolved by numerical feasibility):\n"
                   + "".join(f"  - {line}\n" for line in limitations)) if limitations else ""
    if limitations:
        (case_dir / "KNOWN-LIMITATIONS").write_text(source_text, encoding="utf-8", newline="\n")
    if assessment["status"] == feasibility.STATUS_CANDIDATE:
        return (f"{case_dir.resolve().name}: {feasibility.STATUS_CANDIDATE} -- "
                f"every pre-run criterion is met; not yet a validated result\n" + source_text)
    text = marker_text(assessment) + source_text
    marker = case_dir / MARKER
    marker.write_text(text, encoding="utf-8", newline="\n")
    return text + f"\nmarked {MARKER}: {marker.resolve()}\n"


# --------------------------------------------------------------------------- #
# finished
# --------------------------------------------------------------------------- #

def reached_end(case_dir: Path) -> bool:
    """Whether the case's newest time, reconstructed or not, is at ``endTime``.

    Half a time step of slack, the test ``Time::run()`` applies: ``endTime`` and
    the time directory's name are the same number printed by two formatters.
    """
    case = Path(case_dir)
    control = case / "system" / "controlDict"
    if not control.is_file():
        return False
    entries = {}
    for line in control.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*(endTime|deltaT)\s+([^;]+);", line)
        if match:
            entries[match.group(1)] = float(match.group(2))

    times = []
    for parent in (case, *case.glob("processor[0-9]*")):
        if not parent.is_dir():
            continue
        for entry in parent.iterdir():
            if entry.is_dir():
                try:
                    value = float(entry.name)
                except ValueError:
                    continue
                if value > 0.0:
                    times.append(value)

    latest = max(times, default=None)
    end = entries.get("endTime")
    delta = entries.get("deltaT", 0.0)
    return latest is not None and end is not None and latest >= end - 0.5 * delta


# --------------------------------------------------------------------------- #
# record
# --------------------------------------------------------------------------- #

def family_root(case_dir: Path) -> Path | None:
    """The family directory of a generated case (``<family>/Cases/<case>``)."""
    family = Path(case_dir).resolve().parents[1]
    return family if (family / "study.yaml").is_file() else None


def measurement(case_dir: Path, seconds: float, derived, assessment) -> dict:
    """One ``cost-model.yaml`` entry: the fields the health family records."""
    case_dir = Path(case_dir).resolve()
    if isinstance(derived, interface_case.Derived):
        return {
            "case": case_dir.name,
            "path": case_dir.relative_to(family_root(case_dir)).as_posix(),
            "Kn": float(derived.cfg.exit.knudsen),
            "seconds": float(seconds),
            "parcels": float(derived.parcels),
            "cells": int(derived.mesh_plan.retained_cells),
            "steps": int(derived.run.n_steps),
            "numerical_particle_multiplier": float(
                derived.cfg.dsmc.numerical_particle_multiplier),
            "production_status": assessment["status"],
            "mode": "interface",
        }
    cfg, _, _, plan, run, _, particles, _ = derived
    return {
        "case": case_dir.name,
        "path": case_dir.relative_to(family_root(case_dir)).as_posix(),
        "Kn": float(cfg.exit.knudsen),
        "seconds": float(seconds),
        "parcels": float(particles["total_particles"]),
        "cells": int(plan.n_cells),
        "steps": int(run.n_steps),
        "numerical_particle_multiplier": float(
            cfg.dsmc.numerical_particle_multiplier),
        "production_status": assessment["status"],
    }


@contextlib.contextmanager
def locked(path: Path):
    """An exclusive ``flock`` on ``<path>.lock`` for the duration."""
    import fcntl

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path.with_name(path.name + ".lock"), "a") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def record_entry(path: Path, entry: dict) -> Path:
    """:func:`.feasibility.record_measurement`, under the lock."""
    with locked(path):
        return feasibility.record_measurement(path, entry)


def record(case_dir: Path, seconds: float) -> list[str]:
    """Record the case's wall clock in its family's cost model; return report lines."""
    case_dir = Path(case_dir)
    family = family_root(case_dir)
    if family is None:
        return [f"not timed: {case_dir.resolve()} is not a generated case "
                f"(<family>/Cases/<case>)"]
    health = family.parent / "cai2012-health"
    model, origin = feasibility.load_cost_model(family, health)
    derived, assessment = assess(case_dir, cost_model=model, cost_origin=origin)
    record_entry(family / COST_MODEL, measurement(case_dir, seconds, derived,
                                                  assessment))
    lines = []
    predicted = assessment["cost"].get("predicted_wall_hours")
    if predicted:
        lines.append(f"predicted >= {float(predicted):.2f} h, took "
                     f"{seconds / 3600:.2f} h")
    model, _ = feasibility.load_cost_model(family, health)
    lines.append(f"cost model now: {model.source} ({family / COST_MODEL})")
    return lines


# --------------------------------------------------------------------------- #
# the command line
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.cai2012ext.hpc",
        description="What a case's runCase.sh does around ./Allrun.")
    steps = parser.add_subparsers(dest="step", required=True)
    for name, text in (("guard", "exit 2 if the case is at an excluded Kn"),
                       ("label", "mark the case NON-PRODUCTION if it is "
                                 "unsuitable for production"),
                       ("finished", "exit 0 if the case reached endTime, else 1"),
                       ("record", "record the wall clock in the cost model")):
        step = steps.add_parser(name, help=text)
        step.add_argument("case_dir", type=Path)
        if name == "record":
            step.add_argument("seconds", type=float)
    args = parser.parse_args(argv)

    if not (args.case_dir / "case.yaml").is_file():
        print(f"hpc {args.step}: {args.case_dir.resolve()} has no case.yaml",
              file=sys.stderr)
        return 1
    if args.step == "guard":
        problem = guard(args.case_dir)
        if problem:
            print(f"hpc guard: REFUSED. {problem}", file=sys.stderr)
            return 2
        return 0
    if args.step == "label":
        try:
            print(label(args.case_dir), end="")
        except policy.ExtensionError as exc:
            print(str(exc), end="")
            print("hpc label: REFUSED -- see above", file=sys.stderr)
            return 3
        return 0
    if args.step == "finished":
        return 0 if reached_end(args.case_dir) else 1
    for line in record(args.case_dir, args.seconds):
        print(f"    {line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
