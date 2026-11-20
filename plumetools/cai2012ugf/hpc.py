"""Run legs, cloud-safe restarts, build checks and records for one case.

``Allrun`` runs a case in two legs (spec §6): ``startup`` to
``average_start_s``, then ``sampling`` to ``end_time_s``, each with its own
``controlDict``. :func:`resume_plan` says which leg comes next and whether the
case may be continued at all.

Why the frozen restart plan is not enough
-----------------------------------------
:func:`plumetools.foamio.checkpoints.plan_restart` finds the newest time every
rank wrote completely. uniGasFoam, however, keeps only the newest particle
cloud: after each write it deletes the previous write's ``lagrangian/uniGas``
(``uniGasCloud::cleanLagrangian``). A job killed mid-write therefore leaves
ranks that finished the new write -- and already deleted their older cloud --
beside ranks that did not. The planner refuses the partial time and advises
moving it aside; resuming from the older time would then find no cloud on
the first kind of rank, and uniGasFoam, finding no ``lagrangian/uniGas``,
silently re-fills that subdomain from ``uniGasInitialisationDict``.

A rank that holds no particle legitimately writes no cloud (OpenFOAM's
``Cloud::writeObject`` writes only when ``size() > 0``), and with the
near-vacuum start many ranks hold none for a while. So the test is the rank's
own density: a rank whose ``rhoN_<field>`` is positive at the resume time but
which has no cloud there has lost it, and the case needs a manual restart.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from plumetools.foamio.checkpoints import RestartPlan, plan_restart, rank_directories, scan

CLOUD = "uniGas"
#: Relative. Time directories are named to timePrecision 8 significant digits,
#: so the write that ends a leg can be named up to ~5e-9 relative either side of
#: its endTime; 1e-6 is far above that and far below one time step.
TOLERANCE = 1.0e-6
COST_MODEL = Path("results") / "cost-model.yaml"
PROVENANCE = Path("results") / "provenance-unigas.yaml"
COMMIT_FILE = "uniGasFoam.commit"
SOLVER = "uniGasFoam hybrid USP-DSMC"
HEADER_KEYS = {"version", "format", "class", "location", "object", "note"}


def leg_for(time_s: float | None, average_start_s: float, end_time_s: float,
            slack_s: float = 0.0) -> str:
    """``startup``, ``sampling`` or ``done`` for a case whose newest time is ``time_s``.

    ``slack_s`` is how far short of a leg's end its last write may fall: OpenFOAM
    writes within half a step of a write time and stops within half a step of
    ``endTime``, so pass about half the step that write was taken with.
    """
    if time_s is None:
        return "startup"
    if time_s >= end_time_s * (1.0 - TOLERANCE) - slack_s:
        return "done"
    if time_s >= average_start_s * (1.0 - TOLERANCE) - slack_s:
        return "sampling"
    return "startup"


_DELTA_T = re.compile(r"^\s*deltaT\s+([^;\s]+)\s*;", re.MULTILINE)


def step_at(case_dir: Path, time: str) -> float:
    """The time step recorded in ``<time>/uniform/time`` (a rank's or the case's), or 0."""
    case_dir = Path(case_dir)
    for base in [*rank_directories(case_dir)[:1], case_dir]:
        path = base / time / "uniform" / "time"
        if path.is_file():
            match = _DELTA_T.search(path.read_text(encoding="utf-8", errors="replace"))
            if match:
                try:
                    return float(match.group(1))
                except ValueError:
                    return 0.0
    return 0.0


def _slack(case_dir: Path, time: str | None) -> float:
    return 0.51 * step_at(case_dir, time) if time else 0.0


@dataclass
class UgfRestartPlan:
    base: RestartPlan
    leg: str
    refusal: str | None

    def shell(self) -> str:
        """``KEY=value`` lines for ``Allrun`` to ``eval``."""
        return (self.base.shell() + f"LEG={self.leg}\n"
                + f"REFUSAL={1 if self.refusal else 0}\n")


def _derive(case_dir: Path):
    from plumetools.cai2012ext import direct_case
    return direct_case.derive(case_dir)


def _mean_field(case_dir: Path) -> str:
    from plumetools.cai2012ugf.settings import load_settings
    return load_settings(case_dir).fields.mean_field


def _held_gas(time_dir: Path, field: str) -> bool:
    from plumetools.foamio.volfield import read_internal
    path = time_dir / f"rhoN_{field}"
    if not path.is_file():
        return False
    return bool(np.any(np.asarray(read_internal(path).values, dtype=float) > 0.0))


def lost_clouds(case_dir: Path, time: str, *, decomposed: bool) -> list[str]:
    """Ranks (or the case directory) that held gas at ``time`` but have no cloud."""
    case_dir = Path(case_dir)
    field = _mean_field(case_dir)
    bases = rank_directories(case_dir) if decomposed else [case_dir]
    lost = []
    for base in bases:
        t = base / time
        if t.is_dir() and not (t / "lagrangian" / CLOUD).is_dir() and _held_gas(t, field):
            lost.append(base.name if base != case_dir else "the case directory")
    return lost


def cloud_complete(case_dir: Path, time: str, *, decomposed: bool = True) -> bool:
    """No rank that held gas at ``time`` has lost its cloud there."""
    return not lost_clouds(case_dir, time, decomposed=decomposed)


def resume_plan(case_dir: Path, n_subdomains: int) -> UgfRestartPlan:
    """The frozen restart plan, the next leg, and the cloud-loss refusal."""
    case_dir = Path(case_dir)
    base = plan_restart(case_dir, n_subdomains)
    run = _derive(case_dir).run
    resume = float(base.resume) if base.resume else None
    leg = leg_for(resume, run.average_start_s, run.end_time_s, _slack(case_dir, base.resume))
    reasons = [base.refusal] if base.refusal else []
    if base.resume:
        decomposed = n_subdomains > 1 and base.proc_time == base.resume
        lost = lost_clouds(case_dir, base.resume, decomposed=decomposed)
        if lost:
            reasons.append(
                f"t = {base.resume}: {', '.join(lost)} held gas (rhoN_{_mean_field(case_dir)} > 0) "
                f"but {'have' if len(lost) > 1 else 'has'} no lagrangian/{CLOUD} cloud there. "
                f"uniGasFoam keeps only the newest cloud, and resuming would silently "
                f"re-fill those subdomains from uniGasInitialisationDict.\n"
                f"This needs a manual restart: find the newest time at which every rank "
                f"that held gas still has lagrangian/{CLOUD}, move every later time aside "
                f"on every rank ('mv processorN/<t> processorN/<t>.lost'), and resubmit. "
                f"If there is none, the case must start again (./Allclean). Nothing was changed.")
    return UgfRestartPlan(base=base, leg=leg, refusal="\n".join(reasons) or None)


def finished(case_dir: Path) -> bool:
    """The newest complete time has reached the end of the sampling leg."""
    state = scan(case_dir)
    latest = state.latest
    run = _derive(case_dir).run
    return leg_for(latest[0] if latest else None, run.average_start_s, run.end_time_s,
                   _slack(case_dir, latest[1] if latest else None)) == "done"


# --------------------------------------------------------------------------- #
# the build
# --------------------------------------------------------------------------- #

def _expected_commit(case_dir: Path) -> str | None:
    from plumetools.cai2012.config import load_case_config
    value = load_case_config(case_dir).meta.get("unigas_commit")
    return None if value in (None, "", "pending-task-1") else str(value)


def _built_commit() -> str | None:
    bin_dir = os.environ.get("FOAM_USER_APPBIN")
    if not bin_dir:
        return None
    path = Path(bin_dir) / COMMIT_FILE
    return path.read_text(encoding="utf-8").strip() if path.is_file() else None


def check_build(case_dir: Path) -> str | None:
    """Why the built uniGasFoam is not the one this case expects, or ``None``."""
    expected = _expected_commit(case_dir)
    if expected is None:
        return (f"{case_dir}/case.yaml has no meta.unigas_commit: set it in study.yaml "
                f"and regenerate")
    built = _built_commit()
    if built is None:
        return (f"no $FOAM_USER_APPBIN/{COMMIT_FILE}: build uniGasFoam with runBuild.sh, "
                f"which records the commit it built")
    if not (built.startswith(expected) or expected.startswith(built)):
        return f"uniGasFoam was built from {built}, but this case expects {expected}"
    return None


# --------------------------------------------------------------------------- #
# records
# --------------------------------------------------------------------------- #

def family_root(case_dir: Path) -> Path:
    """``<family>/Cases/<case>`` or ``<family>/Pilot/<case>`` -> ``<family>``."""
    return Path(case_dir).resolve().parent.parent


def record(case_dir: Path, seconds: float) -> Path:
    """Append this run's wall clock to the family's cost model."""
    from plumetools.cai2012ext.hpc import locked
    from plumetools.cai2012ugf.verdict import read_verdict

    case_dir = Path(case_dir).resolve()
    derived = _derive(case_dir)
    root = family_root(case_dir)
    path = root / COST_MODEL
    entry = {
        "case": case_dir.name,
        "path": str(case_dir.relative_to(root)).replace("\\", "/"),
        "Kn": float(derived.cfg.exit.knudsen),
        "seconds": float(seconds),
        "cells": int(derived.plan.n_cells),
        "solver": derived.cfg.meta.get("solver", SOLVER),
        "unigas_commit": derived.cfg.meta.get("unigas_commit"),
        "production_status": read_verdict(case_dir)["status"],
        "recorded_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with locked(path):                       # flock on <path>.lock
        document = (yaml.safe_load(path.read_text(encoding="utf-8")) if path.is_file() else None) or {}
        document.setdefault("measurements", []).append(entry)
        path.write_text("# GENERATED by plumetools.cai2012ugf.hpc record -- do not edit.\n"
                        "# Measured wall clock per uniGasFoam run.\n\n"
                        + yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def dictionary_entries(path: Path) -> dict:
    """``key value;`` entries of an OpenFOAM dictionary, as strings (header excluded)."""
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    text = re.sub(r"//[^\n]*", "", text)
    text = re.sub(r"FoamFile\s*\{[^}]*\}", "", text, count=1)
    entries = {}
    for key, value in re.findall(r"^\s*(\w+)\s+([^;{}\n]+);", text, re.MULTILINE):
        if key not in HEADER_KEYS and key not in entries:
            entries[key] = value.strip()
    return entries


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_provenance(case_dir: Path) -> Path:
    """Write ``results/provenance-unigas.yaml`` (spec §11)."""
    from plumetools.cai2012ext.provenance import git_provenance
    from plumetools.cai2012ugf import convert
    from plumetools.cai2012ugf.policy import repo_root
    from plumetools.cai2012ugf.settings import as_dict, load_settings

    case_dir = Path(case_dir)
    derived = _derive(case_dir)
    meta = derived.cfg.meta
    document = {
        "solver": meta.get("solver", SOLVER),
        "unigas_commit": meta.get("unigas_commit"),
        "unigas_branch": meta.get("unigas_branch"),
        "built_commit": _built_commit(),
        "openfoam_version": os.environ.get("WM_PROJECT_VERSION"),
        "plumetools": git_provenance(repo_root()),
        "converter": {"module": convert.GENERATOR, "sha256": _sha256(Path(convert.__file__))},
        "unigas_settings": as_dict(load_settings(case_dir)),
        "uniGasProperties": dictionary_entries(case_dir / "constant" / "uniGasProperties"),
        "hybridDecompositionDict": dictionary_entries(
            case_dir / "system" / "hybridDecompositionDict"),
        "initial_state": {
            "protocol": "uniGasMeshFill near-vacuum at rest (the dsmcFoam suite's)",
            "number_density_per_m3": float(derived.cfg.dsmc.initial_number_density_per_m3),
            "temperature_K": float(derived.cfg.dsmc.initial_temperature_K),
        },
        "written_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }
    path = case_dir / PROVENANCE
    if path.is_file():
        previous = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for key in ("prerun_verdict", "postrun_verdict"):
            if key in previous:
                document[key] = previous[key]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# GENERATED by plumetools.cai2012ugf.hpc provenance -- do not edit.\n"
                    + yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


def add_verdicts(case_dir: Path, prerun: dict | None, postrun: dict | None) -> Path:
    """Record the pre- and post-run verdicts in the provenance file."""
    path = Path(case_dir) / PROVENANCE
    document = (yaml.safe_load(path.read_text(encoding="utf-8")) if path.is_file() else None) or {}
    if prerun is not None:
        document["prerun_verdict"] = prerun
    if postrun is not None:
        document["postrun_verdict"] = postrun
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# GENERATED by plumetools.cai2012ugf.hpc provenance -- do not edit.\n"
                    + yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# the command line
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m plumetools.cai2012ugf.hpc")
    steps = parser.add_subparsers(dest="step", required=True)
    plan = steps.add_parser("resume-plan")
    plan.add_argument("case_dir", type=Path)
    plan.add_argument("--np", type=int, required=True)
    for name in ("finished", "check-build", "provenance"):
        steps.add_parser(name).add_argument("case_dir", type=Path)
    rec = steps.add_parser("record")
    rec.add_argument("case_dir", type=Path)
    rec.add_argument("seconds", type=float)
    args = parser.parse_args(argv)

    if args.step == "resume-plan":
        result = resume_plan(args.case_dir, args.np)
        sys.stdout.write(result.shell())
        if result.refusal:
            print(result.refusal, file=sys.stderr)
            return 3
        return 0
    if args.step == "finished":
        return 0 if finished(args.case_dir) else 1
    if args.step == "check-build":
        problem = check_build(args.case_dir)
        if problem:
            print(f"check-build: {problem}", file=sys.stderr)
            return 1
        return 0
    if args.step == "provenance":
        print(write_provenance(args.case_dir))
        return 0
    print(record(args.case_dir, args.seconds))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
