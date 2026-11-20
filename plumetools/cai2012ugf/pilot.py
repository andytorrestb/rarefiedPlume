"""The pilot report: is 3-D production on uniGasFoam affordable? (spec §8)

Reads a pilot case (``Pilot/<case>``, run for a few transits on its production
mesh) and reports what the go/no-go decision needs:

* the time-step history, and whether DSMC collision times or Courant set it;
* the DSMC/USP split and where each lies;
* the effective collision spacing over the mean free path in DSMC cells, and
  how many cells sit at the sub-cell cap;
* the particle count, the wall clock per step, the rank imbalance;
* the disk used by the newest write;
* a projected wall time for the full run (startup and sampling).

    python -m plumetools.cai2012ugf.pilot <case_dir>

writes ``results/pilot_report.yaml`` and ``results/pilot_report.md``. The
decision is the owner's; this only measures.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import yaml

_TIME = re.compile(r"^Time = (\S+)\s*$")
_DT = re.compile(r"^\s*Time step\s*=\s*(\S+)")
_PARTICLES = re.compile(r"^\s*Number of particles\s*=\s*(\S+)")
_IMBALANCE = re.compile(r"^\s*Maximum imbalance\s*=\s*(\S+)")
_CLOCK = re.compile(r"ClockTime = (\S+) s")
KEYS = ("time", "dt", "particles", "imbalance", "clock")


def parse_log(path: Path) -> dict:
    """One record per ``Time =`` block of a (possibly appended) uniGasFoam log."""
    records = []
    current = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        match = _TIME.match(line)
        if match:
            try:
                current = {"time": float(match.group(1))}
            except ValueError:
                current = None
                continue
            records.append(current)
            continue
        if current is None:
            continue
        for key, pattern in (("dt", _DT), ("particles", _PARTICLES),
                             ("imbalance", _IMBALANCE), ("clock", _CLOCK)):
            found = pattern.search(line)
            if found and key not in current:
                try:
                    current[key] = float(found.group(1))
                except ValueError:
                    pass
    out = {key: np.array([r.get(key, np.nan) for r in records], dtype=float) for key in KEYS}
    out["particles"] = np.array([int(v) if np.isfinite(v) else -1 for v in out["particles"]])
    return out


def projected_hours(t_now: float, end: float, recent_clock_per_sim_s: float,
                    elapsed_s: float) -> float:
    """Elapsed plus the remaining simulated time at the recent rate, in hours."""
    remaining = max(0.0, float(end) - float(t_now))
    return (float(elapsed_s) + remaining * float(recent_clock_per_sim_s)) / 3600.0


def time_step_setter(max_dt_over_tau: float | None, limit: float) -> str:
    """What the measured time step answers to, from the largest DSMC dt/tau.

    Near the limit, DSMC collision times set the step. Well below it, Courant
    does. ABOVE it, neither is in control: the limit is not being met -- e.g.
    the solver has not yet adapted (it adapts every ``adaptation_interval``
    steps, and only after DSMC cells exist).
    """
    if max_dt_over_tau is None:
        return "Courant (no DSMC cells)"
    if max_dt_over_tau > 1.05 * limit:
        return (f"NEITHER: DSMC dt/tau reaches {max_dt_over_tau:.3g}, above the "
                f"{limit:g} limit -- the time step has not been adapted to it")
    if max_dt_over_tau >= 0.9 * limit:
        return "DSMC collision time"
    return "Courant"


def _recent_rate(log: dict) -> tuple:
    """Wall seconds per simulated second over the last 10 % of records, and elapsed."""
    time, clock = log["time"], log["clock"]
    good = np.isfinite(time) & np.isfinite(clock)
    time, clock = time[good], clock[good]
    if time.size < 2:
        return None, float(clock[-1]) if clock.size else 0.0
    # The clock restarts at every resubmission; use the last monotonic stretch.
    restart = np.where(np.diff(clock) < 0)[0]
    start = int(restart[-1]) + 1 if restart.size else 0
    time, clock = time[start:], clock[start:]
    elapsed = float(clock[-1])
    if time.size < 2:
        return None, elapsed
    k = max(0, int(0.9 * (time.size - 1)))
    span = float(time[-1] - time[k])
    return (float(clock[-1] - clock[k]) / span if span > 0 else None), elapsed


def _bytes(path: Path) -> int:
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())


def report(case_dir: Path) -> dict:
    """Measure the pilot and write ``results/pilot_report.{yaml,md}``."""
    from plumetools.cai2012ext import direct_case
    from plumetools.cai2012ugf.convert import time_directories
    from plumetools.cai2012ugf.settings import load_settings
    from plumetools.cai2012ugf.verdict import (SENTINEL, _read, effective_spacing_over_mfp,
                                               mesh_extents)

    case_dir = Path(case_dir)
    case = direct_case.derive(case_dir)
    settings = load_settings(case_dir)
    run = case.run
    log = parse_log(case_dir / "log.uniGasFoam")
    rate, elapsed = _recent_rate(log)
    t_now = float(np.nanmax(log["time"])) if log["time"].size else 0.0
    dt = log["dt"][np.isfinite(log["dt"])]
    doc = {
        "case": case_dir.name,
        "Kn": float(case.cfg.exit.knudsen),
        "simulated_s": t_now,
        "domain_transits": t_now / run.domain_transit_s,
        "time_step": {"first": float(dt[0]) if dt.size else None,
                      "last": float(dt[-1]) if dt.size else None,
                      "min": float(dt.min()) if dt.size else None,
                      "dsmcfoam_suite_delta_t_s": float(run.delta_t_s)},
        "particles": {"last": int(log["particles"][-1]) if log["particles"].size else None,
                      "max": int(log["particles"].max()) if log["particles"].size else None},
        "imbalance_percent": {
            "max": float(np.nanmax(log["imbalance"])) if np.isfinite(log["imbalance"]).any() else None,
            "p95": (float(np.nanpercentile(log["imbalance"], 95))
                    if np.isfinite(log["imbalance"]).any() else None)},
        "wall_s_per_simulated_s": rate,
        "elapsed_wall_s": elapsed,
        "full_run_end_s": float(run.end_time_s),
        "projected_wall_hours": (None if rate is None else
                                 projected_hours(t_now, run.end_time_s, rate, elapsed)),
    }

    times = time_directories(case_dir)
    if times:
        t = times[-1]
        extents = mesh_extents(case_dir)
        n = extents.shape[0]
        f = settings.fields.mean_field
        ids = np.rint(_read(t, "uniGasCollisionModelId", n))
        levels = _read(t, "uniGasSubCellLevels", n)
        rho = _read(t, f"rhoN_{f}", n)
        counts = _read(t, f"uniGasRhoNMean_{f}", n)
        counts = np.where(rho > 0, counts, 0.0)
        mfp = _read(t, f"variableHardSphereMeanFreePath_{f}", n)
        mct = _read(t, f"timeStepMCTRatio_{f}", n)
        x = _read(t, "C", n)[:, 0] if (t / "C").is_file() else None
        dsmc, usp = ids == 1, ids == 0
        occupied = dsmc & (counts > 0)
        eff = effective_spacing_over_mfp(extents, levels, counts, mfp,
                                         min_particles_per_sub_cell=settings.verdict.min_particles_per_sub_cell)
        mct_ok = occupied & np.isfinite(mct) & (mct < SENTINEL)
        diameter = float(case.cfg.nozzle.diameter_m)

        def extent(mask):
            if x is None or not mask.any():
                return None
            return [float(x[mask].min() / diameter), float(x[mask].max() / diameter)]

        finite = eff[occupied & np.isfinite(eff)]
        doc["fields"] = {
            "time": t.name,
            "dsmc_fraction": float(dsmc.mean()), "usp_fraction": float(usp.mean()),
            "dsmc_x_over_D": extent(dsmc), "usp_x_over_D": extent(usp),
            "effective_spacing_over_mfp": (
                {"max": float(finite.max()), "p50": float(np.percentile(finite, 50)),
                 "p95": float(np.percentile(finite, 95)),
                 "fraction_over_1": float(np.mean(finite > 1.0))} if finite.size else None),
            "cells_at_sub_cell_cap": int(np.sum(levels.max(axis=1) >= settings.adaptive.max_sub_cell_levels)),
            "max_dt_over_tau_dsmc": float(mct[mct_ok].max()) if mct_ok.any() else None,
            "time_step_set_by": time_step_setter(
                float(mct[mct_ok].max()) if mct_ok.any() else None,
                settings.adaptive.max_time_step_mct_ratio),
        }
        newest = [t] + [p / t.name for p in sorted(case_dir.glob("processor*")) if (p / t.name).is_dir()]
        doc["disk_bytes_newest_write"] = int(sum(_bytes(d) for d in newest))

    results = case_dir / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "pilot_report.yaml").write_text(
        "# GENERATED by plumetools.cai2012ugf.pilot -- measurements for the owner's go/no-go.\n"
        + yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
    (results / "pilot_report.md").write_text(_markdown(doc), encoding="utf-8")
    return doc


def _markdown(doc: dict) -> str:
    fields = doc.get("fields") or {}
    hours = doc.get("projected_wall_hours")
    lines = [
        f"# Pilot report: {doc['case']} (Kn = {doc['Kn']:g})", "",
        "Measurements only; going ahead is the owner's decision (spec §8). Not validated.", "",
        f"* simulated: {doc['simulated_s']:.4g} s ({doc['domain_transits']:.2f} domain transits)",
        f"* time step: last {doc['time_step']['last']}, min {doc['time_step']['min']} "
        f"(dsmcFoam suite: {doc['time_step']['dsmcfoam_suite_delta_t_s']:.4g} s); "
        f"set by {fields.get('time_step_set_by', 'unknown')}",
        f"* particles: {doc['particles']['last']} (max {doc['particles']['max']})",
        f"* rank imbalance: max {doc['imbalance_percent']['max']} %",
        f"* DSMC / USP cells: {fields.get('dsmc_fraction')} / {fields.get('usp_fraction')}; "
        f"USP at x/D {fields.get('usp_x_over_D')}",
        f"* effective collision spacing / mfp in DSMC cells: {fields.get('effective_spacing_over_mfp')}",
        f"* cells at the sub-cell cap: {fields.get('cells_at_sub_cell_cap')}",
        f"* newest write on disk: {doc.get('disk_bytes_newest_write')} bytes",
        f"* projected wall time for the full run: "
        + (f"{hours:.1f} h" if hours is not None else "unknown (too few log records)"),
        "",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m plumetools.cai2012ugf.pilot")
    parser.add_argument("case_dir", type=Path)
    args = parser.parse_args(argv)
    try:
        doc = report(args.case_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"pilot: {exc}", file=sys.stderr)
        return 1
    print(_markdown(doc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
