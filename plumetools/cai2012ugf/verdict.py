"""The post-run suitability verdict of a cai2012-unigas case (spec §7).

The dsmcFoam suite judged a case before it ran, by cell/lambda0 and dt/tau0 at
the exit. USP cells fail both by design, so here the verdict is made after the
run, from uniGasFoam's own per-cell output, over the DSMC cells
(``uniGasCollisionModelId == 1``) that hold particles:

* the **effective** collision spacing over the mean free path is at most
  ``checks.max_cell_over_mfp`` (Cai's dx/lambda0 = 1). Not the nominal sub-cell
  size: a sub-cell with one particle draws its partner from the whole cell
  (uniGasFoam ``noTimeCounter.C``), so a cell with N particles gains nothing
  from more than about N/2 sub-cells -- see :func:`effective_spacing_over_mfp`;
* the time step over the mean collision time is at most 0.2;
* a DSMC cell above the occupancy floor (n/n0 >= 1e-3, Cai's lowest contour)
  holds at least ``checks.min_particles_per_cell`` particles;

and, over the case, ``Allconverge.py``'s stationarity verdict is
``stable-window-demonstrated`` and no known source limitation applies.
USP cells are reported (fraction, extent) and never gated.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml

VERDICT_FILE = Path("results") / "unigas_verdict.yaml"
STABLE = "stable-window-demonstrated"
PRODUCTION = "production-candidate"
NON_PRODUCTION = "non-production"
PENDING = {
    "status": "post-run-verdict-pending",
    "production_ready": False,
    "unresolved": ["post-run verdict not computed: run Allverdict.py"],
}
MCT_SLACK = 1.0e-6
#: uniGasFoam writes OpenFOAM's GREAT (1e15) where the mean free path or the
#: collision time is undefined -- no gas, or a collision rate of zero
#: (uniGasVolFields.C). Anything at or above this is "undefined", never a value.
SENTINEL = 1.0e10


def effective_spacing_over_mfp(extents, sub_cell_levels, particles, mfp, *,
                               min_particles_per_sub_cell: float = 2.0) -> np.ndarray:
    """Effective collision-partner spacing over the mean free path, per cell.

    ``S = sx sy sz``; ``f = min(1, N / (p S))^(1/3)`` with ``p`` the particles a
    sub-cell needs to supply its own partners; the usable divisions per
    direction are ``max(1, s_d f)``; the spacing is the largest
    ``extent_d / divisions_d``. NaN where the mean free path is not a finite
    positive number.
    """
    extents = np.asarray(extents, dtype=float).reshape(-1, 3)
    levels = np.maximum(np.asarray(sub_cell_levels, dtype=float).reshape(-1, 3), 1.0)
    particles = np.maximum(np.asarray(particles, dtype=float).reshape(-1), 0.0)
    mfp = np.asarray(mfp, dtype=float).reshape(-1)
    n_sub = levels.prod(axis=1)
    fill = np.minimum(1.0, particles / (float(min_particles_per_sub_cell) * n_sub))
    divisions = np.maximum(1.0, levels * np.cbrt(fill)[:, None])
    spacing = (extents / divisions).max(axis=1)
    out = np.full(spacing.shape, np.nan)
    good = np.isfinite(mfp) & (mfp > 0.0) & (mfp < SENTINEL)
    out[good] = spacing[good] / mfp[good]
    return out


def _stats(values: np.ndarray, limit: float | None = None) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"n": 0, "max": None, "p50": None, "p95": None, "p99": None,
                "fraction_over_limit": None}
    doc = {"n": int(values.size), "max": float(values.max()),
           "p50": float(np.percentile(values, 50)),
           "p95": float(np.percentile(values, 95)),
           "p99": float(np.percentile(values, 99))}
    doc["fraction_over_limit"] = (None if limit is None
                                  else float(np.mean(values > limit)))
    return doc


def assess_fields(fields: dict, *, n0: float, settings, limits: dict,
                  stationarity_verdict: str | None, known_limitations: list) -> dict:
    """The verdict from per-cell arrays. Pure: no files are read."""
    ids = np.rint(np.asarray(fields["collision_model_id"], dtype=float))
    particles = np.asarray(fields["particles"], dtype=float)
    density = np.asarray(fields["number_density"], dtype=float)
    n_cells = ids.size
    dsmc = ids == 1.0
    usp = ids == 0.0
    occupied = particles > 0.0
    judged = dsmc & occupied

    spacing_limit = float(limits["max_effective_spacing_over_mfp"])
    occupancy_min = float(limits["min_particles_per_dsmc_cell"])
    mct_limit = float(settings.max_time_step_mct_ratio)

    effective = effective_spacing_over_mfp(
        fields["extents"], fields["sub_cell_levels"], particles, fields["mfp"],
        min_particles_per_sub_cell=settings.min_particles_per_sub_cell)
    mct = np.asarray(fields["time_step_mct_ratio"], dtype=float)
    mct_defined = np.isfinite(mct) & (mct < SENTINEL)
    floor = density / float(n0) >= float(settings.occupancy_floor_n_over_n0)
    counted = dsmc & floor

    spacing_stats = _stats(effective[judged], spacing_limit)
    mct_stats = _stats(mct[judged & mct_defined], mct_limit)
    occupancy_stats = _stats(particles[counted])
    below = int(np.sum(particles[counted] < occupancy_min))

    unresolved = []
    if spacing_stats["max"] is not None and spacing_stats["max"] > spacing_limit:
        unresolved.append(
            f"collision spacing: effective spacing/mfp reaches {spacing_stats['max']:.4g} "
            f"in DSMC cells (limit {spacing_limit:g}; "
            f"{spacing_stats['fraction_over_limit']:.3%} of them over it)")
    if mct_stats["max"] is not None and mct_stats["max"] > mct_limit * (1.0 + MCT_SLACK):
        unresolved.append(
            f"time step: dt/tau reaches {mct_stats['max']:.4g} in DSMC cells "
            f"(limit {mct_limit:g})")
    if below:
        unresolved.append(
            f"particles: {below} DSMC cell(s) above n/n0 = "
            f"{settings.occupancy_floor_n_over_n0:g} hold fewer than {occupancy_min:g} particles")
    if stationarity_verdict != STABLE:
        unresolved.append(f"stationarity not demonstrated ({stationarity_verdict or 'unknown'})")
    unresolved.extend(str(item) for item in known_limitations)

    centres = np.asarray(fields["centres"], dtype=float).reshape(-1, 3)
    x_usp = centres[usp, 0]
    kn_gll = fields.get("kn_gll")
    ready = not unresolved
    return {
        "status": PRODUCTION if ready else NON_PRODUCTION,
        "production_ready": ready,
        "unresolved": unresolved,
        "n_cells": int(n_cells),
        "dsmc": {
            "n_cells": int(dsmc.sum()),
            "fraction": float(dsmc.mean()) if n_cells else 0.0,
            "excluded_unoccupied": int(np.sum(dsmc & ~occupied)),
            "undefined_mean_free_path": int(np.sum(judged & ~np.isfinite(effective))),
            "undefined_time_step_ratio": int(np.sum(judged & ~mct_defined)),
            "effective_spacing_over_mfp": spacing_stats,
            "time_step_mct_ratio": mct_stats,
            "particles_above_floor": occupancy_stats,
            "cells_below_particle_minimum": below,
        },
        "usp": {
            "n_cells": int(usp.sum()),
            "fraction": float(usp.mean()) if n_cells else 0.0,
            "x_range_m": ([float(x_usp.min()), float(x_usp.max())] if x_usp.size else None),
        },
        "kn_gll": None if kn_gll is None else _stats(np.asarray(kn_gll, dtype=float)),
        "stationarity": stationarity_verdict,
        "known_limitations": list(known_limitations),
        "limits": {"max_effective_spacing_over_mfp": spacing_limit,
                   "max_time_step_mct_ratio": mct_limit,
                   "min_particles_per_dsmc_cell": occupancy_min,
                   "occupancy_floor_n_over_n0": float(settings.occupancy_floor_n_over_n0)},
    }


# --------------------------------------------------------------------------- #
# reading a case
# --------------------------------------------------------------------------- #

def mesh_extents(case_dir: Path) -> np.ndarray:
    """``(n_cells, 3)`` bounding-box sizes of every cell [m]."""
    from plumetools.mesh.polymesh import read_faces_flat, read_points_array
    from plumetools.rarefaction.foam import _mesh_topology, cell_extents

    owner, neighbour, n_cells = _mesh_topology(case_dir)
    points = read_points_array(case_dir)
    offsets, labels = read_faces_flat(case_dir)
    return cell_extents(points, offsets, labels, owner, neighbour, n_cells)


def _read(time_dir: Path, name: str, n_cells: int, required: bool = True):
    from plumetools.foamio.volfield import read_internal

    path = time_dir / name
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"no {path}")
        return None
    return np.asarray(read_internal(path).over(n_cells, str(time_dir)), dtype=float)


def _stationarity(case_dir: Path) -> str | None:
    path = case_dir / "results" / "stationarity.yaml"
    if not path.is_file():
        return None
    document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return (document.get("assessment") or {}).get("verdict")


def assess_run(case_dir: Path, *, time: str | None = None) -> dict:
    """Assess a reconstructed, converted case and write ``results/unigas_verdict.yaml``."""
    from plumetools.cai2012ext import direct_case, policy
    from plumetools.cai2012ugf.convert import time_directories
    from plumetools.cai2012ugf.settings import load_settings

    case_dir = Path(case_dir)
    case = direct_case.derive(case_dir)
    settings = load_settings(case_dir)
    if time is None:
        times = time_directories(case_dir)
        if not times:
            raise FileNotFoundError(f"{case_dir} has no reconstructed time directory")
        time = times[-1].name
    time_dir = case_dir / time
    extents = mesh_extents(case_dir)
    n = extents.shape[0]
    f = settings.fields.mean_field
    particles = _read(time_dir, "dsmcRhoNMean", n, required=False)
    if particles is None:
        rho_n = _read(time_dir, f"rhoN_{f}", n)
        particles = np.where(rho_n > 0, _read(time_dir, f"uniGasRhoNMean_{f}", n), 0.0)
    fields = {
        "collision_model_id": _read(time_dir, "uniGasCollisionModelId", n),
        "sub_cell_levels": _read(time_dir, "uniGasSubCellLevels", n),
        "particles": particles,
        "mfp": _read(time_dir, f"variableHardSphereMeanFreePath_{f}", n),
        "extents": extents,
        "time_step_mct_ratio": _read(time_dir, f"timeStepMCTRatio_{f}", n),
        "kn_gll": _read(time_dir, "KnGLL", n, required=False),
        "number_density": _read(time_dir, f"rhoN_{f}", n),
        "centres": _read(time_dir, "C", n),
    }
    v = settings.verdict
    limits = {
        "max_effective_spacing_over_mfp": (v.max_effective_spacing_over_mfp
                                           if v.max_effective_spacing_over_mfp is not None
                                           else case.cfg.checks.max_cell_over_mfp),
        "min_particles_per_dsmc_cell": (v.min_particles_per_dsmc_cell
                                        if v.min_particles_per_dsmc_cell is not None
                                        else case.cfg.checks.min_particles_per_cell),
    }
    verdict = assess_fields(
        fields, n0=case.exit_state.number_density_per_m3, settings=v, limits=limits,
        stationarity_verdict=_stationarity(case_dir),
        known_limitations=policy.source_limitations(case.cfg))
    verdict["time"] = time
    verdict["diameter_m"] = float(case.cfg.nozzle.diameter_m)
    path = case_dir / VERDICT_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# GENERATED by plumetools.cai2012ugf.verdict -- the post-run "
                    "suitability verdict (spec §7).\n"
                    + yaml.safe_dump(verdict, sort_keys=False), encoding="utf-8")
    return verdict


def read_verdict(case_dir: Path) -> dict:
    """The written verdict, or the pending one when none has been written."""
    path = Path(case_dir) / VERDICT_FILE
    if not path.is_file():
        return dict(PENDING)
    document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not all(k in document for k in ("status", "production_ready", "unresolved")):
        return dict(PENDING)
    return document


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.cai2012ugf.verdict",
        description="Write a cai2012-unigas case's post-run suitability verdict.")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--time", default=None)
    args = parser.parse_args(argv)
    try:
        verdict = assess_run(args.case_dir, time=args.time)
    except (FileNotFoundError, ValueError) as exc:
        print(f"verdict: {exc}", file=sys.stderr)
        return 1
    print(f"  {verdict['status'].upper()} at t = {verdict['time']}")
    for reason in verdict["unresolved"]:
        print(f"    - {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
