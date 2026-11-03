"""Resource and resolution assessment, before anything is launched.

Two estimates, both explicit about their inputs:

* :func:`assess_configuration` -- before the CFD exists. The DSMC density is the
  **collisionless** solution (Cai's disk integral) at ``T = Te``. For a
  collisional case that is a shape estimate, labelled as such.
* :func:`dsmc_estimate` -- after the CFD and the interface. Uses the CFD state on
  every cell the DSMC will actually own.

What is counted
---------------
* cells and parcels, with the weight ``W`` set so the densest DSMC cell holds
  ``dsmc.parcels_per_cell_target`` parcels (standard dsmcFoam has one global
  weight);
* memory = parcels x ``bytes_per_parcel`` + cells x ``bytes_per_cell`` (both
  configurable; refine them from a measured run with :func:`measured_rates`);
* the time step = min(Courant on ``|U| + 3 sigma``, a fraction of the smallest
  nominal collision time), the number of steps for startup + sampling, and the
  runtime at ``seconds_per_parcel_step``;
* collision-cell size over local mean free path, and step over collision time.

A dense full-DSMC case is never called affordable because its mesh fits in RAM:
the parcel storage, the step count and the collision work are all here, and a
smaller weight reduces noise but does not repair ``dx/lambda > 1``.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

from .config import CaseConfig
from .mesh import MeshPlan
from .physics import SourceState

GB = 1024.0 ** 3


def _lattice_sample(plan: MeshPlan, max_points: int = 60_000):
    """Cell centres and volumes of the box lattice, sub-sampled for big meshes.

    Returns ``(centres, volumes, sizes, weight)`` where ``weight`` is how many
    real cells each sample stands for (1 when not sub-sampled).
    """
    edges = [plan.edges(a) for a in range(3)]
    n = [len(e) - 1 for e in edges]
    stride = max(1, int(math.ceil((n[0] * n[1] * n[2] / max_points) ** (1 / 3))))
    idx = [np.arange(0, k, stride) for k in n]
    cen = [0.5 * (e[:-1] + e[1:])[i] for e, i in zip(edges, idx)]
    wid = [np.diff(e)[i] for e, i in zip(edges, idx)]
    gx, gy, gz = np.meshgrid(*cen, indexing="ij")
    wx, wy, wz = np.meshgrid(*wid, indexing="ij")
    centres = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    sizes = np.stack([wx.ravel(), wy.ravel(), wz.ravel()], axis=1)
    volumes = sizes.prod(axis=1)
    weight = (n[0] * n[1] * n[2]) / centres.shape[0]
    return centres, volumes, sizes, weight


def cfd_estimate(cfg: CaseConfig, source: SourceState, plan: MeshPlan) -> dict:
    from .cfd import timing
    tm = timing(cfg, source, plan)
    dt = min(tm["convective_dt_estimate_s"], tm["max_delta_t_s"])
    steps = tm["end_time_s"] / dt
    cost = float(cfg.get("resources.seconds_per_cfd_cell_step", 3e-6))
    ranks = int(cfg.get("resources.max_ranks", 1))
    return {
        "cells": plan.n_cells,
        "memory_GB": plan.n_cells * 4000.0 / GB,
        "dt_estimate_s": dt, "steps_estimate": steps,
        "runtime_estimate_s_at_max_ranks": plan.n_cells * steps * cost / ranks,
        "ranks": ranks, "basis": "explicit Courant on sqrt(U^2+5RT)+c; "
                                 "seconds_per_cfd_cell_step from configuration",
    }


def dsmc_estimate(cfg: CaseConfig, source: SourceState, *, number_density, temperature,
                  speed, volumes, sizes, length_scale_m: float, cell_weight=1.0,
                  basis: str) -> dict:
    """Parcels, memory, step, runtime and resolution for a set of DSMC cells."""
    gas = source.gas
    n = np.asarray(number_density, float)
    t = np.asarray(temperature, float)
    u = np.asarray(speed, float)
    v = np.asarray(volumes, float)
    h = np.asarray(sizes, float)
    w_cells = np.broadcast_to(np.asarray(cell_weight, float), n.shape)
    ok = (n > 0) & (t > 0) & np.isfinite(n) & np.isfinite(t)
    if not np.any(ok):
        raise ValueError("no DSMC cell has a positive density and temperature")
    dx = h.max(axis=1) if h.ndim == 2 else h
    molecules = n * v
    ppc = float(cfg.get("dsmc.parcels_per_cell_target", 20))
    weight = float(np.max(molecules[ok])) / ppc
    parcels = float(np.sum(w_cells[ok] * molecules[ok]) / weight)
    lam = gas.mean_free_path(n, t, cfg.get("breakdown.mean_free_path_convention", "vhs"))
    tau = lam / gas.mean_speed(t)
    collisions = bool(cfg.get("dsmc.collisions", True))
    sigma = np.sqrt(gas.gas_constant * t)
    signal = float(np.max((u + 3.0 * sigma)[ok]))
    dt_courant = float(cfg.get("dsmc.max_courant", 0.2)) * float(np.min(dx[ok])) / signal
    dt_coll = float(cfg.get("dsmc.max_dt_over_collision_time", 0.5)) * float(np.nanmin(tau[ok]))
    dt = min(dt_courant, dt_coll) if collisions else dt_courant
    transit = length_scale_m / source.velocity
    duration = (float(cfg.get("dsmc.startup_transits", 2.0))
                + float(cfg.get("dsmc.sample_transits", 3.0))) * transit
    steps = duration / dt
    ranks = int(cfg.get("resources.max_ranks", 1))
    cost = float(cfg.get("resources.seconds_per_parcel_step", 1e-6))
    mem = (parcels * float(cfg.get("resources.bytes_per_parcel", 300))
           + float(np.sum(w_cells)) * float(cfg.get("resources.bytes_per_cell", 2500))) / GB
    ratio = dx / lam
    occ = molecules / weight
    wsum = float(np.sum(w_cells[ok] * v[ok]))

    def frac(mask):
        return float(np.sum((w_cells * v)[ok & mask]) / wsum)

    densest = int(np.argmax(np.where(ok, n, -np.inf)))
    return {
        "basis": basis,
        "cells": float(np.sum(w_cells)),
        "particle_weight": weight, "parcels": parcels,
        "memory_GB": mem, "memory_budget_GB": float(cfg.get("resources.max_memory_GB", 8)),
        "fits_memory_budget": bool(mem <= float(cfg.get("resources.max_memory_GB", 8))),
        "dt_s": dt, "dt_courant_s": dt_courant,
        "dt_collision_limit_s": dt_coll if collisions else None,
        "dt_over_min_collision_time": dt / float(np.nanmin(tau[ok])) if collisions else None,
        "steps": steps, "duration_s": duration,
        "runtime_estimate_s_at_max_ranks": parcels * steps * cost / ranks,
        "ranks": ranks,
        "cell_over_mfp": {
            "max": float(np.nanmax(ratio[ok])),
            "volume_fraction_above_1": frac(ratio > 1.0),
            "volume_fraction_above_3": frac(ratio > 3.0),
        } if collisions else "collisions disabled: not a criterion",
        "occupancy": {
            "volume_fraction_below_1_parcel": frac(occ < 1.0),
            "volume_fraction_below_5_parcels": frac(occ < 5.0),
            "volume_fraction_below_20_parcels": frac(occ < 20.0),
        },
        "densest_cell": {"number_density": float(n[densest]),
                         "mean_free_path_m": float(lam[densest]),
                         "cell_size_m": float(dx[densest])},
        "transit_s": transit,
    }


def assess_configuration(cfg: CaseConfig, source: SourceState, plan: MeshPlan) -> dict:
    """Pre-CFD: CFD cost, and a full-domain DSMC estimate from the collisionless field."""
    from plumetools.cai2012 import analytical
    centres, volumes, sizes, weight = _lattice_sample(plan)
    n = source.number_density * analytical.density_ratio(
        centres[:, 0], centres[:, 1], centres[:, 2], source.radius_m, source.speed_ratio)
    full = dsmc_estimate(
        cfg, source, number_density=n, temperature=np.full(n.shape, source.temperature_K),
        speed=np.full(n.shape, source.velocity), volumes=volumes, sizes=sizes,
        length_scale_m=plan.x_segments[-1].hi, cell_weight=weight,
        basis="full-domain DSMC, collisionless density at T = Te (shape estimate)")
    cfd = (cfd_estimate(cfg, source, plan) if cfg.get("interface.mode", "automatic") != "none"
           else {"not_run": "interface.mode is none: full DSMC, no CFD precursor"})
    # Cai's startup is 1e4 t0 with t0 the Kne = 0.01 collision time, for EVERY case.
    t0_cai = source.with_knudsen(0.01).collision_time
    startup = float(cfg.get("dsmc.startup_transits", 2.0)) * full["transit_s"]
    return {"cfd": cfd, "full_dsmc": full,
            "source_cell_over_mfp": plan.core_cell_m / source.mean_free_path,
            "cai_startup": {"t0_kn0p01_s": t0_cai, "cai_1e4_t0_s": 1e4 * t0_cai,
                            "configured_startup_s": startup,
                            "note": "reported, not imposed; stationarity is measured"}}


_COUNT = re.compile(r"^\s*(\d+)\s*$")


def lagrangian_count(path: Path) -> int | None:
    """Parcel count from the list header of a lagrangian field file."""
    path = Path(path)
    if not path.is_file():
        return None
    with path.open() as fh:
        for line in fh:
            m = _COUNT.match(line)
            if m:
                return int(m.group(1))
    return None


def rank_imbalance(case_dir: Path, time: str, cloud: str = "dsmc") -> dict | None:
    """Measured parcels per rank at one written time, if the case is decomposed."""
    procs = sorted(Path(case_dir).glob("processor*"),
                   key=lambda p: int(p.name.replace("processor", "") or 0))
    counts = []
    for p in procs:
        c = lagrangian_count(p / time / "lagrangian" / cloud / "U")
        counts.append(0 if c is None else c)
    if not counts:
        return None
    arr = np.asarray(counts, float)
    return {"ranks": len(counts), "parcels_per_rank": counts,
            "max_over_mean": float(arr.max() / arr.mean()) if arr.mean() > 0 else None}


def measured_rates(log_path: Path, parcels: float, cells: int, ranks: int) -> dict | None:
    """Seconds per parcel-step from a dsmcFoam log's ExecutionTime lines."""
    path = Path(log_path)
    if not path.is_file():
        return None
    times = [float(m.group(1)) for m in re.finditer(
        r"ExecutionTime = ([0-9.eE+-]+) s", path.read_text(errors="replace"))]
    if len(times) < 3 or parcels <= 0:
        return None
    steps = len(times) - 1
    per_step = (times[-1] - times[0]) / steps
    return {"steps": steps, "seconds_per_step": per_step,
            "seconds_per_parcel_step_x_ranks": per_step * ranks / parcels,
            "cells": cells, "ranks": ranks}
