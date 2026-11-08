"""Measured cell acceptance, reusing shared kinetic-theory functions."""
from __future__ import annotations

import math
import numpy as np

from plumetools.rarefaction.core import mean_free_path
from .model import Assessment, Geometry, Measurements, Targets

EDGE_COMPARISON_RTOL = 1e-12


def usable_window(window: dict) -> bool:
    """Require an identified window of mean per-cell parcel counts."""
    try:
        start, end = float(window["start_s"]), float(window["end_s"])
        return (math.isfinite(start) and math.isfinite(end) and 0 <= start < end
                and window["basis"] == "mean_per_cell"
                and window["source"] in ("declared", "verified")
                and window.get("stationary") is not False)
    except (KeyError, TypeError, ValueError):
        return False


def assess(geometry: Geometry, measured: Measurements, species, *,
           convention: str = "vhs", targets: Targets = Targets()) -> Assessment:
    """Check measured MFP and occupancy without asking passing cells to refine.

    Edge comparisons allow 1e-12 relative floating-point roundoff. Occupancy is
    inclusive (>=20 by default), independently of the generic diagnostic >20
    confidence mask. Invalid or insufficient coverage cannot produce a pass.
    """
    n_cells = geometry.volumes_m3.size
    if measured.number_density.shape != (n_cells,):
        raise ValueError("measurement cell count does not match the geometry")
    n, t, p = measured.number_density, measured.temperature_K, measured.mean_parcels
    lam = mean_free_path(n, t, species, convention=convention)
    desired = targets.edge_over_mfp * lam
    physical = (np.isfinite(n) & (n > 0) & np.isfinite(t) & (t > 0)
                & np.isfinite(lam) & (lam > 0) & np.isfinite(desired) & (desired > 0))
    positive_counts = np.isfinite(p) & (p > 0)
    window_ok = usable_window(measured.window)
    enough = positive_counts & (p >= targets.min_mean_parcels)
    valid = physical & enough & window_ok
    mesh_status = np.full(n_cells, -1, dtype=np.int8)
    mesh_status[valid] = (geometry.max_edge_m[valid]
                          <= desired[valid] * (1 + EDGE_COMPARISON_RTOL)).astype(np.int8)
    occupancy_status = np.full(n_cells, -1, dtype=np.int8)
    count_known = np.isfinite(p) & (p >= 0) & window_ok
    occupancy_status[count_known] = (p[count_known] >= targets.min_mean_parcels).astype(np.int8)
    if not window_ok or np.any(~physical | ~positive_counts):
        action = "unresolved"
    elif np.any(~enough):
        action = "needs_parcels"
    elif np.any(mesh_status == 0):
        action = "needs_mesh"
    else:
        action = "criteria_satisfied"
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ratio = geometry.max_edge_m / lam
    volume = geometry.volumes_m3
    total = float(volume.sum())
    finite_counts = p[np.isfinite(p)]
    summary = {
        "action": action, "cells": n_cells, "window_usable": window_ok,
        "mesh_sha256": geometry.mesh_sha256, "source_sha256": measured.source_sha256,
        "mean_free_path_convention": convention,
        "mesh_failed_cells": int(np.count_nonzero(mesh_status == 0)),
        "occupancy_failed_cells": int(np.count_nonzero(occupancy_status == 0)),
        "unresolved_cells": int(np.count_nonzero(~valid)),
        "mesh_failed_volume_fraction": float(volume[mesh_status == 0].sum() / total),
        "unresolved_volume_fraction": float(volume[~valid].sum() / total),
        "max_edge_over_mfp": float(np.max(ratio[valid])) if np.any(valid) else None,
        "min_mean_parcels": float(np.min(finite_counts)) if finite_counts.size else None,
        "edge_over_mfp_target": targets.edge_over_mfp,
        "mean_parcels_target": targets.min_mean_parcels,
        "edge_comparison_rtol": EDGE_COMPARISON_RTOL,
        "requires_further_mesh_refinement": action == "needs_mesh",
    }
    return Assessment(targets, lam, desired, mesh_status, occupancy_status, valid, summary)
