"""Exact conservative capsule candidate sweep over full 3D CFD voxels.

A capsule is a cylinder of radius R over 0 <= x <= L capped by a spheroidal
tip with centre (L, 0, 0) and semi-axes (a, R, R). For each (R, L) on a grid
the largest admissible tip length a is computed exactly from the excluded
voxels, reduced to axial slabs; every candidate is then re-certified with
`Capsule.verify`.
"""
from __future__ import annotations

import numpy as np

from .shapes import Capsule


def _excluded_slabs(lower, upper, retained, clearance):
    """Axial slabs of padded excluded voxels with their closest radial distance."""
    bad_lo, bad_hi = lower[~retained]-clearance, upper[~retained]+clearance
    slabs, inverse = np.unique(np.round(np.c_[bad_lo[:, 0], bad_hi[:, 0]], 12),
                               axis=0, return_inverse=True)
    radial = np.maximum(np.maximum(bad_lo[:, 1:], -bad_hi[:, 1:]), 0)
    radii2 = np.full(len(slabs), np.inf)
    np.minimum.at(radii2, np.ravel(inverse), np.sum(radial*radial, axis=1))
    # Capsule.verify ignores footprints entirely upstream of the exit plane.
    keep = slabs[:, 1] >= 0
    return slabs[keep, 0], slabs[keep, 1], np.sqrt(radii2[keep])


def feasible_candidates(lower, upper, retained, *, source_radius: float, clearance: float,
                        radius_samples: int, length_samples: int) -> list[dict]:
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    retained = np.asarray(retained, bool)
    if (lower.ndim != 2 or lower.shape[1] != 3 or upper.shape != lower.shape
            or retained.shape != (len(lower),) or not len(lower)
            or not np.isfinite(lower).all() or not np.isfinite(upper).all()
            or np.any(upper <= lower) or not np.isfinite(source_radius)
            or source_radius <= 0 or not np.isfinite(clearance) or clearance < 0):
        raise ValueError("fit requires finite 3D voxels, positive source and nonnegative clearance")
    if int(radius_samples) != radius_samples or int(length_samples) != length_samples \
            or radius_samples < 1 or length_samples < 1:
        raise ValueError("radius_samples and length_samples must be positive integers")
    if not retained.any():
        return []
    domain_pad = max(clearance, .25*float(np.min(upper-lower)))
    xmax = float(upper[:, 0].max()-domain_pad)
    rmax = float(min(upper[:, 1:].max(axis=0).min(), -lower[:, 1:].min(axis=0).max())-domain_pad)
    r_top = min(source_radius, rmax)-clearance
    if r_top <= 0 or xmax <= 0:
        return []
    x_lo, x_hi, rho = _excluded_slabs(lower, upper, retained, clearance)

    found = []
    for k in range(int(radius_samples)):
        R = r_top*(1-k/radius_samples*.5)
        blocking = rho < R
        l_max = float(np.clip(x_lo[blocking].min(), 0, xmax)) if blocking.any() else xmax
        for L in np.unique(l_max*np.arange(int(length_samples)+1)/length_samples):
            L = float(L)
            a = xmax-L
            ahead = blocking & (x_hi > L)
            if ahead.any():
                reach = np.maximum(x_lo[ahead]-L, 0)/np.sqrt(1-(rho[ahead]/R)**2)
                a = min(a, float(reach.min()))
            if not a > 0:
                continue  # a flat-ended cylinder is not a capsule
            shape = Capsule(float(R), L, a)
            verification = shape.verify(lower, upper, retained, clearance)
            if verification["feasible"]:
                found.append({"R": shape.R, "L": shape.L, "a": shape.a,
                              "volume_m3": shape.volume(), "verification": verification})
    return found
