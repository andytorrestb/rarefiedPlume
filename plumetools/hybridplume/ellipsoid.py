"""Conservative spheroid fitting from full three-dimensional CFD voxels."""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import minimize


def clipped_volume(parameters):
    xc, a, b = parameters
    s = np.clip(xc/a, 0, 1)
    return float(math.pi*a*b*b*(2/3+s-s**3/3))


def verify_ellipsoid(fit, lower, upper, retained):
    """Certify the entire ellipsoid against closed excluded voxel footprints."""
    centre = np.array([fit["xc"], 0, 0])
    axes = np.array([fit["a"], fit["b"], fit["b"]])
    margin = float(fit.get("clearance_m", 0))
    bad, checked = 0, 0
    minimum = float("inf")
    for start in range(0, len(lower), 65536):
        excluded = ~np.asarray(retained[start:start+65536], bool)
        lo = np.asarray(lower[start:start+65536])[excluded]-margin
        hi = np.asarray(upper[start:start+65536])[excluded]+margin
        delta = np.maximum(np.maximum(lo-centre, centre-hi), 0)
        q = np.sum((delta/axes)**2, axis=1)
        bad += int(np.count_nonzero(q < 1-1e-8))
        checked += len(q)
        if len(q):
            minimum = min(minimum, float(q.min()))
    return {"feasible": bad == 0, "intersected_excluded_voxels": bad,
            "excluded_voxels_checked": checked,
            "min_scaled_distance_squared": minimum if checked else None}


def fit_ellipsoid(lower, upper, retained, *, source_radius, clearance):
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    retained = np.asarray(retained, bool)
    if (lower.ndim != 2 or lower.shape[1] != 3 or upper.shape != lower.shape
            or retained.shape != (len(lower),) or not len(lower)
            or not np.isfinite(lower).all() or not np.isfinite(upper).all()
            or np.any(upper <= lower) or not np.isfinite(source_radius)
            or source_radius <= 0 or not np.isfinite(clearance) or clearance < 0):
        raise ValueError("fit requires finite 3D voxels, positive source and nonnegative clearance")
    empty = {"status": "no_feasible_ellipsoid", "reason": "no source-connected conservative ellipsoid"}
    if not retained.any():
        return empty
    scale = 2*source_radius
    lo, hi, pad = lower/scale, upper/scale, clearance/scale
    domain_pad = max(pad, .25*float(np.min(hi-lo)))
    xmax = float(hi[:, 0].max()-domain_pad)
    rmax = float(min(hi[:, 1:].max(axis=0).min(), -lo[:, 1:].min(axis=0).max())-domain_pad)
    bad_lo, bad_hi = lo[~retained]-pad, hi[~retained]+pad
    # For equal transverse axes, only the closest rejected radial footprint in
    # each axial slab can be active. No azimuthal mean is taken.
    slabs, inverse = np.unique(np.round(np.c_[bad_lo[:, 0], bad_hi[:, 0]], 12),
                               axis=0, return_inverse=True)
    radial = np.maximum(np.maximum(bad_lo[:, 1:], -bad_hi[:, 1:]), 0)
    radii2 = np.full(len(slabs), np.inf)
    np.minimum.at(radii2, inverse, np.sum(radial*radial, axis=1))
    origin_dx = np.maximum(np.maximum(slabs[:, 0], -slabs[:, 1]), 0)
    distance = float(np.sqrt(origin_dx**2+radii2).min()) if len(slabs) else min(xmax, rmax)
    seed = .4*min(distance, xmax, rmax, .5)
    if seed < 1e-7:
        return empty

    def constraints(p):
        xc, a, b = p
        dx = np.maximum(np.maximum(slabs[:, 0]-xc, xc-slabs[:, 1]), 0)
        source_half = b*np.sqrt(max(0, 1-(xc/a)**2))
        return np.r_[dx**2/a**2+radii2/b**2-1-1e-8,
                     a-xc-1e-7, .5-source_half, xmax-xc-a, rmax-b]

    result = minimize(lambda p: -np.log(clipped_volume(p)), [0, seed, seed],
                      method="SLSQP", bounds=[(0, xmax), (1e-7, xmax), (1e-7, rmax)],
                      constraints={"type": "ineq", "fun": constraints},
                      options={"maxiter": 250, "ftol": 1e-10})
    xc, a, b = result.x
    if (not np.isfinite(result.x).all() or not 0 <= xc < a or b <= 0
            or constraints(result.x).min() < -1e-8):
        return dict(empty, optimizer_message=str(result.message))
    fit = {"status": "fitted", "xc": float(xc*scale), "a": float(a*scale), "b": float(b*scale),
           "volume_m3": clipped_volume(result.x)*scale**3,
           "source_radius_m": float(b*scale*np.sqrt(1-(xc/a)**2)),
           "clearance_m": float(clearance), "domain_clearance_m": domain_pad*scale,
           "optimizer_success": bool(result.success), "optimizer_message": str(result.message),
           "optimality": "local constrained volume optimum; no global optimum claim"}
    verification = verify_ellipsoid(fit, lower, upper, retained)
    if not verification["feasible"]:
        return dict(empty, verification=verification)
    fit["verification"] = verification
    return fit
