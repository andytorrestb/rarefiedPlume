"""Conservative axis-aligned ellipse fit to accepted rectangular CFD footprints.

All returned lengths are metres. Constraints inspect excluded footprints, including
interior holes; there is no centre-only classification or threshold adjustment.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import minimize


def footprint_clearance(parameters, lower, upper):
    """Minimum squared ellipse-coordinate distance to each closed rectangle."""
    xc, a, b = parameters
    centre = np.array([xc, 0.0])
    delta = np.maximum(np.maximum(lower - centre, centre - upper), 0.0)
    return np.sum((delta / np.array([a, b])) ** 2, axis=1)


def clipped_area(parameters):
    xc, a, b = parameters
    s = np.clip(xc / a, 0.0, 1.0)
    return float(a * b * (math.pi / 2 + np.arcsin(s) + s * np.sqrt(1 - s*s)))


def ellipse_points(fit, count=101):
    theta = np.linspace(np.arccos(-fit["xc"] / fit["a"]),
                        -np.arccos(-fit["xc"] / fit["a"]), count)
    points = np.column_stack((fit["xc"] + fit["a"] * np.cos(theta),
                              fit["b"] * np.sin(theta)))
    points[[0, -1], 0] = 0.0
    return points


def fit_ellipse(lower, upper, retained, *, source_radius, clearance):
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    retained = np.asarray(retained, bool)
    if (lower.ndim != 2 or lower.shape[1] != 2 or upper.shape != lower.shape
            or retained.shape != (len(lower),) or not len(lower)
            or not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper))
            or np.any(upper <= lower) or not np.isfinite(source_radius)
            or source_radius <= 0 or not np.isfinite(clearance) or clearance < 0):
        raise ValueError("fit requires finite rectangular footprints, a positive source and nonnegative clearance")
    empty = {"status": "no_feasible_ellipse",
             "reason": "no symmetric source-connected ellipse fits the accepted cell footprints"}
    if not retained.any():
        return empty
    scale = 2 * source_radius
    lo, hi = lower / scale, upper / scale
    pad = clearance / scale
    bad_lo, bad_hi = lo[~retained] - pad, hi[~retained] + pad
    xmax = hi[:, 0].max() - pad
    ymax = min(hi[:, 1].max(), -lo[:, 1].min()) - pad
    if min(xmax, ymax) <= 0:
        return empty
    if len(bad_lo):
        distance = np.linalg.norm(np.maximum(np.maximum(bad_lo, -bad_hi), 0), axis=1).min()
    else:
        distance = min(xmax, ymax)
    seed_radius = min(distance, xmax, ymax, 0.5) * 0.4
    if seed_radius < 1e-7:
        return empty
    seed = np.array([0.0, seed_radius, seed_radius])

    def constraints(p):
        xc, a, b = p
        source_half = b * np.sqrt(max(0.0, 1 - (xc/a)**2))
        return np.r_[footprint_clearance(p, bad_lo, bad_hi) - 1.0 - 1e-8,
                     a - xc - 1e-7, 0.5 - source_half,
                     xmax - xc - a, ymax - b]

    result = minimize(lambda p: -np.log(clipped_area(p)), seed, method="SLSQP",
                      bounds=[(0, xmax), (1e-7, xmax), (1e-7, ymax)],
                      constraints={"type": "ineq", "fun": constraints},
                      options={"ftol": 1e-10, "maxiter": 250})
    candidate = np.asarray(result.x)
    # Independent verification of the actual rectangle intersections and bounds.
    xc, a, b = candidate
    nearest_x = np.clip(xc, bad_lo[:, 0], bad_hi[:, 0])
    nearest_y = np.clip(0.0, bad_lo[:, 1], bad_hi[:, 1])
    minimum = ((nearest_x-xc)/a)**2 + (nearest_y/b)**2
    feasible = (np.all(np.isfinite(candidate)) and 0 <= xc < a and b > 0
                and np.all(minimum >= 1 - 1e-9)
                and xc + a <= xmax + 1e-9 and b <= ymax + 1e-9
                and b*np.sqrt(1-(xc/a)**2) <= 0.5 + 1e-9)
    if not feasible:
        return dict(empty, optimizer_message=str(result.message))
    physical = candidate * scale
    return {"status": "fitted", "xc": float(physical[0]), "a": float(physical[1]),
            "b": float(physical[2]), "area_m2": clipped_area(physical),
            "source_half_height_m": float(physical[2]*np.sqrt(1-(xc/a)**2)),
            "clearance_m": float(clearance), "feasible": True,
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message), "optimizer_iterations": int(result.nit),
            "optimality": "local constrained optimum; no global optimality claim",
            "conservative_basis": "accepted CFD cell footprints including interior exclusions"}
