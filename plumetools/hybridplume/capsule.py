"""Exact conservative capsule candidate sweep over full 3D CFD voxels.

A capsule is a cylinder of radius R over 0 <= x <= L capped by a spheroidal
tip with centre (L, 0, 0) and semi-axes (a, R, R). For each (R, L) on a grid
the largest admissible tip length a is computed exactly from the excluded
voxels, reduced to axial slabs; every candidate is then re-certified with
`Capsule.verify`. Each candidate is then screened with the half-range
Maxwellian return-mass fraction and the subsonic normal-Mach area fraction of
the CFD state on its analytic surface, and the largest passing one selected.
"""
from __future__ import annotations

import math

import numpy as np

from . import halfrange
from . import workflow as wf
from .shapes import Capsule

SCREEN_DEFAULTS = {"axial_points_per_D": 40, "azimuthal_points": 90}


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


def _positive_int(value, name):
    if isinstance(value, bool) or int(value) != value or value < 1:
        raise ValueError(f"{name} must be a positive integer, got {value!r}")
    return int(value)


def surface_samples(shape, axial_points_per_D, azimuthal_points, D):
    """Midpoint samples of the cylinder wall and spheroidal tip (no base disc).

    Returns ``(points, normals, areas)``: points on the surface, unit outward
    normals and the area each sample represents. The cylinder is sampled at
    midpoints in (x, theta), the tip at midpoints in (polar angle, theta) with
    ``x = L + a sin(phi)``, ``r = R cos(phi)``; ``L = 0`` gives a pure tip.
    """
    per_d = _positive_int(axial_points_per_D, "axial_points_per_D")
    n_theta = _positive_int(azimuthal_points, "azimuthal_points")
    D = float(D)
    if not (math.isfinite(D) and D > 0):
        raise ValueError("source diameter D must be finite and positive")
    R, L, a = float(shape.R), float(shape.L), float(shape.a)
    theta = (np.arange(n_theta)+.5)*(2*math.pi/n_theta)
    d_theta = 2*math.pi/n_theta
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    points, normals, areas = [], [], []
    if L > 0:
        n_x = max(1, math.ceil(per_d*L/D))
        x, _ = np.meshgrid((np.arange(n_x)+.5)*(L/n_x), theta, indexing="ij")
        c, s = np.broadcast_to(cos_t, x.shape), np.broadcast_to(sin_t, x.shape)
        points.append(np.column_stack((x.ravel(), R*c.ravel(), R*s.ravel())))
        normals.append(np.column_stack((np.zeros(x.size), c.ravel(), s.ravel())))
        areas.append(np.full(x.size, R*d_theta*L/n_x))
    n_phi = max(8, math.ceil(per_d*max(a, R)/D))
    d_phi = .5*math.pi/n_phi
    phi, _ = np.meshgrid((np.arange(n_phi)+.5)*d_phi, theta, indexing="ij")
    c, s = np.broadcast_to(cos_t, phi.shape), np.broadcast_to(sin_t, phi.shape)
    x, r = L+a*np.sin(phi), R*np.cos(phi)
    dx, dr = a*np.cos(phi), -R*np.sin(phi)
    ds = np.hypot(dx, dr)
    nx, nr = -dr/ds, dx/ds
    points.append(np.column_stack((x.ravel(), (r*c).ravel(), (r*s).ravel())))
    normals.append(np.column_stack((nx.ravel(), (nr*c).ravel(), (nr*s).ravel())))
    areas.append((r*ds*d_phi*d_theta).ravel())
    return np.vstack(points), np.vstack(normals), np.concatenate(areas)


def screen(candidates, locator, state, retained, gas, D, cfg) -> list[dict]:
    """Half-range flux screen of each candidate on its analytic surface.

    Adds ``return_mass_fraction`` (reverse/forward Maxwellian mass through the
    surface), ``subsonic_area_fraction`` (area with normal Mach < 1) and
    ``screen_feasible``. A candidate whose samples leave the CFD domain, land
    on cells outside ``retained`` or meet an invalid CFD state is marked
    screen-infeasible with the other two fields ``None``; the sweep goes on.
    """
    params = {**SCREEN_DEFAULTS, **(cfg.get("interface.fit.screen", {}) or {})}
    per_d, n_theta = params["axial_points_per_D"], params["azimuthal_points"]
    out = []
    for cand in candidates:
        row = {**cand, "return_mass_fraction": None, "subsonic_area_fraction": None,
               "screen_feasible": False}
        out.append(row)
        points, normals, areas = surface_samples(
            Capsule(float(cand["R"]), float(cand["L"]), float(cand["a"])), per_d, n_theta, D)
        try:
            cells = locator.locate(points, retained)
        except wf.WorkflowError:
            continue
        rho, velocity, temperature = state.rho[cells], state.U[cells], state.T[cells]
        if not ((rho > 0) & (temperature > 0) & np.isfinite(rho) & np.isfinite(temperature)
                & np.isfinite(velocity).all(axis=1)).all():
            continue
        forward = halfrange.integrate(halfrange.incoming_flux(
            rho, velocity, temperature, normals, gas.gas_constant), areas)["mass"]
        reverse = halfrange.integrate(halfrange.incoming_flux(
            rho, velocity, temperature, -normals, gas.gas_constant), areas)["mass"]
        if not (forward > 0 and np.isfinite(reverse)):
            continue
        mach = np.einsum("ij,ij->i", velocity, normals)/gas.sound_speed(temperature)
        row.update(return_mass_fraction=float(reverse/forward),
                   subsonic_area_fraction=float(areas[mach < 1].sum()/areas.sum()),
                   screen_feasible=True)
    return out


def select(screened, crit) -> dict:
    """Largest-volume screen-passing capsule, else the lowest-return feasible one."""
    max_return = float(crit.get("max_return_mass_fraction", 0.05))
    max_subsonic = float(crit.get("max_subsonic_normal_area_fraction", 0.25))
    feasible = [c for c in screened if c.get("screen_feasible")]
    if not feasible:
        reason = ("no certified-feasible capsule candidate" if not screened else
                  "no capsule candidate could be screened on the CFD state "
                  "(samples outside the domain, on excluded cells or on invalid states)")
        return {"status": "no_feasible_capsule", "reason": reason}
    passing = [c for c in feasible if c["return_mass_fraction"] <= max_return
               and c["subsonic_area_fraction"] <= max_subsonic]
    if passing:
        best, status = max(passing, key=lambda c: c["volume_m3"]), "passed"
    else:
        best, status = min(feasible, key=lambda c: c["return_mass_fraction"]), "no_passing_candidate"
    shape = Capsule(float(best["R"]), float(best["L"]), float(best["a"]))
    return {"status": "fitted", "shape": shape.to_json(),
            "screen": {"status": status, "candidates": list(screened)}}
