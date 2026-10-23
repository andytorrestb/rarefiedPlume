r"""Compare a sampled 3-D field with Cai's collisionless solution, off the axis.

Outputs 1 and 2 of the extension -- the centreline ``n/n0``, ``U sqrt(beta0)``,
``T/T0`` and the ``n/n0`` contours -- are :mod:`plumetools.cai2012.post`'s own,
produced by every case's unchanged ``postProcess.py``. This module adds the rest
and nothing it adds is a new physical model: every analytical value comes from
:func:`plumetools.cai2012.analytical.moments`, the free-molecular disk integral
already in the repository, evaluated at new points.

Conventions
-----------
**Coordinates.** ``x`` along the plume axis (+x downstream), the exit disk
centred on the origin in ``x = 0``. ``r = sqrt(y^2 + z^2)`` is the cylindrical
radius and ``R = sqrt(x^2 + r^2)`` the spherical one, both from the centre of
the exit; ``theta = atan2(r, x)`` is the polar angle from the axis.

**Velocity, 3-D to meridional.** The Cartesian sampled velocity is projected,
cell by cell, onto

.. code-block:: text

    u_x                                axial
    u_r     = (y u_y + z u_z) / r      CYLINDRICAL radial: away from the axis
    u_phi   = (y u_z - z u_y) / r      swirl: zero for this flow, reported as a
                                       symmetry diagnostic
    u_R     = cos(theta) u_x + sin(theta) u_r     SPHERICAL radial: away from
    u_theta = -sin(theta) u_x + cos(theta) u_r    the exit centre; polar

``u_r`` and ``u_R`` are different quantities and are never interchanged: the
meridional plane plots ``u_r``, the r/D = 10 profile uses ``u_R`` because that
is the component that carries mass through the sphere. On the axis ``r = 0``
the projections are 0/0; the limit for an axisymmetric field is ``u_r = u_phi
= 0``, which is what is written there. Cartesian cell centres never lie exactly
on the axis here (an even number of core cells spans it), so in practice only
the analytical grid meets the limit.

**Normalisation.** Cai's: ``n/n0``, ``u sqrt(beta0)`` with
``beta0 = 1/(2 R T0)`` (so velocities are local speed ratios), ``T/T0``.

**Pressure** is the scalar kinetic pressure, one third of the trace of the
pressure tensor, ``p = n k T`` with ``T`` the LOCAL translational temperature
from the sampled moments, ``3 R T = <v^2> - |u|^2``. So ``p/p0 = (n/n0)(T/T0)``
-- never ``n k T0``, which would be the exit temperature carried through the
plume. In a collisionless expansion the tensor is anisotropic and this is its
mean; the same definition is applied to the analytical solution.

**Mass flux through the r/D = 10 sphere**, per unit solid angle, normalised by
the mass flow out of the exit:

.. math::

    F(\theta) = \frac{n\,u_R\,R^2}{\Gamma_0\,\pi R_0^2},\qquad
    \Gamma_0 = \frac{n_0}{2\sqrt\pi\sqrt{\beta_0}}
               \left[e^{-S_0^2} + \sqrt\pi S_0 (1 + \mathrm{erf}\,S_0)\right]

with ``Gamma0`` the one-way number flux of the exit Maxwellian (the molecular
mass cancels). Integrating needs the axisymmetric area factor,
``dOmega = 2 pi sin(theta) dtheta``, and in the collisionless solution

.. math::

    \int_0^{\pi/2} F(\theta)\, 2\pi \sin\theta\, d\theta = 1

exactly: every molecule leaves the disk with ``v_x > 0`` from inside the sphere
and crosses it once, at ``theta < 90 deg``. :func:`angular_profile` reports the
integral for the DSMC field, for the analytical solution through the same cells,
and for the analytical solution on the exact sphere.

Matched estimators
------------------
A DSMC value is an average over cells -- a tube on the axis, a shell around the
sphere, two cell layers around the centre plane -- and the analytical solution
has gradients across all of them. As in :func:`plumetools.cai2012.post.centerline`,
the analytical solution is therefore evaluated **at the same cell centres and
averaged with the same weights**, and errors are computed against that. The
analytical value on the exact line, sphere or plane is carried alongside for the
figures. Densities, fluxes and pressures are volume-weighted; velocities and
temperatures are molecule-weighted.

Errors
------
:func:`masked_metrics` reports a relative error only where it is well
conditioned and an absolute one everywhere it is defined:

* a point enters the comparison only where the analytical ``n/n0 >=``
  :data:`DENSITY_FLOOR_OVER_N0` (Cai's lowest contour) -- outside it the plume
  is vacuum and relative errors measure nothing. On the r/D = 10 sphere, which
  lies below that contour beyond ~40 deg, the floor is instead
  :data:`ANGULAR_FLOOR_FRACTION` of the sphere's own peak, and the flux floor
  :data:`FLUX_FLOOR_FRACTION` of the peak flux;
* within :data:`LIP_EXCLUSION_OVER_D` of the nozzle lip -- the circle
  ``x = 0, r = R0``, where the analytical density is discontinuous and a
  cell-centre value is not a cell average -- points are excluded and counted;
* velocity, temperature and pressure also need DSMC molecules in the cell; an
  empty cell is counted as excluded, not scored as zero. An empty cell inside
  the plume IS scored for the density, as a 100% deficit;
* ``u_r`` has no relative error at all: it passes through zero on the axis;
* absolute errors are in the normalised units, divided by a stated scale (1 for
  ``n/n0``, ``T/T0``, ``p/p0``; ``S0`` for velocities; the peak analytical ``F``
  for the mass flux).

Every metric carries its valid and excluded counts.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import erf

from plumetools.cai2012 import analytical
from plumetools.cai2012.audit import PLUME_FLOOR_OVER_N0

#: Radius of the angular comparison sphere, in nozzle diameters.
ANGULAR_RADIUS_OVER_D = 10.0

#: Half-thickness of the shell of cells averaged onto that sphere, in D.
#:
#: The outer cells of the graded mesh are up to ~0.5 D across at r/D = 10, so a
#: thinner shell would have gaps. The matched analytical estimator removes the
#: bias a thick shell would otherwise introduce.
SHELL_HALF_WIDTH_OVER_D = 0.5

#: Width of the polar-angle bins on the sphere [deg].
ANGULAR_BIN_DEG = 2.5

#: The analytical ``n/n0`` below which a point is outside the compared plume.
DENSITY_FLOOR_OVER_N0 = PLUME_FLOOR_OVER_N0

#: Radius of the zone around the nozzle lip excluded from error metrics, in D.
LIP_EXCLUSION_OVER_D = 0.1

#: Fraction of the peak analytical F below which the flux is not compared
#: relatively.
FLUX_FLOOR_FRACTION = 1.0e-3

#: Fraction of the profile's peak analytical ``n/n0`` below which an angular
#: bin is not compared.
#:
#: Relative to the PEAK ON THE SPHERE, not to ``n0``: at r/D = 10 the whole
#: plume is below Cai's lowest contour (``n/n0 = 1e-3``) beyond about 40 deg,
#: so the absolute floor used in the plane would leave the angular comparison
#: almost nothing to compare.
ANGULAR_FLOOR_FRACTION = 1.0e-3


# --------------------------------------------------------------------------- #
# geometry and projections
# --------------------------------------------------------------------------- #

def projections(centres, vectors) -> dict:
    """Cylindrical and spherical components of a Cartesian vector field.

    Args:
        centres: ``(n, 3)`` positions [m] relative to the exit centre.
        vectors: ``(n, 3)`` vectors at those positions, in any units.

    Returns:
        ``{"r", "R", "theta", "ux", "ur", "uphi", "uR", "utheta"}``, each
        ``(n,)``. ``theta`` in radians. On the axis ``ur = uphi = 0``.
    """
    centres = np.asarray(centres, dtype=np.float64).reshape(-1, 3)
    vectors = np.asarray(vectors, dtype=np.float64).reshape(-1, 3)
    x, y, z = centres[:, 0], centres[:, 1], centres[:, 2]
    ux, uy, uz = vectors[:, 0], vectors[:, 1], vectors[:, 2]

    r = np.hypot(y, z)
    big_r = np.hypot(x, r)
    scale = max(float(np.max(big_r, initial=0.0)), 1.0e-300)
    on_axis = r <= 1.0e-12 * scale
    safe_r = np.where(on_axis, 1.0, r)
    ur = np.where(on_axis, 0.0, (y * uy + z * uz) / safe_r)
    uphi = np.where(on_axis, 0.0, (y * uz - z * uy) / safe_r)

    origin = big_r <= 0.0
    safe_big = np.where(origin, 1.0, big_r)
    cos_t = np.where(origin, 1.0, x / safe_big)
    sin_t = np.where(origin, 0.0, r / safe_big)
    return {
        "r": r,
        "R": big_r,
        "theta": np.arctan2(r, x),
        "ux": ux,
        "ur": ur,
        "uphi": uphi,
        "uR": cos_t * ux + sin_t * ur,
        "utheta": -sin_t * ux + cos_t * ur,
    }


def exit_flux_factor(speed_ratio: float) -> float:
    r"""``Gamma0 sqrt(beta0) / n0`` -- the exit one-way flux, normalised.

    ``[e^{-S0^2} + sqrt(pi) S0 (1 + erf S0)] / (2 sqrt(pi))``. 1/(2 sqrt(pi))
    for effusion, and -> ``S0`` for a fast drift.
    """
    s = float(speed_ratio)
    return (math.exp(-s * s) + math.sqrt(math.pi) * s * (1.0 + math.erf(s))) / (
        2.0 * math.sqrt(math.pi))


def mass_flux_per_steradian(n_over_n0, uR_sqrt_beta0, big_r, radius_m: float,
                            speed_ratio: float):
    """``F = n u_R R^2 / (Gamma0 pi R0^2)``, from normalised quantities."""
    return (np.asarray(n_over_n0) * np.asarray(uR_sqrt_beta0)
            * (np.asarray(big_r) / float(radius_m)) ** 2
            / (math.pi * exit_flux_factor(speed_ratio)))


def lip_zone(centres, radius_m: float, diameter_m: float,
             exclusion_over_D: float = LIP_EXCLUSION_OVER_D):
    """Points within ``exclusion_over_D * D`` of the nozzle lip circle."""
    centres = np.asarray(centres, dtype=np.float64).reshape(-1, 3)
    r = np.hypot(centres[:, 1], centres[:, 2])
    distance = np.hypot(centres[:, 0], r - float(radius_m))
    return distance < float(exclusion_over_D) * float(diameter_m)


# --------------------------------------------------------------------------- #
# the two fields, in the same normalised form
# --------------------------------------------------------------------------- #

def dsmc_fields(sampled, exit_state, *, mask=None) -> dict:
    """Normalised DSMC moments at (a subset of) the sampled cells.

    Args:
        sampled: a :class:`~plumetools.cai2012.post.SampledField` -- the
            time-averaged fields.
        exit_state: for ``n0``, ``T0`` and ``beta0``.
        mask: optional boolean selection of cells.

    Returns:
        ``{"centres", "volumes", "n_over_n0", "T_over_T0", "p_over_p0",
        "occupied", <projections of u sqrt(beta0)>}``.
    """
    select = slice(None) if mask is None else np.asarray(mask, dtype=bool)
    centres = sampled.centres[select]
    density = sampled.number_density[select]
    n = density / float(exit_state.number_density_per_m3)
    temperature = sampled.temperature[select] / float(exit_state.T0_K)
    velocity = sampled.velocity[select] * math.sqrt(float(exit_state.beta0))
    out = {
        "centres": centres,
        "volumes": sampled.volumes[select],
        "n_over_n0": n,
        "T_over_T0": temperature,
        "p_over_p0": n * temperature,
        "occupied": density > 0.0,
    }
    out.update(projections(centres, velocity))
    return out


def analytical_fields(centres, radius_m: float, speed_ratio: float, *,
                      n_radial: int = 48, n_azimuthal: int = 96) -> dict:
    """The collisionless solution at arbitrary points, normalised as above.

    :func:`plumetools.cai2012.analytical.moments`, called in bands of ``x``.
    That function raises its quadrature order for the point nearest the exit
    plane and applies it to every point it is given; banding keeps the fine
    rule for the few points that need it. Points at ``x <= 0`` are NaN -- the
    disk integral is singular in the exit plane.
    """
    centres = np.asarray(centres, dtype=np.float64).reshape(-1, 3)
    x, y, z = centres[:, 0], centres[:, 1], centres[:, 2]
    count = x.size
    n = np.full(count, np.nan)
    velocity = np.full((count, 3), np.nan)
    temperature = np.full(count, np.nan)

    radius = float(radius_m)
    # Below this x the moments routine refines (two nodes across the 1/R^2
    # peak); at and above it the base order is already enough.
    threshold = max(2.0 * math.pi * radius / n_radial,
                    4.0 * math.pi * radius / n_azimuthal)
    floor = 4.0 * math.pi * radius / analytical.MAX_AZIMUTHAL
    edges = [math.inf, threshold]
    while edges[-1] / 2.0 > floor:
        edges.append(edges[-1] / 2.0)
    edges.append(0.0)

    for upper, lower in zip(edges[:-1], edges[1:]):
        band = (x > lower) & (x <= upper) if math.isfinite(upper) else (x > lower)
        if not np.any(band):
            continue
        result = analytical.moments(x[band], y[band], z[band], radius,
                                    speed_ratio, n_radial=n_radial,
                                    n_azimuthal=n_azimuthal)
        n[band] = result["n_over_n0"]
        velocity[band] = result["U_sqrt_beta0"]
        temperature[band] = result["T_over_T0"]

    out = {
        "centres": centres,
        "n_over_n0": n,
        "T_over_T0": temperature,
        "p_over_p0": n * temperature,
    }
    out.update(projections(centres, velocity))
    return out


# --------------------------------------------------------------------------- #
# error metrics
# --------------------------------------------------------------------------- #

def masked_metrics(measured, exact, *, name: str, valid, scale: float,
                   relative: bool = True) -> dict:
    """Relative and normalised absolute errors over the valid points only.

    Args:
        measured: DSMC values.
        exact: analytical values at the same points, same estimator.
        name: key prefix.
        valid: boolean mask of points that enter the comparison; the caller
            decides it and records why the others are out.
        scale: what absolute errors are divided by.
        relative: also report ``|m - e| / |e|``. Off for quantities that pass
            through zero.

    Returns:
        ``{name}_n_total``, ``_n_valid``, ``_n_excluded``, ``_abs_error_scale``,
        ``_max/_mean/_rms_abs_error`` and, when ``relative``, the same three
        ``_rel_error`` keys plus ``_max_rel_error_at_index`` -- fractions, not
        percents. NaN where nothing is valid.
    """
    measured = np.asarray(measured, dtype=np.float64).ravel()
    exact = np.asarray(exact, dtype=np.float64).ravel()
    valid = np.asarray(valid, dtype=bool).ravel()
    if not measured.shape == exact.shape == valid.shape:
        raise ValueError(f"{name}: {measured.shape} measured, {exact.shape} "
                         f"exact, {valid.shape} mask")
    use = valid & np.isfinite(measured) & np.isfinite(exact)
    if relative:
        use &= np.abs(exact) > 0.0
    out = {
        f"{name}_n_total": int(measured.size),
        f"{name}_n_valid": int(np.count_nonzero(use)),
        f"{name}_n_excluded": int(measured.size - np.count_nonzero(use)),
        f"{name}_abs_error_scale": float(scale),
    }
    keys = ["abs"] + (["rel"] if relative else [])
    if not np.any(use):
        for kind in keys:
            for stat in ("max", "mean", "rms"):
                out[f"{name}_{stat}_{kind}_error"] = float("nan")
        if relative:
            out[f"{name}_max_rel_error_at_index"] = -1
        return out

    difference = np.abs(measured[use] - exact[use])
    errors = {"abs": difference / float(scale)}
    if relative:
        errors["rel"] = difference / np.abs(exact[use])
    for kind, values in errors.items():
        out[f"{name}_max_{kind}_error"] = float(values.max())
        out[f"{name}_mean_{kind}_error"] = float(values.mean())
        out[f"{name}_rms_{kind}_error"] = float(math.sqrt(float((values ** 2).mean())))
    if relative:
        out[f"{name}_max_rel_error_at_index"] = int(
            np.flatnonzero(use)[int(np.argmax(errors["rel"]))])
    return out


# --------------------------------------------------------------------------- #
# output 1, extended: the centreline pressure
# --------------------------------------------------------------------------- #

def centerline_extended(profile: dict, speed_ratio: float) -> tuple[dict, dict]:
    """Add the pressure to :func:`plumetools.cai2012.post.centerline`'s profile.

    In a bin, ``n`` is volume-weighted and ``T`` molecule-weighted, so their
    product is exactly the volume-weighted mean of ``n T`` -- the same estimator
    as every other volume quantity, with no new averaging. The same identity
    holds for the matched analytical columns.

    Returns:
        ``(profile with p_over_p0 columns, metrics)``. The metrics are the
        four centreline quantities through :func:`masked_metrics`, so they carry
        absolute errors and counts beside the relative ones
        ``post.centerline_metrics`` reports.
    """
    out = dict(profile)
    for suffix in ("", "_analytical", "_analytical_axis"):
        out[f"p_over_p0{suffix}"] = (np.asarray(profile[f"n_over_n0{suffix}"])
                                     * np.asarray(profile[f"T_over_T0{suffix}"]))
    exact_n = np.asarray(profile["n_over_n0_analytical"])
    valid = exact_n >= DENSITY_FLOOR_OVER_N0
    metrics = {}
    for key, label, scale in (("n_over_n0", "density", 1.0),
                              ("U_sqrt_beta0", "velocity", float(speed_ratio)),
                              ("T_over_T0", "temperature", 1.0),
                              ("p_over_p0", "pressure", 1.0)):
        metrics.update(masked_metrics(out[key], out[f"{key}_analytical"],
                                      name=label, valid=valid, scale=scale))
    return out, metrics


# --------------------------------------------------------------------------- #
# outputs 3 and 4: the angular profiles at r/D = 10
# --------------------------------------------------------------------------- #

def footprint_fractions(theta, big_r, volumes, edges_rad):
    """How much of each cell's angular footprint falls in each polar bin.

    A cell of volume ``V`` at radius ``R`` subtends about ``V^(1/3) / R`` of
    polar angle. Out at r/D = 10 the graded mesh's cells are 0.4-0.5 D across --
    as wide as a 2.5 degree bin there -- so assigning each cell to the one bin
    its CENTRE falls in leaves some bins empty and the profile full of holes.
    Each cell is instead shared among the bins its footprint
    ``[theta - h, theta + h]`` overlaps, in proportion to the overlap.

    The same fractions are applied to the analytical values at the same cells,
    so the matched comparison is unaffected by the smoothing this implies.

    Returns:
        ``(n_cells, n_bins)`` fractions; each row sums to 1, or to 0 for a cell
        whose footprint misses every bin.
    """
    theta = np.asarray(theta, dtype=np.float64)
    half = 0.5 * np.cbrt(np.asarray(volumes, dtype=np.float64)) / np.maximum(
        np.asarray(big_r, dtype=np.float64), 1.0e-300)
    lo = (theta - half)[:, None]
    hi = (theta + half)[:, None]
    overlap = np.clip(np.minimum(hi, edges_rad[None, 1:])
                      - np.maximum(lo, edges_rad[None, :-1]), 0.0, None)
    total = overlap.sum(axis=1, keepdims=True)
    return np.divide(overlap, total, out=np.zeros_like(overlap),
                     where=total > 0.0)


def _binned(fractions, weights, values):
    """Per-bin weighted mean of ``values`` over shared cells; NaN if empty."""
    shared = fractions * np.asarray(weights, dtype=np.float64)[:, None]
    total = shared.sum(axis=0)
    weighted = (shared * np.nan_to_num(np.asarray(values, dtype=np.float64))[:, None]
                ).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0.0, weighted / np.where(total > 0.0, total, 1.0),
                        np.nan)


def _solid_angle(edges_rad):
    """``2 pi (cos a - cos b)`` per bin -- the axisymmetric area factor."""
    return 2.0 * math.pi * (np.cos(edges_rad[:-1]) - np.cos(edges_rad[1:]))


def angular_profile(sampled, geom, exit_state, *,
                    radius_over_D: float = ANGULAR_RADIUS_OVER_D,
                    half_width_over_D: float = SHELL_HALF_WIDTH_OVER_D,
                    bin_deg: float = ANGULAR_BIN_DEG,
                    parcel_count=None) -> dict:
    """Density, velocity, temperature, pressure and mass flux on a sphere.

    Args:
        sampled: the time-averaged DSMC field.
        geom: the geometry, for ``D`` and ``R0``.
        exit_state: for the normalisation.
        radius_over_D: the sphere, ``R/D``. 10 by default.
        half_width_over_D: cells whose radial extent reaches
            ``|R - R*| <= this * D`` are averaged, each shared among the polar
            bins its angular footprint overlaps (:func:`footprint_fractions`).
        bin_deg: polar-angle bin width.
        parcel_count: optional ``dsmcRhoNMean`` per cell, to report how many
            parcels stand behind each bin.

    Returns:
        A dict of per-bin arrays -- ``theta_deg`` (bin centres), ``n_cells``,
        and for each of ``n_over_n0``, ``uR``, ``utheta``, ``T_over_T0``,
        ``p_over_p0``, ``F``: the DSMC value, ``*_analytical`` (same cells,
        same weights) and ``*_exact`` (on the sphere at the bin centre) -- plus
        ``integrals``: the solid-angle integral of ``F`` for each of the three.
    """
    d = float(geom.diameter_m)
    radius = float(geom.radius_m)
    s0 = float(exit_state.speed_ratio)
    target = float(radius_over_D) * d
    half = float(half_width_over_D) * d

    x = sampled.centres[:, 0]
    big_r = np.linalg.norm(sampled.centres, axis=1)
    # A cell belongs to the shell if its radial extent reaches it, not only its
    # centre: on the axis the domain ends at x = 10 D, and the last outer cell's
    # centre can sit short of r/D = 9.5 on a coarse mesh.
    reach = half + 0.5 * np.cbrt(np.asarray(sampled.volumes, dtype=np.float64))
    shell = (np.abs(big_r - target) <= reach) & (x > 0.0)

    measured = dsmc_fields(sampled, exit_state, mask=shell)
    exact_cells = analytical_fields(measured["centres"], radius, s0)

    edges_deg = np.arange(0.0, 90.0 + 0.5 * bin_deg, float(bin_deg))
    edges_deg[-1] = min(edges_deg[-1], 90.0)
    fractions = footprint_fractions(measured["theta"], measured["R"],
                                    measured["volumes"], np.radians(edges_deg))

    out = {"theta_deg": 0.5 * (edges_deg[:-1] + edges_deg[1:]),
           "theta_edges_deg": edges_deg,
           # Cells sharing a bin count fractionally: this is the effective
           # number of cells behind each bin.
           "n_cells": fractions.sum(axis=0)}

    volumes = measured["volumes"]
    for label, fields in (("", measured), ("_analytical", exact_cells)):
        flux = mass_flux_per_steradian(fields["n_over_n0"], fields["uR"],
                                       fields["R"], radius, s0)
        molecules = volumes * np.nan_to_num(fields["n_over_n0"])
        for key, weights, values in (
                ("n_over_n0", volumes, fields["n_over_n0"]),
                ("p_over_p0", volumes, fields["p_over_p0"]),
                ("F", volumes, flux),
                ("uR", molecules, fields["uR"]),
                ("utheta", molecules, fields["utheta"]),
                ("T_over_T0", molecules, fields["T_over_T0"]),
                ("uphi", molecules, fields["uphi"])):
            out[f"{key}{label}"] = _binned(fractions, weights, values)

    # On the exact sphere, at the bin centres. phi = 0 is enough: the solution
    # is axisymmetric (asserted in the parent package's tests).
    theta = np.radians(out["theta_deg"])
    sphere = np.column_stack([target * np.cos(theta), np.zeros_like(theta),
                              target * np.sin(theta)])
    exact = analytical_fields(sphere, radius, s0)
    for key in ("n_over_n0", "uR", "utheta", "T_over_T0", "p_over_p0", "uphi"):
        out[f"{key}_exact"] = exact[key]
    out["F_exact"] = mass_flux_per_steradian(exact["n_over_n0"], exact["uR"],
                                             exact["R"], radius, s0)

    omega = _solid_angle(np.radians(edges_deg))
    occupied_bins = out["n_cells"] > 0.0
    out["integrals"] = {
        key: float(np.nansum(np.where(occupied_bins, out[key], 0.0) * omega))
        for key in ("F", "F_analytical")
    }
    out["integrals"]["F_exact"] = float(np.sum(out["F_exact"] * omega))
    # Over the bins that HAVE cells only, so a partly covered sphere is
    # compared like for like. Equal to F_exact when nothing is missing, which
    # is the case on every production mesh of this study.
    out["integrals"]["F_exact_occupied_bins"] = float(
        np.sum(np.where(occupied_bins, out["F_exact"], 0.0) * omega))
    out["integrals"]["bins_without_cells"] = int(np.count_nonzero(~occupied_bins))
    out["solid_angle"] = omega

    if parcel_count is not None:
        parcels = np.asarray(parcel_count, dtype=np.float64)[shell]
        out["parcels"] = (fractions * parcels[:, None]).sum(axis=0)
    out["shell"] = {"radius_over_D": float(radius_over_D),
                    "half_width_over_D": float(half_width_over_D),
                    "bin_deg": float(bin_deg),
                    "n_cells": int(np.count_nonzero(shell))}
    return out


def angular_metrics(profile: dict, speed_ratio: float) -> dict:
    """Errors of the r/D = 10 profiles against the matched analytical ones."""
    exact_n = np.asarray(profile["n_over_n0_analytical"])
    has_cells = np.asarray(profile["n_cells"]) > 0
    peak_n = float(np.nanmax(exact_n)) if np.any(np.isfinite(exact_n)) else 0.0
    in_plume = has_cells & (exact_n >= ANGULAR_FLOOR_FRACTION * peak_n)
    flux = np.asarray(profile["F_analytical"])
    peak = float(np.nanmax(flux)) if np.any(np.isfinite(flux)) else float("nan")
    flux_valid = has_cells & (flux >= FLUX_FLOOR_FRACTION * peak)
    occupied = np.isfinite(np.asarray(profile["uR"]))

    metrics = {}
    metrics.update(masked_metrics(profile["n_over_n0"], exact_n,
                                  name="angular_density", valid=in_plume,
                                  scale=1.0))
    metrics.update(masked_metrics(profile["F"], flux, name="angular_mass_flux",
                                  valid=flux_valid, scale=peak))
    metrics.update(masked_metrics(profile["uR"], profile["uR_analytical"],
                                  name="angular_radial_velocity",
                                  valid=in_plume & occupied,
                                  scale=float(speed_ratio)))
    metrics.update(masked_metrics(profile["T_over_T0"],
                                  profile["T_over_T0_analytical"],
                                  name="angular_temperature",
                                  valid=in_plume & occupied, scale=1.0))
    integrals = profile["integrals"]
    metrics.update({
        "angular_mass_flux_integral": integrals["F"],
        "angular_mass_flux_integral_analytical": integrals["F_analytical"],
        "angular_mass_flux_integral_exact": integrals["F_exact"],
        "angular_mass_flux_integral_exact_occupied_bins":
            integrals["F_exact_occupied_bins"],
        "angular_bins_without_cells": integrals["bins_without_cells"],
    })
    return metrics


# --------------------------------------------------------------------------- #
# outputs 5 and 6: the meridional plane
# --------------------------------------------------------------------------- #

def meridional_cells(sampled, cfg, geom):
    """The cells :func:`plumetools.cai2012.post.density_plane` takes.

    The same rule -- the cell layers within 1.5 of the finest normal offset of
    the centre plane, inside the plotted window -- so the extension's off-axis
    fields and the parent's density contour are drawn from the same cells.
    """
    d = float(geom.diameter_m)
    x = sampled.centres[:, 0]
    if cfg.post.plane == "xz":
        normal, lateral = sampled.centres[:, 1], sampled.centres[:, 2]
    else:
        normal, lateral = sampled.centres[:, 2], sampled.centres[:, 1]
    finest = float(np.min(np.abs(normal))) if normal.size else 0.0
    tolerance = max(finest * 1.5, 1.0e-12)
    half = float(cfg.post.plane_half_over_D) * d
    x_max = float(cfg.post.centerline_x_over_D_max) * d
    return ((np.abs(normal) <= tolerance) & (x >= 0.0) & (x <= x_max)
            & (np.abs(lateral) <= half))


def meridional_plane(sampled, cfg, geom, exit_state) -> dict:
    """Every compared quantity on the centre-plane cells, DSMC and analytical.

    The analytical values are at the SAME cell centres -- the full 3-D position,
    including the half-cell offset from the plane -- so a cell is compared with
    the solution where it actually is.
    """
    mask = meridional_cells(sampled, cfg, geom)
    measured = dsmc_fields(sampled, exit_state, mask=mask)
    exact = analytical_fields(measured["centres"], geom.radius_m,
                              exit_state.speed_ratio)
    lateral_index = 2 if cfg.post.plane == "xz" else 1
    d = float(geom.diameter_m)
    return {
        "x_over_D": measured["centres"][:, 0] / d,
        "lateral_over_D": measured["centres"][:, lateral_index] / d,
        "lateral_axis": "Z" if cfg.post.plane == "xz" else "Y",
        "lip": lip_zone(measured["centres"], geom.radius_m, d),
        "measured": measured,
        "exact": exact,
    }


def plane_metrics(plane: dict, speed_ratio: float) -> dict:
    """Errors over the centre plane, with every exclusion counted."""
    measured, exact = plane["measured"], plane["exact"]
    exact_n = exact["n_over_n0"]
    finite = np.isfinite(exact_n)
    above = finite & (exact_n >= DENSITY_FLOOR_OVER_N0)
    base = above & ~plane["lip"]
    occupied = measured["occupied"]

    metrics = {
        "plane_cells": int(exact_n.size),
        "plane_excluded_below_density_floor": int(np.count_nonzero(finite & ~above)),
        "plane_excluded_nozzle_lip": int(np.count_nonzero(above & plane["lip"])),
        "plane_excluded_no_dsmc_sample": int(np.count_nonzero(base & ~occupied)),
        "plane_swirl_max_abs": float(np.max(np.abs(measured["uphi"][base & occupied]),
                                            initial=0.0)),
    }
    metrics.update(masked_metrics(measured["n_over_n0"], exact_n,
                                  name="plane_density", valid=base, scale=1.0))
    s0 = float(speed_ratio)
    for key, label, scale, relative in (
            ("ux", "plane_axial_velocity", s0, True),
            ("ur", "plane_radial_velocity", s0, False),
            ("T_over_T0", "plane_temperature", 1.0, True),
            ("p_over_p0", "plane_pressure", 1.0, True)):
        metrics.update(masked_metrics(measured[key], exact[key], name=label,
                                      valid=base & occupied, scale=scale,
                                      relative=relative))
    return metrics


def analytical_meridional_grid(cfg, geom, exit_state) -> dict:
    """The analytical fields on a regular grid in the centre plane, for figures.

    The same grid :func:`plumetools.cai2012.post.analytical_plane` uses, with
    every moment rather than the density alone; the non-negative half is
    evaluated and mirrored.
    """
    d = float(geom.diameter_m)
    x_max = float(cfg.post.centerline_x_over_D_max) * d
    half = float(cfg.post.plane_half_over_D) * d
    nx = int(cfg.post.plane_points_x)
    n_half = max(2, (int(cfg.post.plane_points_lateral) + 1) // 2)
    x = np.linspace(0.0, x_max, nx + 1)[1:] - 0.5 * x_max / nx
    half_lateral = np.linspace(0.0, half, n_half)
    xh, lh = np.meshgrid(x, half_lateral, indexing="ij")
    points = np.column_stack([xh.ravel(), np.zeros(xh.size), lh.ravel()])
    if cfg.post.plane != "xz":
        points = points[:, [0, 2, 1]]
    fields = analytical_fields(points, geom.radius_m, exit_state.speed_ratio)

    lateral = np.concatenate([-half_lateral[:0:-1], half_lateral])
    xx, ll = np.meshgrid(x, lateral, indexing="ij")
    grid = {"x_over_D": xx / d, "lateral_over_D": ll / d}
    for key in ("n_over_n0", "ux", "ur", "T_over_T0", "p_over_p0"):
        values = fields[key].reshape(xh.shape)
        grid[key] = np.concatenate([values[:, :0:-1], values], axis=1)
    return grid
