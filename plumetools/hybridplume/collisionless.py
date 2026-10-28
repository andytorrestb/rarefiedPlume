r"""An independent collisionless reference for the circular drifting-Maxwellian source.

The configuration is Cai & Wang (2012), DOI 10.2514/1.A32046. This module does
**not** reuse :mod:`plumetools.cai2012.analytical`: it integrates the raw kinetic
moments numerically, in source area *and* in molecular speed, so agreement
between the two is a check rather than a tautology. It also returns the whole
velocity covariance tensor, which the closed forms do not.

Formulation
-----------
At the source, the number-density distribution in velocity space is

.. math::

    f_e(\mathbf v) = \frac{n_e}{\pi^{3/2} a_e^3}
        \exp\!\left[-|\mathbf v - U_e\hat x|^2/a_e^2\right].

For a field point :math:`\mathbf r` with :math:`x > 0` and a source point
:math:`\mathbf r_s = (0, y_s, z_s)` on the disk, :math:`\ell = |\mathbf r-\mathbf r_s|`,
:math:`\hat q = (\mathbf r - \mathbf r_s)/\ell`, and the solid angle of a source
element is :math:`d\Omega = x\,dA_s/\ell^3`. With no collisions only molecules
travelling along :math:`\hat q` reach :math:`\mathbf r` from that element, so

.. math::

    n = \int_{disk}\frac{x}{\ell^3}\int_0^\infty s^2 f_e(s\hat q)\,ds\,dA_s,\quad
    \mathbf M_1 = \int \frac{x}{\ell^3}\hat q\int_0^\infty s^3 f_e\,ds\,dA_s,\quad
    \mathsf M_2 = \int \frac{x}{\ell^3}\hat q\hat q\int_0^\infty s^4 f_e\,ds\,dA_s,

and :math:`U = M_1/n`, :math:`C = M_2/n - UU`, :math:`T = \mathrm{tr}\,C/(3R)`,
:math:`P = m n C`. The plume is not isotropic: directional temperatures
:math:`C_{ii}/R` are returned alongside the trace temperature.

Quadrature
----------
* source area: Gauss--Legendre in radius, the periodic trapezoid rule in
  azimuth (spectrally accurate for a periodic integrand), with the azimuthal
  nodes offset by half a step so that no node sits on a coordinate axis;
* speed: Gauss--Legendre on ``[0, (S0 + 10) ae]`` -- the integrand beyond is
  below ``exp(-100)`` of its peak.

Every evaluation can be repeated at doubled order (:func:`converged`), and the
relative change is returned as the quadrature-error estimate. Rotational
symmetry is **not** imposed anywhere: the azimuthal nodes break it at the
quadrature-error level, which is what the symmetry test measures.

The source lip
--------------
The integrand is singular where the field point approaches the disk boundary
``x -> 0, r -> Re`` (the disk edge is a discontinuity in the source). Points
closer to the lip circle than ``lip_exclusion`` are refused rather than
evaluated as if smooth; see :func:`lip_distance`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import roots_legendre

from .physics import SourceState

DEFAULT_ORDER = (48, 96, 64)   # (radial, azimuthal, speed)


def lip_distance(points, radius_m: float) -> np.ndarray:
    """Distance from each point to the lip circle ``x = 0, r = Re`` [m]."""
    p = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    r = np.hypot(p[:, 1], p[:, 2])
    return np.hypot(p[:, 0], r - radius_m)


@dataclass(frozen=True)
class Quadrature:
    nodes_y: np.ndarray
    nodes_z: np.ndarray
    weights: np.ndarray       # area weights, sum = pi Re^2
    w: np.ndarray             # speed nodes, in units of ae
    w_weights: np.ndarray


def quadrature(radius_m: float, speed_ratio: float, order=DEFAULT_ORDER) -> Quadrature:
    n_r, n_phi, n_w = (int(v) for v in order)
    x, wx = roots_legendre(n_r)
    rho = 0.5 * radius_m * (x + 1.0)
    w_rho = 0.5 * radius_m * wx * rho                    # includes the Jacobian rho
    phi = (np.arange(n_phi) + 0.5) * (2.0 * math.pi / n_phi)
    w_phi = np.full(n_phi, 2.0 * math.pi / n_phi)
    ny = (rho[:, None] * np.cos(phi)[None, :]).ravel()
    nz = (rho[:, None] * np.sin(phi)[None, :]).ravel()
    weights = (w_rho[:, None] * w_phi[None, :]).ravel()
    w_max = float(speed_ratio) + 10.0
    t, wt = roots_legendre(n_w)
    return Quadrature(ny, nz, weights, 0.5 * w_max * (t + 1.0), 0.5 * w_max * wt)


def raw_moments(points, source: SourceState, order=DEFAULT_ORDER, *,
                lip_exclusion: float = 0.0) -> dict:
    """Number-density moments ``n``, ``M1``, ``M2`` at each point.

    Returns ``{"n": (p,), "M1": (p, 3), "M2": (p, 3, 3)}`` in SI units
    (``M1`` [1/(m^2 s)], ``M2`` [1/(m s^2)]).

    Raises:
        ValueError: for a point with ``x <= 0`` or within ``lip_exclusion`` of
            the lip.
    """
    p = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    if np.any(p[:, 0] <= 0.0):
        raise ValueError("the collisionless reference is defined for x > 0 only")
    if lip_exclusion > 0.0 and np.any(lip_distance(p, source.radius_m) < lip_exclusion):
        raise ValueError(f"point(s) within {lip_exclusion} m of the source lip; "
                         "the integrand is singular there")
    q = quadrature(source.radius_m, source.speed_ratio, order)
    s0 = float(source.speed_ratio)
    ae = source.thermal_speed
    ne = source.number_density
    pref = ne / math.pi ** 1.5

    n = np.empty(p.shape[0])
    m1 = np.empty((p.shape[0], 3))
    m2 = np.empty((p.shape[0], 3, 3))
    for i, (x, y, z) in enumerate(p):
        d = np.stack([np.full_like(q.nodes_y, x), y - q.nodes_y, z - q.nodes_z], axis=1)
        ell = np.linalg.norm(d, axis=1)
        qhat = d / ell[:, None]
        c = qhat[:, 0]
        kernel = q.weights * x / ell ** 3
        # Speed integrals I_k(c) = int w^k exp(-(w - S0 c)^2 - S0^2 (1 - c^2)) dw.
        expo = np.exp(-(q.w[None, :] - s0 * c[:, None]) ** 2
                      - (s0 * s0) * (1.0 - c[:, None] ** 2))
        i2 = (expo * q.w ** 2 * q.w_weights).sum(axis=1)
        i3 = (expo * q.w ** 3 * q.w_weights).sum(axis=1)
        i4 = (expo * q.w ** 4 * q.w_weights).sum(axis=1)
        n[i] = pref * np.dot(kernel, i2)
        m1[i] = pref * ae * (kernel * i3) @ qhat
        m2[i] = pref * ae * ae * np.einsum("k,ki,kj->ij", kernel * i4, qhat, qhat)
    return {"n": n, "M1": m1, "M2": m2}


def macroscopic(raw: dict, source: SourceState) -> dict:
    """Density, velocity, covariance, temperatures and pressure tensor."""
    gas = source.gas
    n = raw["n"]
    u = raw["M1"] / n[:, None]
    cov = raw["M2"] / n[:, None, None] - np.einsum("pi,pj->pij", u, u)
    r = gas.gas_constant
    return {
        "n": n,
        "rho": gas.mass_kg * n,
        "U": u,
        "C": cov,
        "T": np.trace(cov, axis1=1, axis2=2) / (3.0 * r),
        "T_directional": np.diagonal(cov, axis1=1, axis2=2) / r,
        "P": gas.mass_kg * n[:, None, None] * cov,
    }


def reference(points, source: SourceState, order=DEFAULT_ORDER, *,
              lip_exclusion: float = 0.0) -> dict:
    """:func:`raw_moments` then :func:`macroscopic`, plus Cai's normalisations."""
    out = macroscopic(raw_moments(points, source, order, lip_exclusion=lip_exclusion),
                      source)
    out["n_over_ne"] = out["n"] / source.number_density
    out["U_over_ae"] = out["U"] / source.thermal_speed
    out["T_over_Te"] = out["T"] / source.temperature_K
    return out


def converged(points, source: SourceState, order=DEFAULT_ORDER, *,
              lip_exclusion: float = 0.0) -> dict:
    """Evaluate at ``order`` and at twice it; report the larger-order values.

    ``quadrature_rel_change`` holds, per point, the largest relative change of
    ``n``, ``|U|`` and ``T`` between the two orders -- the error estimate of
    the lower order, and a conservative bound for the returned one.
    """
    lo = reference(points, source, order, lip_exclusion=lip_exclusion)
    hi = reference(points, source, tuple(2 * int(v) for v in order),
                   lip_exclusion=lip_exclusion)
    change = np.max(np.stack([
        np.abs(hi["n"] - lo["n"]) / np.abs(hi["n"]),
        np.linalg.norm(hi["U"] - lo["U"], axis=1) / np.linalg.norm(hi["U"], axis=1),
        np.abs(hi["T"] - lo["T"]) / np.abs(hi["T"]),
    ]), axis=0)
    hi["quadrature_rel_change"] = change
    hi["order"] = tuple(2 * int(v) for v in order)
    return hi


def centreline_and_angular_points(source: SourceState, *, x_over_d=None,
                                  radius_over_d: float = 10.0, n_angles: int = 19,
                                  max_angle_deg: float = 90.0,
                                  azimuth_deg: float = 0.0) -> dict:
    """The Cai stations: centreline points and an arc at spherical radius ``r/D``.

    The arc is centred on the source centre and lies in the plane at azimuth
    ``azimuth_deg`` about the x axis; angle 0 is the axis. Points with
    ``x <= 0`` (angle 90 deg) are dropped -- the reference is not defined there.
    """
    d = source.diameter_m
    if x_over_d is None:
        x_over_d = np.array([0.05, 0.1, 0.25, 0.5, 1, 2, 3, 5, 7, 10])
    centre = np.stack([np.asarray(x_over_d, float) * d,
                       np.zeros(len(x_over_d)), np.zeros(len(x_over_d))], axis=1)
    theta = np.radians(np.linspace(0.0, max_angle_deg, n_angles))
    az = math.radians(azimuth_deg)
    rr = radius_over_d * d
    arc = np.stack([rr * np.cos(theta), rr * np.sin(theta) * math.cos(az),
                    rr * np.sin(theta) * math.sin(az)], axis=1)
    keep = arc[:, 0] > 1e-12 * rr
    return {"centreline": centre, "arc": arc[keep],
            "arc_angle_deg": np.degrees(theta[keep])}
