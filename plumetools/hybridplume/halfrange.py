r"""Exact half-range fluxes of a drifting monatomic Maxwellian through a surface.

Used three ways, all independent of the particle solver: to predict what the
source and the interface inject, to predict what returns across the interface,
and to check the injection a DSMC run actually measured.

For a local Maxwellian ``(rho, U, T)`` and a unit normal ``n`` pointing INTO
the region being fed::

    Un = U.n,   Ut = U - Un n,   a = sqrt(2RT),   s = Un/a
    H  = (1 + erf s)/2,           G = exp(-s^2) / (2 sqrt(pi))

    J1 = Un H + a G
    J2 = (Un^2 + a^2/2) H + Un a G
    J3 = (Un^3 + 3 Un a^2/2) H + (Un^2 + a^2) a G

``Jk`` is the k-th moment of the normal velocity over molecules with
``v.n > 0``. The incoming flux densities are::

    mass       rho J1                         (number: n J1)
    momentum   rho (J2 n + J1 Ut)
    energy     (rho/2) (J3 + (|Ut|^2 + 2RT) J1)

The reverse half-range is the same expression with ``n -> -n``. Forward minus
reverse reproduces the *inviscid* equilibrium fluxes ``rho Un``,
``rho U Un + p n`` and ``rho Un (cp T + |U|^2/2)``: a Maxwellian carries no
viscous stress and no conductive heat flux, so matching ``(n, U, T)`` at an
interface does not match Navier--Stokes fluxes there.

Numerical stability
-------------------
For ``s`` well below zero ``Un H`` and ``a G`` nearly cancel. There the moments
are evaluated as

.. math::

    J_k = \frac{a^k e^{-s^2}}{\sqrt\pi}\,(2b)^{-(k+1)}
          \int_0^\infty t^k e^{-t} e^{-t^2/(4b^2)}\,dt, \qquad b = -s,

by Gauss--Laguerre quadrature of an integrand that is positive and smooth, so
nothing cancels. The switch is at ``s = -2``, where the closed form still keeps
about 15 significant figures less one.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import erfc, roots_laguerre

_SQRT_PI = math.sqrt(math.pi)
#: Below this speed ratio the closed form is replaced by the cancellation-free form.
S_SWITCH = -2.0
_LAG_T, _LAG_W = roots_laguerre(80)


def _scaled_moments(s: np.ndarray) -> tuple:
    """Return ``(J1/a, J2/a^2, J3/a^3)`` for arrays of ``s``."""
    s = np.asarray(s, dtype=np.float64)
    j1 = np.empty_like(s)
    j2 = np.empty_like(s)
    j3 = np.empty_like(s)

    closed = s >= S_SWITCH
    if np.any(closed):
        sc = s[closed]
        h = 0.5 * erfc(-sc)          # (1 + erf s)/2 without cancellation
        g = np.exp(-sc * sc) / (2.0 * _SQRT_PI)
        j1[closed] = sc * h + g
        j2[closed] = (sc * sc + 0.5) * h + sc * g
        j3[closed] = (sc ** 3 + 1.5 * sc) * h + (sc * sc + 1.0) * g

    far = ~closed
    if np.any(far):
        b = -s[far]
        damp = np.exp(-(_LAG_T[None, :] ** 2) / (4.0 * b[:, None] ** 2))
        pref = np.exp(-b * b) / _SQRT_PI
        for k, out in ((1, j1), (2, j2), (3, j3)):
            integral = (_LAG_W[None, :] * _LAG_T[None, :] ** k * damp).sum(axis=1)
            out[far] = pref * integral / (2.0 * b) ** (k + 1)
    return j1, j2, j3


def normal_moments(normal_velocity, thermal_speed) -> tuple:
    """``(J1, J2, J3)`` [m/s, m^2/s^2, m^3/s^3] for ``v.n > 0``."""
    un = np.asarray(normal_velocity, dtype=np.float64)
    a = np.asarray(thermal_speed, dtype=np.float64)
    if np.any(~(a > 0)):
        raise ValueError("half-range moments need a positive thermal speed")
    un, a = np.broadcast_arrays(un, a)
    j1, j2, j3 = _scaled_moments(un / a)
    return j1 * a, j2 * a * a, j3 * a ** 3


def incoming_flux(density, velocity, temperature, normal, gas_constant) -> dict:
    """Incoming half-range flux densities through faces with unit ``normal``.

    Args:
        density: ``(n,)`` mass (or number) density.
        velocity: ``(n, 3)`` bulk velocity [m/s].
        temperature: ``(n,)`` translational temperature [K].
        normal: ``(n, 3)`` unit normals pointing into the receiving region.
        gas_constant: ``R`` [J/(kg K)].

    Returns:
        ``{"mass": (n,), "momentum": (n, 3), "energy": (n,)}`` per unit area.
        With a number density in, ``mass`` is a number flux and the other two
        are per molecule mass (multiply by ``m``).
    """
    rho = np.asarray(density, dtype=np.float64)
    u = np.asarray(velocity, dtype=np.float64).reshape(-1, 3)
    t = np.asarray(temperature, dtype=np.float64)
    nhat = np.asarray(normal, dtype=np.float64).reshape(-1, 3)
    if np.any(np.abs(np.linalg.norm(nhat, axis=1) - 1.0) > 1e-9):
        raise ValueError("normals must be unit vectors")
    un = np.einsum("ij,ij->i", u, nhat)
    ut = u - un[:, None] * nhat
    a = np.sqrt(2.0 * gas_constant * t)
    j1, j2, j3 = normal_moments(un, a)
    ut2 = np.einsum("ij,ij->i", ut, ut)
    return {
        "mass": rho * j1,
        "momentum": rho[:, None] * (j2[:, None] * nhat + j1[:, None] * ut),
        "energy": 0.5 * rho * (j3 + (ut2 + 2.0 * gas_constant * t) * j1),
    }


def net_flux(density, velocity, temperature, normal, gas_constant) -> dict:
    """Forward minus reverse half-range transport across ``normal``, per area.

    For a Maxwellian this equals the inviscid continuum flux; see the module
    docstring. Kept as a separate function so tests can check the identity.
    """
    nhat = np.asarray(normal, dtype=np.float64).reshape(-1, 3)
    fwd = incoming_flux(density, velocity, temperature, nhat, gas_constant)
    rev = incoming_flux(density, velocity, temperature, -nhat, gas_constant)
    return {k: fwd[k] - rev[k] for k in fwd}


def integrate(flux: dict, areas) -> dict:
    """Area-integrate per-face flux densities: rates [per s]."""
    area = np.asarray(areas, dtype=np.float64)
    return {
        "mass": float(np.dot(flux["mass"], area)),
        "momentum": (flux["momentum"] * area[:, None]).sum(axis=0).tolist(),
        "energy": float(np.dot(flux["energy"], area)),
    }


def gross_source_number_flux(number_density: float, speed_ratio: float,
                             thermal_speed: float) -> float:
    """Number flux density ``n J1`` a uniform normal-drifting source injects."""
    j1, _, _ = normal_moments(speed_ratio * thermal_speed, thermal_speed)
    return float(number_density * j1)


def bird_eq_4_22_number_flux(number_density, normal_velocity, thermal_speed):
    """Bird eqn 4.22 exactly as ``plumeFieldInflow``/``FreeStream`` accumulate it.

    ``n a (exp(-s^2) + sqrt(pi) s (1 + erf s)) / (2 sqrt(pi))`` -- algebraically
    ``n J1``. Kept as a second, solver-shaped expression so a test can show the
    two agree and the injected count can be predicted from the field files.
    """
    n = np.asarray(number_density, dtype=np.float64)
    a = np.asarray(thermal_speed, dtype=np.float64)
    s = np.asarray(normal_velocity, dtype=np.float64) / a
    return n * a * (np.exp(-s * s) + _SQRT_PI * s * (2.0 - erfc(s))) / (2.0 * _SQRT_PI)
