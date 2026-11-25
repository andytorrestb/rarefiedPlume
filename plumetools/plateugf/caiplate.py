r"""Collisionless loads on a plate normal to Cai's plume (spec §9).

Cai's model: the exit disk ``r <= R`` at ``x = 0`` emits a uniform drifting
Maxwellian (``n0``, ``U0`` along +x, ``T0``) into vacuum, with no collisions. A
point ``P`` then sees, along each direction ``e`` that points back onto the disk,
the exit distribution itself, and nothing along any other direction. With
``beta = m / (2 k T0)`` and ``mu = e_x``,

.. math::

    \langle Q \rangle(P) = n_0 (\beta/\pi)^{3/2}
        \int_{disk} \frac{\mu\, dA}{r^2}\; Q(e)\, M_k(\mu),
    \qquad
    M_k(\mu) = \int_0^\infty c^k e^{-\beta (c^2 - 2 c U_0 \mu + U_0^2)}\, dc ,

so the density is ``k = 2`` with ``Q = 1``, the +x number flux ``k = 3`` with
``Q = mu``, the momentum fluxes ``k = 4`` with ``Q = m mu e_i``, and the energy
flux ``k = 5`` with ``Q = m mu / 2``. Every molecule that reaches a plate at
``x > 0`` moves in +x, so these are the INCIDENT fluxes on its upstream face.

A fully diffuse wall at ``T_w`` re-emits the incident number flux ``Gamma`` as
a half-range Maxwellian, which pushes back with ``m Gamma sqrt(pi k T_w/2m)``,
carries no tangential momentum and takes ``2 k T_w`` per molecule away::

    p   = Pi_xx + m Gamma sqrt(pi k T_w / (2 m))
    tau = (Pi_xy, Pi_xz)
    q   = E - 2 k T_w Gamma

This module integrates the disk numerically (Gauss-Legendre in radius, the
midpoint rule in angle) rather than transcribing Cai and He's closed forms; the
density it gives is held to :func:`plumetools.cai2012.analytical.density_ratio`
by the tests. Gas reflected off the plate never returns to it here: that is the
collisionless limit, which the simulation should approach at large Kn_e.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

BOLTZMANN_J_PER_K = 1.380649e-23
_MU = np.linspace(0.0, 1.0, 4001)


def _radial_table(s: float) -> dict:
    """``M_k(mu) beta^{(k+1)/2}`` for k = 2..5, tabulated on ``_MU``; ``s = sqrt(beta) U0``."""
    t = np.linspace(0.0, s + 12.0, 6000)[None, :]
    mu = _MU[:, None]
    weight = np.exp(-(t * t - 2.0 * t * s * mu + s * s))
    return {k: np.trapezoid(t ** k * weight, t, axis=1) for k in (2, 3, 4, 5)}


def _disk(radius_m: float, n_radial: int, n_azimuthal: int):
    nodes, weights = np.polynomial.legendre.leggauss(int(n_radial))
    rho = 0.5 * radius_m * (nodes + 1.0)
    w_rho = 0.5 * radius_m * weights * rho
    phi = (np.arange(int(n_azimuthal)) + 0.5) * 2.0 * math.pi / n_azimuthal
    w_phi = np.full(int(n_azimuthal), 2.0 * math.pi / n_azimuthal)
    qy = (rho[:, None] * np.cos(phi)[None, :]).ravel()
    qz = (rho[:, None] * np.sin(phi)[None, :]).ravel()
    return qy, qz, (w_rho[:, None] * w_phi[None, :]).ravel()


def moments(points: np.ndarray, *, radius_m: float, number_density: float, speed: float,
            temperature_K: float, mass_kg: float, n_radial: int = 48,
            n_azimuthal: int = 96) -> dict:
    """``n``, ``flux`` (+x number flux), ``pxx``, ``pxy``, ``pxz`` (momentum
    fluxes, Pa) and ``energy`` (W/m^2) of the incident gas at ``points`` (n, 3)."""
    points = np.atleast_2d(np.asarray(points, float))
    if np.any(points[:, 0] <= 0.0):
        raise ValueError("every point must lie downstream of the exit plane (x > 0)")
    beta = mass_kg / (2.0 * BOLTZMANN_J_PER_K * temperature_K)
    table = _radial_table(math.sqrt(beta) * speed)
    pref = number_density * (beta / math.pi) ** 1.5
    qy, qz, w = _disk(radius_m, n_radial, n_azimuthal)
    out = {key: np.empty(len(points)) for key in ("n", "flux", "pxx", "pxy", "pxz", "energy")}
    chunk = max(1, 4_000_000 // w.size)
    for start in range(0, len(points), chunk):
        p = points[start:start + chunk]
        dx = p[:, 0:1]
        dy = p[:, 1:2] - qy[None, :]
        dz = p[:, 2:3] - qz[None, :]
        r2 = dx * dx + dy * dy + dz * dz
        r = np.sqrt(r2)
        mu, ey, ez = dx / r, dy / r, dz / r
        d_omega = mu * w[None, :] / r2
        m = {k: np.interp(mu, _MU, table[k]) * beta ** (-(k + 1) / 2.0) for k in table}
        sl = slice(start, start + len(p))
        out["n"][sl] = pref * (d_omega * m[2]).sum(axis=1)
        out["flux"][sl] = pref * (d_omega * mu * m[3]).sum(axis=1)
        out["pxx"][sl] = mass_kg * pref * (d_omega * mu * mu * m[4]).sum(axis=1)
        out["pxy"][sl] = mass_kg * pref * (d_omega * mu * ey * m[4]).sum(axis=1)
        out["pxz"][sl] = mass_kg * pref * (d_omega * mu * ez * m[4]).sum(axis=1)
        out["energy"][sl] = 0.5 * mass_kg * pref * (d_omega * mu * m[5]).sum(axis=1)
    return out


@dataclass(frozen=True)
class CaiPlate:
    """Collisionless loads on the upstream face of a diffuse plate, per point."""

    n: np.ndarray
    flux: np.ndarray
    pressure: np.ndarray
    shear: np.ndarray        # (n, 2): the y and z components [Pa]
    heat_flux: np.ndarray

    @property
    def shear_magnitude(self) -> np.ndarray:
        return np.hypot(self.shear[:, 0], self.shear[:, 1])


def plate_loads(points: np.ndarray, *, radius_m: float, number_density: float, speed: float,
                temperature_K: float, mass_kg: float, wall_temperature_K: float,
                n_radial: int = 48, n_azimuthal: int = 96) -> CaiPlate:
    """Collisionless pressure, shear and heat flux at plate points (x = stand-off)."""
    m = moments(points, radius_m=radius_m, number_density=number_density, speed=speed,
                temperature_K=temperature_K, mass_kg=mass_kg,
                n_radial=n_radial, n_azimuthal=n_azimuthal)
    kt = BOLTZMANN_J_PER_K * wall_temperature_K
    reemitted = mass_kg * m["flux"] * math.sqrt(math.pi * kt / (2.0 * mass_kg))
    return CaiPlate(
        n=m["n"], flux=m["flux"],
        pressure=m["pxx"] + reemitted,
        shear=np.column_stack([m["pxy"], m["pxz"]]),
        heat_flux=m["energy"] - 2.0 * kt * m["flux"],
    )
