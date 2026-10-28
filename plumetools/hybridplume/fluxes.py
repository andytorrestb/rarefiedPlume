r"""Continuum (Navier--Stokes--Fourier) flux densities through oriented faces.

The CFD precursor solves, for a calorically perfect monatomic gas,

.. math::

    \partial_t\rho + \nabla\cdot(\rho U) = 0, \qquad
    \partial_t(\rho U) + \nabla\cdot(\rho U\otimes U + pI - \tau) = 0,

    \partial_t(\rho E) + \nabla\cdot[(\rho E + p)U - \tau\cdot U + q] = 0,

with ``p = rho R T``, ``E = cv T + |U|^2/2``, ``q = -kappa grad T`` and

.. math::

    \tau = \mu\left[\nabla U + (\nabla U)^T - \tfrac23(\nabla\cdot U) I\right].

Through a face with unit normal ``n`` and ``Un = U.n`` the flux densities are::

    mass      Fm = rho Un
    momentum  Fp = rho U Un + p n - tau.n
    energy    FE = rho Un (cp T + |U|^2/2) - (tau.U).n + q.n

Tensor convention: ``grad_u[..., i, j] = dU_j/dx_i`` -- OpenFOAM's ``grad(U)``
layout. ``tau`` is symmetric, so ``tau.n = n.tau`` and the layout only matters
for building ``tau``; :func:`viscous_stress` transposes explicitly.

Orientation: the caller supplies ``n``. For a CFD-to-DSMC transfer it points
from the CFD region into the DSMC region; for a control-volume budget it is the
outward normal of every boundary face.
"""

from __future__ import annotations

import numpy as np


def viscous_stress(grad_u, viscosity) -> np.ndarray:
    """``tau`` [Pa], ``(n, 3, 3)``, from ``grad_u[i, j] = dU_j/dx_i``."""
    g = np.asarray(grad_u, dtype=np.float64).reshape(-1, 3, 3)
    mu = np.asarray(viscosity, dtype=np.float64).reshape(-1)
    sym = g + np.transpose(g, (0, 2, 1))
    div = np.trace(g, axis1=1, axis2=2)
    tau = sym - (2.0 / 3.0) * div[:, None, None] * np.eye(3)[None, :, :]
    return mu[:, None, None] * tau


def continuum_flux(density, velocity, temperature, normal, gas, *,
                   grad_u=None, grad_t=None, inviscid: bool = False) -> dict:
    """Navier--Stokes--Fourier flux densities through faces with unit ``normal``.

    Args:
        density, velocity, temperature: face states, ``(n,)``, ``(n, 3)``, ``(n,)``.
        normal: ``(n, 3)`` unit normals.
        gas: a :class:`plumetools.hybridplume.physics.Gas`.
        grad_u: ``(n, 3, 3)`` with ``[i, j] = dU_j/dx_i``; needed unless inviscid.
        grad_t: ``(n, 3)`` temperature gradient; needed unless inviscid.
        inviscid: drop ``tau`` and ``q`` (the Euler part only).

    Returns:
        ``{"mass", "momentum", "energy", "pressure", "viscous_momentum",
        "viscous_energy", "conductive_energy"}`` per unit area. The split terms
        are what a Maxwellian inlet cannot carry.
    """
    rho = np.asarray(density, dtype=np.float64).reshape(-1)
    u = np.asarray(velocity, dtype=np.float64).reshape(-1, 3)
    t = np.asarray(temperature, dtype=np.float64).reshape(-1)
    nhat = np.asarray(normal, dtype=np.float64).reshape(-1, 3)
    un = np.einsum("ij,ij->i", u, nhat)
    p = rho * gas.gas_constant * t
    u2 = np.einsum("ij,ij->i", u, u)

    out = {
        "mass": rho * un,
        "pressure": p,
    }
    momentum = rho[:, None] * u * un[:, None] + p[:, None] * nhat
    energy = rho * un * (gas.cp * t + 0.5 * u2)

    zeros3 = np.zeros_like(u)
    zeros = np.zeros_like(rho)
    if inviscid:
        out.update(viscous_momentum=zeros3, viscous_energy=zeros, conductive_energy=zeros)
    else:
        if grad_u is None or grad_t is None:
            raise ValueError("viscous fluxes need grad_u and grad_t (or inviscid=True)")
        tau = viscous_stress(grad_u, gas.viscosity(t))
        tau_n = np.einsum("fij,fj->fi", tau, nhat)
        q = -gas.conductivity(t)[:, None] * np.asarray(grad_t, dtype=np.float64).reshape(-1, 3)
        out["viscous_momentum"] = -tau_n
        out["viscous_energy"] = -np.einsum("fi,fi->f", tau_n, u)
        out["conductive_energy"] = np.einsum("fi,fi->f", q, nhat)
        momentum = momentum + out["viscous_momentum"]
        energy = energy + out["viscous_energy"] + out["conductive_energy"]
    out["momentum"] = momentum
    out["energy"] = energy
    return out
