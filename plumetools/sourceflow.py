"""The analytical plume source-flow model, equations E1-E8.

Pure functions: no I/O, no globals, no hard-coded physical values. Every constant
the model needs is an argument, so a case supplies them from ``case.yaml``.

PROVENANCE -- READ THIS BEFORE TRUSTING THE MODEL
-------------------------------------------------
The source of these equations could not be established. An exhaustive search of
every comment, commit message, and document in this repository turned up exactly
one literature reference, in ``cases/2d-wedge/mesh/generateBlockMeshDict.m``,
covering the *archived* 1d/2d-wedge isentropic-expansion validation case
(Stewart & Lumpkin 2012, via Jonathan S. Pitt, Aegis Aerospace, for LENS).

That reference does **not** cover this model. The limiting angle, the
``(gamma + 0.41)/(gamma - 1)`` angular exponent, and the normalisation integral
appear in no citation anywhere. In particular the ``0.41`` is an unexplained
empirical constant. It is exposed as ``exponent_offset`` so a future correction is
a configuration change rather than a code change.

The equations below were reconstructed from the code, not from a paper. See
``docs/source-flow-model.md`` for the full reconstruction, the confidence
classification, and the open questions.

Expression forms are preserved verbatim from the original -- e.g. ``1 / (ga - 1)``
is not simplified -- because ``ga - 1`` for ``ga = 1.4`` is ``0.3999999999999999``
in binary floating point and rearranging changes the last digits.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate

from plumetools.constants import AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K


def limiting_angle(gamma: float) -> float:
    """E1 -- the limiting turning angle of the plume.

    Args:
        gamma: ratio of specific heats [-].

    Returns:
        theta_l [rad]. 2.27686 rad (130.45 deg) at gamma = 1.4.

    ``theta_l = (pi/2) * (sqrt((gamma+1)/(gamma-1)) - 1)``, the Prandtl-Meyer
    maximum turning angle for expansion to vacuum.

    Note the parenthesisation: commit 22a5fdb corrected this from
    ``sqrt((ga+1)/(ga-1) - 1)`` -- the ``- 1`` was inside the square root. The
    earlier form gave 3.512 rad at gamma = 1.4, so any result produced before
    2022-06-08 used a different limiting angle.
    """
    if gamma <= 1.0:
        raise ValueError(f"gamma must exceed 1, got {gamma}")
    return 0.5 * math.pi * (np.sqrt((gamma + 1) / (gamma - 1)) - 1)


def limiting_velocity(gamma: float, T0_K: float, mass_kg: float) -> float:
    """E2 -- the limiting (vacuum-expansion) speed.

    Args:
        gamma: ratio of specific heats [-].
        T0_K: stagnation temperature [K].
        mass_kg: mass of one molecule [kg].

    Returns:
        v_l [m/s]. 788.164 m/s at gamma = 1.4, T0 = 300 K, N2
        (m = 4.6672e-26 kg); 558.885 m/s for argon at the same temperature.

    ``v_l = sqrt(2 * gamma * k_B * T0 / ((gamma - 1) * m))``.

    Commit 22a5fdb fixed ``k = 1`` here -- Boltzmann's constant was literally 1
    before that date.
    """
    if gamma <= 1.0:
        raise ValueError(f"gamma must exceed 1, got {gamma}")
    if T0_K <= 0.0 or mass_kg <= 0.0:
        raise ValueError(f"T0_K and mass_kg must be positive, got {T0_K}, {mass_kg}")
    return np.sqrt((2 * gamma * BOLTZMANN_J_PER_K * T0_K) / ((gamma - 1) * (mass_kg)))


def angular_theta_only(gamma: float, theta, exponent_offset: float = 0.41):
    """E3 -- the single-angle angular dependence.

    Args:
        gamma: ratio of specific heats [-].
        theta: polar angle [rad].
        exponent_offset: the ``0.41`` in ``(gamma + 0.41)/(gamma - 1)``. Uncited;
            see the module docstring.

    Returns:
        f(theta) [-], 1 on axis and 0 at the limiting angle. Exponent 4.525 at
        gamma = 1.4.

    Used for the normalisation integral (E5). The density itself uses
    :func:`angular_separable` -- which is the inconsistency recorded as SM-05.
    """
    theta_l = limiting_angle(gamma)
    pi = math.pi
    return np.cos(0.5 * pi * theta / theta_l) ** ((gamma + exponent_offset) / (gamma - 1))


def angular_separable(gamma: float, theta, phi, exponent_offset: float = 0.41):
    """E4 -- the separable two-angle dependence actually used for density.

    Args:
        gamma: ratio of specific heats [-].
        theta: polar angle from +z [rad].
        phi: azimuth in the x-y plane from +x [rad].
        exponent_offset: see :func:`angular_theta_only`.

    Returns:
        f(theta, phi) [-].

    ``f = |cos(pi/2 * |theta - pi/2| / theta_l)|^n * cos(pi/2 * phi / theta_l)^n``

    Three deliberate quirks are preserved here:

    * The ``|theta - pi/2|`` shift compensates for theta being measured from +z
      rather than from the plume axis. The result is a product of two angles in
      *different planes*, not a function of the single off-axis angle (SM-04).
    * ``f_theta`` takes ``abs()`` of the cosine, so past the limiting angle the
      profile folds and rises again instead of going to zero (SM-06).
    * ``f_phi`` does **not** take ``abs()``. For ``|phi| > theta_l`` the base is
      negative and the fractional exponent yields **NaN**. That is latent on the
      committed meshes (0 of 2044 faces exceed it) but fires on every face of the
      archived wake-cylinder meshes (99 of 99). Callers should assert finiteness
      before writing -- :func:`number_density` does not do it for them, because
      doing so would change behaviour.
    """
    theta_l = limiting_angle(gamma)
    pi = math.pi
    exponent = (gamma + exponent_offset) / (gamma - 1)
    f_theta = np.abs(np.cos(0.5 * pi * np.abs(theta - 0.5 * pi) / theta_l)) ** exponent
    f_phi = np.cos(0.5 * pi * phi / theta_l) ** exponent
    return f_theta * f_phi


def normalization_coefficient(gamma: float, theta_l: float,
                              integrand_gamma: float = 1.4,
                              exponent_offset: float = 0.41,
                              n_points: int = 500) -> float:
    """E5 -- the normalisation constant A.

    Args:
        gamma: ratio of specific heats [-], used for the numerator.
        theta_l: limiting angle [rad], the upper integration bound.
        integrand_gamma: gamma used *inside* the integrand. See below.
        exponent_offset: see :func:`angular_theta_only`.
        n_points: trapezoid samples over ``[0, theta_l]``.

    Returns:
        A [-].

    ``A = 0.5 * sqrt((gamma-1)/(gamma+1)) / integral(sin(x) * f(x) dx, 0, theta_l)``

    ``integrand_gamma`` defaults to 1.4 because that value was hard-coded inside
    the original integrand, so the ``ga`` argument passed to the legacy
    ``calculateNormCoeff`` was silently ignored. Worse, the integrand normalises
    :func:`angular_theta_only` while the density applies
    :func:`angular_separable` -- so A does not normalise the function it scales,
    and the absolute number density is off by the ratio of the two integrals.
    Both halves of that are finding SM-05, frozen here to preserve behaviour.

    Commit 22a5fdb also fixed the numerator from ``sqrt((ga-1) + (ga+1))`` to
    ``sqrt((ga-1) / (ga+1))``.
    """
    theta = np.linspace(0, theta_l, n_points)
    y = np.sin(theta) * angular_theta_only(integrand_gamma, theta, exponent_offset)
    denom = integrate.trapezoid(y, theta)
    if denom == 0.0:
        raise ValueError("normalisation integral vanished")
    return (0.5 * np.sqrt((gamma - 1) / (gamma + 1))) / denom


def number_density(theta, phi, *, gamma: float, T0_K: float, mass_kg: float,
                   molar_mass_g_per_mol: float, p0_pa: float,
                   throat_radius_m: float, sphere_radius_m: float,
                   exponent_offset: float = 0.41,
                   normalization_gamma: float = 1.4,
                   quadrature_points: int = 500):
    """E6 -- inflow number density.

    Args:
        theta: polar angle from +z [rad].
        phi: azimuth from +x in the x-y plane [rad].
        gamma: ratio of specific heats [-].
        T0_K: stagnation temperature [K], used for the limiting velocity.
        mass_kg: mass of one molecule [kg].
        molar_mass_g_per_mol: molar mass [g/mol], for the mass -> number conversion.
        p0_pa: stagnation (chamber) pressure [Pa].
        throat_radius_m: nozzle throat radius r_e [m]. Named "exit" in the
            original; f2 below is the *sonic* density ratio, which implies the
            throat (finding SM-10, unresolved).
        sphere_radius_m: radius r of the inflow surface [m].
        exponent_offset: see :func:`angular_theta_only`.
        normalization_gamma: see :func:`normalization_coefficient`.
        quadrature_points: see :func:`normalization_coefficient`.

    Returns:
        n [molecules/m^3].

    ``n = (2*A*p0 / v_l^2) * (2/(gamma+1))^(1/(gamma-1)) * (r_e/r)^2
         * f(theta,phi) * (1000*N_A/M_w)``

    Dimensionally: the first factor is Pa/(m^2 s^-2) = kg/m^3, a mass density;
    the middle factors are dimensionless; the last converts kg/m^3 to
    molecules/m^3 for M_w in g/mol. Consistent.

    The result is **not** checked for finiteness -- see :func:`angular_separable`
    for when it can be NaN.
    """
    theta_l = limiting_angle(gamma)
    vel_l = limiting_velocity(gamma, T0_K, mass_kg)
    A = normalization_coefficient(
        gamma, theta_l,
        integrand_gamma=normalization_gamma,
        exponent_offset=exponent_offset,
        n_points=quadrature_points,
    )

    f1 = (2 * A * p0_pa) / (vel_l ** 2)
    f2 = (2 / (gamma + 1)) ** (1 / (gamma - 1))
    f3 = (throat_radius_m / sphere_radius_m) ** 2
    f4 = angular_separable(gamma, theta, phi, exponent_offset)
    conv = (AVOGADRO_PER_MOL * 1000) / molar_mass_g_per_mol
    return f1 * f2 * f3 * f4 * conv


def velocity(v_limit: float, theta, phi) -> np.ndarray:
    """E7 -- inflow velocity components.

    Args:
        v_limit: limiting speed v_l [m/s], from :func:`limiting_velocity`.
        theta: polar angle from +z [rad].
        phi: azimuth from +x in the x-y plane [rad].

    Returns:
        ``(n, 3)`` velocity [m/s].

    ``Ux = v_l cos(phi) sin(theta)``, ``Uy = v_l sin(phi) sin(theta)``,
    ``Uz = v_l cos(theta)`` -- the standard spherical-to-Cartesian conversion with
    polar axis +z.

    Because theta and phi are derived from the same face centroid, this yields
    ``U = v_l * r_hat`` exactly: the magnitude is v_l on every face and the
    direction is exactly radially outward. This part of the model is
    self-consistent by construction, and it is the one piece unaffected by the
    theta-convention problem in SM-04.
    """
    theta = np.asarray(theta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    return np.stack(
        [
            v_limit * np.cos(phi) * np.sin(theta),
            v_limit * np.sin(phi) * np.sin(theta),
            v_limit * np.cos(theta),
        ],
        axis=-1,
    )
