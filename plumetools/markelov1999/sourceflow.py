r"""The source-flow model as printed in AIAA 99-3455.

Model identifier: ``markelov1999_axisymmetric``.

Pure functions. No I/O, no globals, no hard-coded physical values -- every
constant is an argument, supplied by a case's ``case.yaml``.

The equations
-------------

With ``gamma`` the ratio of specific heats, ``k`` Boltzmann's constant, ``T0``
the stagnation temperature, ``m`` the mass of one molecule, ``p0`` the reservoir
pressure, ``r_e`` the physical orifice radius, ``r`` the radial distance from the
orifice centre and ``theta`` the angle from the ``+x`` plume centreline::

    P1   theta_L = pi/2 * (sqrt((gamma + 1)/(gamma - 1)) - 1)

    P2   V       = sqrt(2*gamma*k*T0 / ((gamma - 1)*m))

    P3   f(theta) = cos(pi*theta / (2*theta_L)) ** ((gamma + 0.41)/(gamma - 1))

    P4   A       = 0.5*sqrt((gamma - 1)/(gamma + 1))
                   / integral_0^theta_L [ sin(theta) * f(theta) dtheta ]

    P5   rho(r,theta) = (2*A*p0 / V**2)
                        * (2/(gamma + 1))**(1/(gamma - 1))
                        * (r_e/r)**2
                        * f(theta)

    P6   U(r,theta) = V * r_hat

At ``gamma = 1.4``: ``theta_L`` = 2.276853 rad = 130.454 deg, the angular exponent
is 4.525, the normalisation integral is 0.357656358 and ``A`` = 0.570727014.

``rho`` is a **mass** density [kg/m^3]. DSMC wants a number density, and the
conversion is applied in exactly one place, :func:`number_density`, as
``n = rho / m`` with the same ``m`` that P2 uses. There is no second path.

How this differs from :mod:`plumetools.sourceflow`
--------------------------------------------------

The legacy module is frozen bit-for-bit to reproduce the pre-refactor results,
defects included. It is **not** the same model, and none of the differences below
is a rounding detail:

======  ==========================================================================
SF-1    **Angular dependence.** Legacy ``angular_separable`` returns a *product*
        ``f(theta) * f(phi)`` of two angles measured in different planes. The
        paper prints one angle. This package implements P3 on the single off-axis
        angle and nothing else.
SF-2    **Angle convention.** Legacy ``geometry.spherical`` measures ``theta``
        from ``+z`` and compensates with an ``abs(theta - pi/2)`` shift inside the
        angular function. The plume axis is ``+x``, so this package measures the
        angle directly from ``+x`` -- see :func:`off_axis_angle`.
SF-3    **Orifice radius.** The legacy ``case.yaml`` carries
        ``throat_radius_m: 0.0041275``, ten times the paper's value. The paper's
        orifice diameter is 0.8255 mm, so ``r_e`` = 4.1275e-4 m. Since P5 scales
        as ``r_e**2``, the legacy density is 100x too large.
SF-4    **Reservoir pressure.** Legacy ``calculateRhoN`` hard-coded 475 psi and
        ignored each case's own setting (finding SM-01). Here ``p0_pa`` is an
        argument with no default, so it cannot be silently ignored, and
        :func:`mass_density` is exactly linear in it.
SF-5    **Normalisation gamma.** Legacy ``normalization_coefficient`` takes a
        separate ``integrand_gamma`` defaulting to 1.4, because the original
        hard-coded 1.4 inside the integrand and ignored its own argument (SM-05).
        Worse, it normalised the *theta-only* function while the density applied
        the *separable* one, so ``A`` did not normalise what it scaled. Here P4
        takes one ``gamma`` and normalises exactly the ``f`` that P5 applies.
SF-6    **Beyond the limiting angle.** Legacy ``f_theta`` takes ``abs()`` of the
        cosine so the profile folds and rises again past ``theta_L``, while
        ``f_phi`` does not and yields NaN (SM-06). Here the density is exactly
        zero for ``theta >= theta_L``, which is what "limiting angle" means.
======  ==========================================================================

A seventh difference is not in this module but is worth stating with them: the
reference case mixes N2 source-flow properties with an argon solver setup
(SM-03). This case family is N2 throughout, and
:func:`plumetools.markelov1999.foamdicts.render_dsmc_properties` writes the same
species the model was evaluated for.

The ``0.41``
------------

``(gamma + 0.41)/(gamma - 1)`` is the exponent printed in the paper. The ``0.41``
is an empirical constant with no derivation given there either. It is exposed as
``exponent_offset`` so that a future correction is a configuration change rather
than a code change -- the same treatment the legacy module gives it.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate

from plumetools.markelov1999.constants import BOLTZMANN_J_PER_K

#: The ``model:`` value in ``case.yaml`` that selects this implementation.
MODEL_NAME = "markelov1999_axisymmetric"

#: The paper's orifice diameter [m]: 0.8255 mm.
ORIFICE_DIAMETER_M = 0.8255e-3

#: The paper's orifice radius [m]: half of :data:`ORIFICE_DIAMETER_M`.
ORIFICE_RADIUS_M = ORIFICE_DIAMETER_M / 2.0


def limiting_angle(gamma: float) -> float:
    """P1 -- the limiting turning angle of the plume [rad].

    Args:
        gamma: ratio of specific heats [-].

    Returns:
        ``theta_L`` [rad]; 2.276853 rad = 130.454 deg at ``gamma = 1.4``.

    Raises:
        ValueError: if ``gamma <= 1``, where the expression is undefined.

    This is the Prandtl-Meyer maximum turning angle for expansion to vacuum.
    ``theta_L > pi/2`` exactly when ``gamma < 5/3``: the condition reduces to
    ``(gamma + 1)/(gamma - 1) > 4``. So for the diatomic nitrogen of this case
    family the cone reaches 130 deg and a hemispherical inflow surface -- which
    spans at most 90 deg off axis -- never reaches the clip in :func:`angular`.
    For a monatomic gas the two coincide at exactly 90 deg, and for a full-sphere
    source boundary the clip does fire. It is therefore not dead code, and
    :mod:`plumetools.markelov1999.checks` verifies the property on the actual
    meshed patch rather than assuming it.
    """
    if gamma <= 1.0:
        raise ValueError(f"gamma must exceed 1, got {gamma}")
    return 0.5 * math.pi * (math.sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0)


def limiting_velocity(gamma: float, T0_K: float, mass_kg: float) -> float:
    """P2 -- the limiting (vacuum-expansion) speed [m/s].

    Args:
        gamma: ratio of specific heats [-].
        T0_K: stagnation temperature [K].
        mass_kg: mass of one molecule [kg].

    Returns:
        ``V`` [m/s]; 789.4849 m/s for N2 at ``gamma = 1.4``, ``T0 = 300 K``.

    Raises:
        ValueError: for a non-physical ``gamma``, ``T0_K`` or ``mass_kg``.

    The legacy module reports 788.164 m/s for the same inputs. The difference is
    entirely its 5 s.f. Boltzmann constant and its ``amu``-based molecular mass;
    see :mod:`plumetools.markelov1999.constants`.
    """
    if gamma <= 1.0:
        raise ValueError(f"gamma must exceed 1, got {gamma}")
    if T0_K <= 0.0:
        raise ValueError(f"T0_K must be positive, got {T0_K}")
    if mass_kg <= 0.0:
        raise ValueError(f"mass_kg must be positive, got {mass_kg}")
    return math.sqrt(2.0 * gamma * BOLTZMANN_J_PER_K * T0_K / ((gamma - 1.0) * mass_kg))


def angular(gamma: float, theta, exponent_offset: float = 0.41):
    """P3 -- the angular density profile ``f(theta)``, clipped at ``theta_L``.

    Args:
        gamma: ratio of specific heats [-].
        theta: off-axis angle from the ``+x`` plume centreline [rad]. Scalar or
            array. Negative values are treated by magnitude, so the profile is
            symmetric about the axis as an axisymmetric model requires.
        exponent_offset: the ``0.41`` in ``(gamma + 0.41)/(gamma - 1)``.

    Returns:
        ``f(theta)`` [-], an array of the same shape as ``theta``. Exactly 1 on
        the axis, falling monotonically to exactly 0 at and beyond ``theta_L``.

    One angle, measured from the plume axis (SF-1, SF-2). Zero -- not a folded
    cosine, and never NaN -- outside the limiting angle (SF-6).
    """
    theta_l = limiting_angle(gamma)
    t = np.abs(np.asarray(theta, dtype=np.float64))
    exponent = (gamma + exponent_offset) / (gamma - 1.0)

    # Outside the cone the cosine base turns negative, and a negative base with a
    # fractional exponent is NaN -- precisely the legacy failure mode (SF-6). So
    # the base is forced to zero there before the power is taken, and 0**4.525 is
    # a clean 0. Inside the cone the argument stays in [0, pi/2) and the cosine is
    # non-negative by construction.
    inside = t < theta_l
    base = np.where(inside, np.cos(0.5 * math.pi * t / theta_l), 0.0)
    return np.where(inside, base ** exponent, 0.0)


def normalization_integral(gamma: float, exponent_offset: float = 0.41,
                           quadrature_points: int = 0) -> float:
    """The denominator of P4: ``integral_0^theta_L sin(theta) f(theta) dtheta``.

    Args:
        gamma: ratio of specific heats [-].
        exponent_offset: see :func:`angular`.
        quadrature_points: 0 (the default) uses adaptive Gauss-Kronrod, which is
            accurate to ~1e-14 here. A positive value uses that many uniformly
            spaced trapezoid samples instead, which exists only so a test can
            show the two agree.

    Returns:
        The integral [-]; 0.357656358 at ``gamma = 1.4``.

    Raises:
        ValueError: if the integral vanishes, which would make ``A`` infinite.

    The integrand uses the *same* ``gamma`` as the caller and the *same* ``f`` the
    density applies. Both were untrue in the legacy code (SF-5).
    """
    theta_l = limiting_angle(gamma)

    def integrand(t: float) -> float:
        return math.sin(t) * float(angular(gamma, t, exponent_offset))

    if quadrature_points > 0:
        theta = np.linspace(0.0, theta_l, int(quadrature_points))
        y = np.sin(theta) * angular(gamma, theta, exponent_offset)
        value = float(integrate.trapezoid(y, theta))
    else:
        value, _err = integrate.quad(integrand, 0.0, theta_l, limit=200)

    if value <= 0.0:
        raise ValueError(
            f"normalisation integral is {value}, so A would not be finite; "
            f"gamma={gamma}, exponent_offset={exponent_offset}"
        )
    return float(value)


def normalization_coefficient(gamma: float, exponent_offset: float = 0.41,
                              quadrature_points: int = 0) -> float:
    """P4 -- the normalisation constant ``A`` [-].

    Args:
        gamma: ratio of specific heats [-].
        exponent_offset: see :func:`angular`.
        quadrature_points: see :func:`normalization_integral`.

    Returns:
        ``A`` [-]; 0.570727014 at ``gamma = 1.4``.

    ``A = 0.5*sqrt((gamma-1)/(gamma+1)) / integral_0^theta_L sin(t) f(t) dt``.
    """
    numerator = 0.5 * math.sqrt((gamma - 1.0) / (gamma + 1.0))
    return numerator / normalization_integral(gamma, exponent_offset, quadrature_points)


def sonic_density_ratio(gamma: float) -> float:
    """``(2/(gamma + 1))**(1/(gamma - 1))``, the sonic-throat density ratio [-].

    0.633938145 at ``gamma = 1.4``. Broken out of P5 because it is the one factor
    there that depends on nothing but ``gamma``, so it can be checked on its own.
    """
    if gamma <= 1.0:
        raise ValueError(f"gamma must exceed 1, got {gamma}")
    return (2.0 / (gamma + 1.0)) ** (1.0 / (gamma - 1.0))


def mass_density(theta, radius_m, *, gamma: float, T0_K: float, mass_kg: float,
                 p0_pa: float, orifice_radius_m: float = ORIFICE_RADIUS_M,
                 exponent_offset: float = 0.41, quadrature_points: int = 0):
    """P5 -- source-flow mass density [kg/m^3].

    Args:
        theta: off-axis angle from ``+x`` [rad]; scalar or array.
        radius_m: distance from the orifice centre [m]; scalar or array
            broadcastable against ``theta``.
        gamma: ratio of specific heats [-].
        T0_K: stagnation temperature [K].
        mass_kg: mass of one molecule [kg].
        p0_pa: reservoir (stagnation) pressure [Pa]. Required -- there is no
            default to fall back to (SF-4).
        orifice_radius_m: the **physical** orifice radius ``r_e`` [m]. Defaults to
            the paper's 4.1275e-4 m. This is *not* the radius of the hemispherical
            DSMC inflow boundary; see
            :class:`plumetools.markelov1999.geometry.MarkelovGeometry`.
        exponent_offset: see :func:`angular`.
        quadrature_points: see :func:`normalization_integral`.

    Returns:
        ``rho`` [kg/m^3], zero outside the limiting angle.

    Raises:
        ValueError: for non-positive ``p0_pa``, ``orifice_radius_m`` or
            ``radius_m``.

    Dimensions: ``2*A*p0/V**2`` is Pa / (m/s)^2 = kg/m^3; every other factor is
    dimensionless. So the result is a mass density and P5 is dimensionally
    consistent as printed.
    """
    if p0_pa <= 0.0:
        raise ValueError(f"p0_pa must be positive, got {p0_pa}")
    if orifice_radius_m <= 0.0:
        raise ValueError(f"orifice_radius_m must be positive, got {orifice_radius_m}")
    r = np.asarray(radius_m, dtype=np.float64)
    if np.any(r <= 0.0):
        raise ValueError("radius_m must be positive everywhere (the source is singular at r=0)")

    V = limiting_velocity(gamma, T0_K, mass_kg)
    A = normalization_coefficient(gamma, exponent_offset, quadrature_points)

    return (
        (2.0 * A * p0_pa / (V * V))
        * sonic_density_ratio(gamma)
        * (orifice_radius_m / r) ** 2
        * angular(gamma, theta, exponent_offset)
    )


def number_density(theta, radius_m, *, gamma: float, T0_K: float, mass_kg: float,
                   p0_pa: float, orifice_radius_m: float = ORIFICE_RADIUS_M,
                   exponent_offset: float = 0.41, quadrature_points: int = 0):
    """Source-flow number density [molecules/m^3].

    Arguments are those of :func:`mass_density`.

    THE conversion, applied in exactly one place::

        n = rho / m

    with the same ``m`` (``mass_kg``) that :func:`limiting_velocity` uses. The
    legacy model instead multiplies by ``1000 * N_A / M_w``, an algebraically
    equivalent but independently-rounded second path; keeping one path means a
    mass-density check and a number-density check cannot disagree.
    """
    return mass_density(
        theta, radius_m,
        gamma=gamma, T0_K=T0_K, mass_kg=mass_kg, p0_pa=p0_pa,
        orifice_radius_m=orifice_radius_m,
        exponent_offset=exponent_offset,
        quadrature_points=quadrature_points,
    ) / mass_kg


def off_axis_angle(points, plume_axis=(1.0, 0.0, 0.0)) -> np.ndarray:
    """The single angle from the plume centreline [rad], for each point.

    Args:
        points: ``(n, 3)`` positions relative to the orifice centre [m].
        plume_axis: the plume centreline direction; ``+x`` for this case family.
            Normalised internally, so it need not be a unit vector.

    Returns:
        ``(n,)`` angles in ``[0, pi]``, where 0 is straight down the axis.

    Raises:
        ValueError: for a degenerate axis, a badly shaped array, or a point at
            the origin (where the angle is undefined).

    ``theta = arccos((p . a_hat) / |p|)``. This is measured directly from the
    plume axis, which is what P3 and P5 want -- the legacy convention measures it
    from ``+z`` and then patches the discrepancy inside the angular function
    (SF-2).
    """
    p = np.asarray(points, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"expected (n, 3) points, got {p.shape}")

    axis = np.asarray(plume_axis, dtype=np.float64)
    axis_norm = np.linalg.norm(axis)
    if axis_norm == 0.0:
        raise ValueError("plume_axis must be non-zero")
    axis = axis / axis_norm

    r = np.linalg.norm(p, axis=1)
    if np.any(r <= 0.0):
        raise ValueError("the off-axis angle is undefined at the orifice centre")

    # Clip before arccos: |cos| can exceed 1 by ~1e-16 from round-off.
    return np.arccos(np.clip((p @ axis) / r, -1.0, 1.0))


def radial_velocity(v_limit: float, points) -> np.ndarray:
    """P6 -- inflow velocity [m/s], radial and of magnitude ``v_limit``.

    Args:
        v_limit: the limiting speed ``V`` [m/s], from :func:`limiting_velocity`.
        points: ``(n, 3)`` positions relative to the orifice centre [m].

    Returns:
        ``(n, 3)`` velocities [m/s], each ``v_limit * r_hat``.

    Raises:
        ValueError: for a badly shaped array or a point at the origin.

    A source flow is radial by definition, so this is ``V`` times the unit
    position vector -- no spherical round trip, and no dependence on how the angle
    happens to be parameterised.
    """
    p = np.asarray(points, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"expected (n, 3) points, got {p.shape}")
    r = np.linalg.norm(p, axis=1)
    if np.any(r <= 0.0):
        raise ValueError("the radial direction is undefined at the orifice centre")
    return float(v_limit) * p / r[:, None]
