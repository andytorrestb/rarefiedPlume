r"""Cai's collisionless solution for a circular exit, as the comparison solution.

Cai, C., and Wang, L., "Numerical Validations for a Set of Collisionless Rocket
Plume Solutions," *J. Spacecraft and Rockets*, Vol. 49, No. 1, 2012.

.. warning::

   The closed forms below are **re-derived from the underlying free-molecular
   integral**, not transcribed from the paper. Equation numbers are cited as
   pointers to where the same results appear in Cai (density field, Eq. 5;
   centreline quantities, Eqs. 18-21), *not* as a claim that the algebra here is
   character-for-character his. The derivation is written out in full so the
   agreement can be checked, and every limit that can be verified independently
   is asserted in ``tests/unit/test_cai2012_analytical.py``.

   Nothing here is digitised from a plot. §10 of the case specification forbids
   it, and a digitised curve cannot be differentiated, integrated or extended.

The physical problem
--------------------
A circular exit of radius :math:`R_0` occupies the disk :math:`r \le R_0` in the
plane :math:`x = 0`. It emits a uniform drifting Maxwellian

.. math::

    f_0(\mathbf{v}) = n_0\left(\frac{\beta_0}{\pi}\right)^{3/2}
                      \exp\!\left[-\beta_0\,|\mathbf{v} - U_0\hat{x}|^2\right],
    \qquad \beta_0 = \frac{1}{2RT_0},\quad S_0 = U_0\sqrt{\beta_0},

into a vacuum, and there are no collisions. Collisionless means the distribution
at a field point :math:`P` is simply :math:`f_0` restricted to those velocity
directions whose backward ray strikes the disk:

.. math::

    f(P,\mathbf{v}) = \begin{cases}
        f_0(\mathbf{v}) & \text{if } P - t\hat{v} \text{ meets the disk for some } t>0\\
        0 & \text{otherwise.}
    \end{cases}

Every moment is therefore an integral of :math:`f_0` over the **solid angle the
exit disk subtends at** :math:`P`.

The centreline
--------------
On the axis the subtended solid angle is a right circular cone of half-angle
:math:`\theta_{max} = \arctan(R_0/x)`, so with :math:`t = \cos\theta` and

.. math::

    t_0 \equiv \cos\theta_{max} = \frac{x}{\sqrt{x^2 + R_0^2}}
    \tag{C0}

every moment collapses to a one-dimensional integral in :math:`t` over
:math:`[t_0, 1]`. Writing :math:`s = S_0 t` and
:math:`G(t) = e^{S_0^2t^2}\left(1+\mathrm{erf}(S_0 t)\right)`, the three
velocity-moment integrals

.. math::

    \int_0^\infty w^2 e^{-w^2+2ws}\,dw &= \tfrac12 s
        + \tfrac{\sqrt\pi}{4}e^{s^2}(1+2s^2)(1+\mathrm{erf}\,s)\\
    \int_0^\infty w^3 e^{-w^2+2ws}\,dw &= \tfrac12(1+s^2)
        + \tfrac{\sqrt\pi}{4}(3s+2s^3)e^{s^2}(1+\mathrm{erf}\,s)\\
    \int_0^\infty w^4 e^{-w^2+2ws}\,dw &= \tfrac12 s^3 + \tfrac54 s
        + \sqrt\pi\left(\tfrac38+\tfrac32 s^2+\tfrac12 s^4\right)
          e^{s^2}(1+\mathrm{erf}\,s)

integrate in :math:`t` **exactly**, because in each case the polynomial part
cancels against the :math:`2a/\sqrt\pi` term produced by differentiating
:math:`t^kG(t)`:

.. math::

    \frac{d}{dt}\left[t^k G(t)\right]
        = \left(k t^{k-1} + 2S_0^2 t^{k+1}\right)G(t)
          + \frac{2S_0}{\sqrt\pi}t^k .

What survives is, with :math:`E_0 \equiv e^{-S_0^2(1-t_0^2)}`,
:math:`A \equiv 1+\mathrm{erf}\,S_0` and :math:`A_0 \equiv 1+\mathrm{erf}(S_0t_0)`:

.. math::

    \frac{n}{n_0} &= \tfrac12\left[A - t_0 E_0 A_0\right]
    \tag{C1}\\
    \frac{n}{n_0}\,U_1\sqrt{\beta_0}
        &= \frac{e^{-S_0^2}(1-t_0^2)}{2\sqrt\pi}
           + \frac{S_0}{2}\left[A - t_0^3 E_0 A_0\right]
    \tag{C2}\\
    \frac{n}{n_0}\,\beta_0\langle v^2\rangle
        &= \frac{S_0 e^{-S_0^2}(1-t_0^2)}{2\sqrt\pi}
           + 2\left(\tfrac38+\tfrac{S_0^2}{4}\right)A
           - 2\left(\tfrac38 t_0+\tfrac{S_0^2}{4}t_0^3\right)E_0 A_0
    \tag{C3}

and the temperature follows from its definition
:math:`3RT = \langle v^2\rangle - U_1^2`:

.. math::

    \frac{T}{T_0} = \frac23\left[\beta_0\langle v^2\rangle
                                 - \left(U_1\sqrt{\beta_0}\right)^2\right].
    \tag{C4}

Limits worth knowing, all asserted in the tests:

======================  ==========================================================
:math:`x \to 0`         :math:`n/n_0 \to \tfrac12(1+\mathrm{erf}\,S_0)`. **Not
                        exactly 1**: only the outgoing half of the distribution
                        is present at the exit plane. It approaches 1 as
                        :math:`S_0\to\infty` and is 0.99766 at :math:`S_0=2`,
                        and exactly 1/2 at :math:`S_0=0`.
:math:`x\to\infty`      :math:`n/n_0 \to 0` like :math:`x^{-2}`;
                        :math:`U_1\sqrt{\beta_0}` **rises** to a finite limit,
                        because the slow molecules fall out of the shrinking
                        visible cone first.
:math:`S_0=0,\;x=0`     :math:`n/n_0 = 1/2`, :math:`U_1\sqrt{\beta_0}=1/\sqrt\pi`
                        (the half-Maxwellian mean :math:`\sqrt{2kT/\pi m}`), and
                        :math:`T/T_0 = 1-2/(3\pi) = 0.78779`.
======================  ==========================================================

The full field
--------------
Off the axis the subtended solid angle is not a cone and no such collapse
happens. The density is instead integrated over the exit disk directly, which is
the same solution Cai writes in closed form as his Eq. 5:

.. math::

    \frac{n(P)}{n_0} = \pi^{-3/2}\int_{disk}
        \frac{\cos\psi}{|P-Q|^2}
        \left[\tfrac12 s\,e^{-S_0^2}
              + \tfrac{\sqrt\pi}{4}e^{-S_0^2(1-\cos^2\psi)}
                (1+2s^2)(1+\mathrm{erf}\,s)\right] dA,
    \tag{C5}

with :math:`\hat e = (P-Q)/|P-Q|`, :math:`\cos\psi = \hat e\cdot\hat x = x/|P-Q|`
and :math:`s = S_0\cos\psi`. Evaluated by Gauss-Legendre quadrature on the disk.
Cai's closed form and this quadrature are the same function; the test that they
agree is :func:`density_ratio` against (C1) on the axis, to 1e-10.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import erf

__all__ = [
    "cos_theta_max",
    "centerline_density_ratio",
    "centerline_axial_speed_ratio",
    "centerline_temperature_ratio",
    "centerline_profile",
    "density_ratio",
    "moments",
]

_SQRT_PI = math.sqrt(math.pi)


def _as_array(value):
    return np.asarray(value, dtype=np.float64)


def cos_theta_max(x, nozzle_radius_m: float):
    r"""Equation (C0): ``t0 = x / sqrt(x^2 + R0^2)``.

    The cosine of the half-angle of the cone the exit disk subtends on the axis.
    ``0`` at the exit plane, ``-> 1`` far downstream.

    Args:
        x: axial position [m]; scalar or array. Must be ``>= 0``.
        nozzle_radius_m: ``R0`` [m].

    Raises:
        ValueError: for a non-positive radius or any ``x < 0``. Upstream of the
            exit plane there is no plume: the solution is one-sided by
            construction and continuing it backwards would return a number with
            no meaning.
    """
    if not nozzle_radius_m > 0.0:
        raise ValueError(f"nozzle radius must be positive, got {nozzle_radius_m}")
    x = _as_array(x)
    if np.any(x < 0.0):
        raise ValueError(
            "the collisionless solution is defined for x >= 0 only; the exit "
            "plane is at x = 0 and there is no upstream branch")
    return x / np.sqrt(x * x + float(nozzle_radius_m) ** 2)


def _centerline_moments(x, nozzle_radius_m: float, speed_ratio: float):
    """``(D, Mx, E)`` -- the right-hand sides of (C1), (C2) and (C3).

    Grouped in one function because all three share ``t0``, ``E0`` and ``A0``,
    and because the ratios that form the velocity and the temperature must be
    taken from moments evaluated at exactly the same ``t0``.
    """
    s0 = float(speed_ratio)
    t0 = cos_theta_max(x, nozzle_radius_m)
    t0sq = t0 * t0
    one_minus = 1.0 - t0sq

    a_full = 1.0 + erf(s0)
    a_cone = 1.0 + erf(s0 * t0)
    # exp(-S0^2 (1 - t0^2)) rather than exp(S0^2 t0^2) * exp(-S0^2): the second
    # form overflows for a large speed ratio while the product does not.
    e0 = np.exp(-(s0 * s0) * one_minus)

    density = 0.5 * (a_full - t0 * e0 * a_cone)
    momentum = (math.exp(-s0 * s0) * one_minus / (2.0 * _SQRT_PI)
                + 0.5 * s0 * (a_full - t0 * t0sq * e0 * a_cone))
    energy = (s0 * math.exp(-s0 * s0) * one_minus / (2.0 * _SQRT_PI)
              + 2.0 * (0.375 + 0.25 * s0 * s0) * a_full
              - 2.0 * (0.375 * t0 + 0.25 * s0 * s0 * t0 * t0sq) * e0 * a_cone)
    return density, momentum, energy


def centerline_density_ratio(x, nozzle_radius_m: float, speed_ratio: float):
    """Equation (C1): ``n/n0`` on the plume axis. Cai Fig. 19-21, Eq. 18.

    Args:
        x: axial position [m], scalar or array.
        nozzle_radius_m: ``R0`` [m].
        speed_ratio: ``S0``.

    Returns:
        ``n/n0``, the same shape as ``x``.
    """
    density, _, _ = _centerline_moments(x, nozzle_radius_m, speed_ratio)
    return density


def centerline_axial_speed_ratio(x, nozzle_radius_m: float, speed_ratio: float):
    """``U1 sqrt(beta0)`` on the axis -- (C2) divided by (C1). Cai Eq. 19-20.

    The quantity Cai plots: a velocity normalised by ``sqrt(2 R T0)``, so it is
    the *local speed ratio*. It starts near ``S0`` at the exit and **increases**
    downstream, because the visible cone narrows and the molecules that leave it
    first are the slow ones.
    """
    density, momentum, _ = _centerline_moments(x, nozzle_radius_m, speed_ratio)
    return momentum / density


def centerline_temperature_ratio(x, nozzle_radius_m: float, speed_ratio: float):
    """Equation (C4): ``T/T0`` on the axis. Cai Eq. 21.

    ``T`` is the full (translational) temperature from
    ``3 R T = <v^2> - U1^2``, not a per-component one. A collisionless expansion
    has no mechanism to equilibrate the components, so ``T_parallel`` and
    ``T_perpendicular`` diverge downstream and only their average is this.
    """
    density, momentum, energy = _centerline_moments(x, nozzle_radius_m, speed_ratio)
    speed = momentum / density
    return (2.0 / 3.0) * (energy / density - speed * speed)


def centerline_profile(x, nozzle_radius_m: float, speed_ratio: float) -> dict:
    """All three centreline quantities at once, plus ``X/D``.

    Returns:
        ``{"x_m", "x_over_D", "n_over_n0", "U_sqrt_beta0", "T_over_T0"}``, each a
        float64 array. One call rather than three so the moments are evaluated
        once and the three curves are guaranteed to come from the same ``t0``.
    """
    density, momentum, energy = _centerline_moments(x, nozzle_radius_m, speed_ratio)
    speed = momentum / density
    return {
        "x_m": _as_array(x),
        "x_over_D": _as_array(x) / (2.0 * float(nozzle_radius_m)),
        "n_over_n0": density,
        "U_sqrt_beta0": speed,
        "T_over_T0": (2.0 / 3.0) * (energy / density - speed * speed),
    }


# --------------------------------------------------------------------------- #
# the full density field
# --------------------------------------------------------------------------- #

#: Ceilings on the adaptive quadrature order, so a field point arbitrarily close
#: to the exit plane cannot ask for an arbitrarily expensive integral.
MAX_RADIAL = 512
MAX_AZIMUTHAL = 2048


def _refined_order(nozzle_radius_m: float, min_x: float,
                   n_radial: int, n_azimuthal: int) -> tuple:
    r"""Raise the quadrature order for field points close to the exit plane.

    The integrand carries a :math:`1/R^2` peak at the disk point directly
    opposite the field point, and that peak is only as wide as the field point's
    distance from the plane. At :math:`x = 0.05 R_0` a fixed 48 x 96 rule
    straddles it and returns :math:`n/n_0 = 1.000001` -- above the physical bound
    of :math:`\tfrac12(1+\mathrm{erf}\,S_0)`, which is how the problem announces
    itself.

    Resolving the peak with about two nodes across its half-width needs

    .. math::

        n_\rho \gtrsim \frac{\pi R_0}{x}, \qquad
        n_\phi \gtrsim \frac{4\pi R_0}{x}

    (the azimuthal arc length at :math:`\rho = R_0` is :math:`2\pi R_0/n_\phi`).
    The order is never *lowered*, and is capped by :data:`MAX_RADIAL` and
    :data:`MAX_AZIMUTHAL`.
    """
    # The requested order is validated here, before refinement can hide a
    # nonsensical one behind a bump.
    if int(n_radial) < 2 or int(n_azimuthal) < 4:
        raise ValueError(
            f"quadrature too coarse: n_radial={n_radial} (>= 2), "
            f"n_azimuthal={n_azimuthal} (>= 4)")
    if not min_x > 0.0:
        return int(n_radial), int(n_azimuthal)
    radius = float(nozzle_radius_m)
    needed_radial = int(math.ceil(2.0 * math.pi * radius / min_x))
    needed_azimuthal = int(math.ceil(4.0 * math.pi * radius / min_x))
    return (min(MAX_RADIAL, max(int(n_radial), needed_radial)),
            min(MAX_AZIMUTHAL, max(int(n_azimuthal), needed_azimuthal)))


def _disk_quadrature(nozzle_radius_m: float, n_radial: int, n_azimuthal: int):
    """Gauss-Legendre nodes and weights on the exit disk.

    Radially Gauss-Legendre in ``rho`` with the Jacobian ``rho`` folded into the
    weight; azimuthally the midpoint rule, which is spectrally accurate for a
    periodic integrand and needs no endpoint handling.
    """
    if n_radial < 2 or n_azimuthal < 4:
        raise ValueError(
            f"quadrature too coarse: n_radial={n_radial} (>= 2), "
            f"n_azimuthal={n_azimuthal} (>= 4)")
    nodes, weights = np.polynomial.legendre.leggauss(int(n_radial))
    radius = float(nozzle_radius_m)
    rho = 0.5 * radius * (nodes + 1.0)
    w_rho = 0.5 * radius * weights * rho

    phi = (np.arange(int(n_azimuthal)) + 0.5) * (2.0 * math.pi / int(n_azimuthal))
    w_phi = np.full(int(n_azimuthal), 2.0 * math.pi / int(n_azimuthal))

    return rho, w_rho, phi, w_phi


def density_ratio(x, y, z, nozzle_radius_m: float, speed_ratio: float, *,
                  n_radial: int = 64, n_azimuthal: int = 128):
    """Equation (C5): ``n/n0`` anywhere in the plume. Cai Eq. 5.

    Args:
        x, y, z: field-point coordinates [m]. Broadcast against each other, so a
            plane can be passed as two meshgrids and a scalar.
        nozzle_radius_m: ``R0`` [m].
        speed_ratio: ``S0``.
        n_radial: Gauss-Legendre points across the disk radius.
        n_azimuthal: midpoint-rule points around the disk.

    Returns:
        ``n/n0`` with the broadcast shape of ``(x, y, z)``.

    Raises:
        ValueError: if any ``x <= 0``. The integrand has an inverse-square
            singularity when the field point lies **in** the exit plane, where
            the quadrature is meaningless; the exit plane itself is covered
            exactly by (C1) on the axis and is a boundary condition everywhere
            else. In practice the first cell centre of any mesh is at
            ``x = dx/2 > 0``.

    Memory: the quadrature array is ``points x n_radial x n_azimuthal``, so a
    400 x 400 plane at the default order is ~2.6e9 entries. It is evaluated in
    chunks of field points to keep that bounded.
    """
    x = _as_array(x)
    y = _as_array(y)
    z = _as_array(z)
    x, y, z = np.broadcast_arrays(x, y, z)
    if np.any(x <= 0.0):
        raise ValueError(
            "density_ratio needs x > 0: in the exit plane the disk integrand "
            "has an inverse-square singularity at the field point itself. Use "
            "centerline_density_ratio for the axis, which is closed form.")

    s0 = float(speed_ratio)
    n_radial, n_azimuthal = _refined_order(
        nozzle_radius_m, float(np.min(x)), n_radial, n_azimuthal)
    rho, w_rho, phi, w_phi = _disk_quadrature(
        nozzle_radius_m, n_radial, n_azimuthal)

    # Source points on the disk, and the outer-product weight.
    qy = (rho[:, None] * np.cos(phi)[None, :]).ravel()
    qz = (rho[:, None] * np.sin(phi)[None, :]).ravel()
    weight = (w_rho[:, None] * w_phi[None, :]).ravel()

    shape = x.shape
    xf, yf, zf = x.ravel(), y.ravel(), z.ravel()
    out = np.empty(xf.shape, dtype=np.float64)

    # ~4e6 (field point, quadrature node) pairs per chunk: enough to keep numpy
    # busy, small enough that the temporaries stay well under a gigabyte.
    chunk = max(1, int(4_000_000 // max(1, weight.size)))
    for start in range(0, xf.size, chunk):
        stop = min(start + chunk, xf.size)
        dx = xf[start:stop, None]
        dy = yf[start:stop, None] - qy[None, :]
        dz = zf[start:stop, None] - qz[None, :]

        r2 = dx * dx + dy * dy + dz * dz
        r = np.sqrt(r2)
        cos_psi = dx / r
        s = s0 * cos_psi

        # exp(s^2) * exp(-S0^2) written as one exponential, as in the centreline
        # moments: the factors overflow separately for a large speed ratio.
        bracket = (0.5 * s * math.exp(-s0 * s0)
                   + 0.25 * _SQRT_PI * np.exp(-(s0 * s0) * (1.0 - cos_psi * cos_psi))
                   * (1.0 + 2.0 * s * s) * (1.0 + erf(s)))
        out[start:stop] = (
            math.pi ** -1.5 * ((cos_psi / r2) * bracket * weight[None, :]).sum(axis=1))

    return out.reshape(shape)


def moments(x, y, z, nozzle_radius_m: float, speed_ratio: float, *,
            n_radial: int = 48, n_azimuthal: int = 96) -> dict:
    r"""Density, velocity **and** temperature anywhere in the plume.

    Args:
        x, y, z: field-point coordinates [m], broadcast against each other.
        nozzle_radius_m: ``R0`` [m].
        speed_ratio: ``S0``.
        n_radial, n_azimuthal: disk quadrature order.

    Returns:
        ``{"n_over_n0", "U_sqrt_beta0", "T_over_T0"}``. ``U_sqrt_beta0`` has a
        trailing axis of length 3 -- off the plume axis the velocity has
        transverse components, and dropping them would make the temperature
        wrong, because ``T`` is defined about the **local** mean velocity.

    The same disk integral as :func:`density_ratio`, taken with the next two
    velocity moments. From the three :math:`w`-integrals in the module docstring,
    with :math:`\hat e` the unit vector from the source point to the field point,
    :math:`\cos\psi = \hat e\cdot\hat x` and :math:`s = S_0\cos\psi`:

    .. math::

        \frac{n}{n_0} &= \pi^{-3/2}\!\int \frac{\cos\psi}{R^2}
            \left[\tfrac12 s\,e^{-S_0^2} + \tfrac{\sqrt\pi}{4}
                  \mathcal{E}(1+2s^2)(1+\mathrm{erf}\,s)\right] dA \\
        \frac{n}{n_0}\mathbf{u}\sqrt{\beta_0} &= \pi^{-3/2}\!\int
            \frac{\cos\psi}{R^2}\,\hat e
            \left[\tfrac12 (1+s^2)e^{-S_0^2} + \tfrac{\sqrt\pi}{4}
                  \mathcal{E}(3s+2s^3)(1+\mathrm{erf}\,s)\right] dA \\
        \frac{n}{n_0}\beta_0\langle v^2\rangle &= \pi^{-3/2}\!\int
            \frac{\cos\psi}{R^2}
            \left[\left(\tfrac12 s^3+\tfrac54 s\right)e^{-S_0^2}
                  + \sqrt\pi\,\mathcal{E}
                  \left(\tfrac38+\tfrac32 s^2+\tfrac12 s^4\right)
                  (1+\mathrm{erf}\,s)\right] dA

    with :math:`\mathcal{E} = e^{-S_0^2(1-\cos^2\psi)}`, and
    :math:`T/T_0 = \tfrac23\left[\beta_0\langle v^2\rangle
    - |\mathbf{u}\sqrt{\beta_0}|^2\right]` as before.

    On the axis this reproduces (C1), (C2) and (C4) to quadrature precision,
    which is what the tests check. Off the axis it is what makes a *fair*
    comparison with a sampled DSMC field possible: a DSMC "centreline" is an
    average over a tube of finite radius, and averaging the analytical solution
    over the same cells removes a bias that is otherwise mistaken for physics --
    at a tube radius of 0.25 D it reaches 6.8% in the density at ``X/D = 1``.
    """
    x = _as_array(x)
    y = _as_array(y)
    z = _as_array(z)
    x, y, z = np.broadcast_arrays(x, y, z)
    if np.any(x <= 0.0):
        raise ValueError(
            "moments needs x > 0: in the exit plane the disk integrand has an "
            "inverse-square singularity at the field point itself.")

    s0 = float(speed_ratio)
    n_radial, n_azimuthal = _refined_order(
        nozzle_radius_m, float(np.min(x)), n_radial, n_azimuthal)
    rho, w_rho, phi, w_phi = _disk_quadrature(
        nozzle_radius_m, n_radial, n_azimuthal)

    qy = (rho[:, None] * np.cos(phi)[None, :]).ravel()
    qz = (rho[:, None] * np.sin(phi)[None, :]).ravel()
    weight = (w_rho[:, None] * w_phi[None, :]).ravel()

    shape = x.shape
    xf, yf, zf = x.ravel(), y.ravel(), z.ravel()
    density = np.empty(xf.shape, dtype=np.float64)
    momentum = np.empty(xf.shape + (3,), dtype=np.float64)
    energy = np.empty(xf.shape, dtype=np.float64)

    exp_s0 = math.exp(-s0 * s0)
    chunk = max(1, int(2_000_000 // max(1, weight.size)))
    for start in range(0, xf.size, chunk):
        stop = min(start + chunk, xf.size)
        dx = np.broadcast_to(xf[start:stop, None], (stop - start, qy.size))
        dy = yf[start:stop, None] - qy[None, :]
        dz = zf[start:stop, None] - qz[None, :]

        r2 = dx * dx + dy * dy + dz * dz
        r = np.sqrt(r2)
        cos_psi = dx / r
        s = s0 * cos_psi
        s2 = s * s
        decay = np.exp(-(s0 * s0) * (1.0 - cos_psi * cos_psi)) * (1.0 + erf(s))
        kernel = (cos_psi / r2) * weight[None, :]

        density[start:stop] = math.pi ** -1.5 * (kernel * (
            0.5 * s * exp_s0 + 0.25 * _SQRT_PI * decay * (1.0 + 2.0 * s2)
        )).sum(axis=1)

        radial = 0.5 * (1.0 + s2) * exp_s0 + 0.25 * _SQRT_PI * decay * s * (
            3.0 + 2.0 * s2)
        for axis, component in enumerate((dx, dy, dz)):
            momentum[start:stop, axis] = math.pi ** -1.5 * (
                kernel * radial * (component / r)).sum(axis=1)

        energy[start:stop] = math.pi ** -1.5 * (kernel * (
            (0.5 * s2 * s + 1.25 * s) * exp_s0
            + _SQRT_PI * decay * (0.375 + 1.5 * s2 + 0.5 * s2 * s2)
        )).sum(axis=1)

    velocity = momentum / density[:, None]
    temperature = (2.0 / 3.0) * (energy / density - (velocity ** 2).sum(axis=1))

    return {
        "n_over_n0": density.reshape(shape),
        "U_sqrt_beta0": velocity.reshape(shape + (3,)),
        "T_over_T0": temperature.reshape(shape),
    }


def describe(nozzle_radius_m: float, speed_ratio: float,
             x_over_d=(0.0, 1.0, 2.0, 5.0, 10.0)) -> list[str]:
    """Report lines: the analytical centreline at a few stations."""
    diameter = 2.0 * float(nozzle_radius_m)
    x = np.asarray(x_over_d, dtype=np.float64) * diameter
    profile = centerline_profile(x, nozzle_radius_m, speed_ratio)
    lines = [
        f"Cai collisionless centreline (S0 = {float(speed_ratio):g}, "
        f"D = {diameter:g} m)",
        f"  {'X/D':>6} {'n/n0':>12} {'U sqrt(beta0)':>15} {'T/T0':>10}",
    ]
    for i, station in enumerate(x_over_d):
        lines.append(
            f"  {station:>6.2f} {profile['n_over_n0'][i]:>12.6f} "
            f"{profile['U_sqrt_beta0'][i]:>15.6f} "
            f"{profile['T_over_T0'][i]:>10.6f}")
    return lines
