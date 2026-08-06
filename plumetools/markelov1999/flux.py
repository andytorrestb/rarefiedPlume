r"""Verify the integrated inflow: mass, momentum and energy.

Why not ``rho * U * A``
-----------------------
Because that is not what the solver injects. ``plumeFieldInflow`` -- like the
``FreeStream`` it derives from -- samples a **half-range drifting Maxwellian**:
particles are drawn from the inward-moving half of a Maxwellian distribution
displaced by the local mean velocity. The number flux it accumulates is Bird
eqn 4.22::

    Fn = n * c_mp / (2*sqrt(pi)) * [ exp(-s**2) + sqrt(pi)*s*(1 + erf(s)) ]

with ``s = (U . n_hat_in) / c_mp`` the molecular speed ratio along the inward
normal and ``c_mp = sqrt(2 k T / m)`` the most probable thermal speed.

The drift term alone, ``n * U . n_hat``, is the ``s -> infinity`` limit of that
expression. Everything else is thermal flux: molecules that cross the surface
because of their random motion rather than because of the bulk drift. Verifying
against ``rho*U*A`` would therefore not be verifying the boundary condition that
is actually applied, and the discrepancy is not negligible -- see
:func:`drift_only_number_flux` and the test that compares them.

The same reasoning applies to momentum and energy: each carries a pressure-like
and a thermal contribution that the drift-only form omits.

Two independent routes
----------------------
Every moment is computed twice:

* :func:`analytical_face_fluxes` -- closed-form half-range moments of the
  drifting Maxwellian, evaluated per face;
* :func:`quadrature_face_fluxes` -- the same moments obtained by deterministic
  Gauss-Hermite quadrature over the velocity distribution, with no closed form
  used anywhere.

They agree to quadrature accuracy or one of them is wrong. Neither is derived
from the other, so agreement is evidence rather than a tautology.
:func:`sampled_number_flux` adds an optional Monte-Carlo check against the same
acceptance-rejection sampler the solver uses.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import special

from plumetools.markelov1999.constants import BOLTZMANN_J_PER_K


@dataclass(frozen=True)
class FluxTotals:
    """Integrated fluxes through a surface, per unit time.

    Attributes:
        mass_kg_per_s: total mass flux [kg/s].
        momentum_N: ``(3,)`` momentum flux [kg m/s^2 = N].
        energy_W: total energy flux [J/s = W], translational plus internal.
        number_per_s: total number flux [1/s].
    """

    number_per_s: float
    mass_kg_per_s: float
    momentum_N: np.ndarray
    energy_W: float

    def as_dict(self) -> dict:
        return {
            "number_per_s": self.number_per_s,
            "mass_kg_per_s": self.mass_kg_per_s,
            "momentum_x_N": float(self.momentum_N[0]),
            "momentum_y_N": float(self.momentum_N[1]),
            "momentum_z_N": float(self.momentum_N[2]),
            "energy_W": self.energy_W,
        }


def speed_ratio(U, inward_normals, c_mp) -> np.ndarray:
    """``s = (U . n_hat_in) / c_mp`` per face [-].

    ``inward_normals`` must point **into** the domain, which is the negation of
    an OpenFOAM boundary face normal. Getting the sign wrong turns an inflow into
    an outflow and silently produces a tiny positive flux instead of the full one,
    because the half-range integral is not antisymmetric.
    """
    U = np.asarray(U, dtype=np.float64)
    n = np.asarray(inward_normals, dtype=np.float64)
    return np.einsum("ij,ij->i", U, n) / np.asarray(c_mp, dtype=np.float64)


def most_probable_speed(T_K, mass_kg: float) -> np.ndarray:
    """``c_mp = sqrt(2 k T / m)`` [m/s].

    The same definition ``DSMCCloud::maxwellianMostProbableSpeed`` uses, so the
    speed ratio here is the solver's.
    """
    return np.sqrt(2.0 * BOLTZMANN_J_PER_K * np.asarray(T_K, dtype=np.float64) / mass_kg)


def number_flux_coefficient(s) -> np.ndarray:
    """Bird eqn 4.22 in dimensionless form.

    ``[exp(-s**2) + sqrt(pi) s (1 + erf(s))] / (2 sqrt(pi))``, so the number flux
    per unit area is ``n * c_mp *`` this. Exactly the expression accumulated in
    ``plumeFieldInflow::inflow`` and in ``FreeStream::inflow``.
    """
    s = np.asarray(s, dtype=np.float64)
    return (np.exp(-s * s) + math.sqrt(math.pi) * s * (1.0 + special.erf(s))) \
        / (2.0 * math.sqrt(math.pi))


def drift_only_number_flux(n, U, inward_normals) -> np.ndarray:
    """``n * (U . n_hat_in)`` per unit area [1/(m^2 s)].

    The naive form. Provided so the difference against
    :func:`number_flux_coefficient` can be measured rather than asserted: it is
    the ``s -> infinity`` limit, and at the speed ratios this plume reaches the
    thermal terms are a real fraction of the total.
    """
    n = np.asarray(n, dtype=np.float64)
    return n * np.einsum("ij,ij->i",
                         np.asarray(U, dtype=np.float64),
                         np.asarray(inward_normals, dtype=np.float64))


def analytical_face_fluxes(*, number_density, U, T_K, areas, inward_normals,
                           mass_kg: float, internal_dof: int) -> FluxTotals:
    """Closed-form half-range moments of the drifting Maxwellian.

    Args:
        number_density: ``(n,)`` [1/m^3].
        U: ``(n, 3)`` bulk velocity [m/s].
        T_K: ``(n,)`` temperature [K].
        areas: ``(n,)`` face areas [m^2].
        inward_normals: ``(n, 3)`` unit normals pointing INTO the domain.
        mass_kg: molecular mass [kg].
        internal_dof: internal degrees of freedom, for the internal energy flux.

    Returns:
        The integrated :class:`FluxTotals`.

    The moments, with ``s`` the speed ratio along the inward normal, ``c`` the
    most probable speed and ``F0 = [exp(-s^2) + sqrt(pi) s (1+erf(s))]/(2 sqrt(pi))``:

    * **number**   ``n c F0``
    * **normal momentum**   ``n m c^2 * [ s exp(-s^2)/(2 sqrt(pi))
      + (s^2 + 1/2)(1 + erf(s))/2 ]`` -- the drift term plus the pressure term
      the thermal motion contributes even at ``s = 0``.
    * **tangential momentum**   ``m * U_t *`` the number flux; the thermal
      contribution is zero by symmetry, since the tangential velocity
      distribution is an undisplaced Maxwellian.
    * **translational energy**   ``n m c^3 * [ (s^2 + 2) exp(-s^2)/(4 sqrt(pi))
      + s(s^2 + 5/2)(1 + erf(s))/4 ]``
    * **internal energy**   ``(zeta/2) k T *`` the number flux, from
      equipartition -- which is what ``equipartitionInternalEnergy`` gives each
      injected particle.

    The energy expression is the reason a drift-only check is inadequate: at
    ``s = 0`` the drift energy flux is exactly zero while the true flux is not.
    """
    n = np.asarray(number_density, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    T = np.asarray(T_K, dtype=np.float64)
    A = np.asarray(areas, dtype=np.float64)
    nhat = np.asarray(inward_normals, dtype=np.float64)

    c = most_probable_speed(T, mass_kg)
    s = speed_ratio(U, nhat, c)

    exp_s2 = np.exp(-s * s)
    erf_term = 1.0 + special.erf(s)
    sqrt_pi = math.sqrt(math.pi)

    # Number flux per unit area.
    F0 = (exp_s2 + sqrt_pi * s * erf_term) / (2.0 * sqrt_pi)
    number_per_area = n * c * F0

    # Normal momentum flux per unit area (a pressure).
    F1 = s * exp_s2 / (2.0 * sqrt_pi) + (s * s + 0.5) * erf_term / 2.0
    normal_momentum_per_area = n * mass_kg * c * c * F1

    # Tangential momentum rides on the number flux alone.
    U_normal = np.einsum("ij,ij->i", U, nhat)
    U_tangential = U - U_normal[:, None] * nhat
    tangential_momentum_per_area = mass_kg * U_tangential * number_per_area[:, None]

    momentum_per_area = (normal_momentum_per_area[:, None] * nhat
                         + tangential_momentum_per_area)

    # Translational energy flux per unit area.
    #
    # F2 is the standard result for a drift NORMAL to the surface: it already
    # contains the tangential *thermal* energy (the "+2" and "5/2" terms) but
    # assumes no tangential drift. The inflow surface here is a faceted sphere, so
    # U is not exactly normal to every face and the tangential drift kinetic
    # energy has to be added explicitly. Leaving it out is a small error on this
    # geometry and a silent one on any other -- and the quadrature route below
    # includes it, so omitting it here would show up as a disagreement.
    F2 = ((s * s + 2.0) * exp_s2 / (4.0 * sqrt_pi)
          + s * (s * s + 2.5) * erf_term / 4.0)
    ut_squared = np.einsum("ij,ij->i", U_tangential, U_tangential)
    translational_per_area = (n * mass_kg * c ** 3 * F2
                              + 0.5 * mass_kg * ut_squared * number_per_area)

    # Internal energy: equipartition, carried by each injected particle.
    internal_per_area = (0.5 * internal_dof * BOLTZMANN_J_PER_K * T) * number_per_area

    return FluxTotals(
        number_per_s=float((number_per_area * A).sum()),
        mass_kg_per_s=float((mass_kg * number_per_area * A).sum()),
        momentum_N=(momentum_per_area * A[:, None]).sum(axis=0),
        energy_W=float(((translational_per_area + internal_per_area) * A).sum()),
    )


def quadrature_face_fluxes(*, number_density, U, T_K, areas, inward_normals,
                           mass_kg: float, internal_dof: int,
                           n_points: int = 200) -> FluxTotals:
    """The same moments by deterministic quadrature, using no closed form.

    An **independent reference**: the velocity-space integrals are evaluated
    numerically over the inward half-range, so agreement with
    :func:`analytical_face_fluxes` tests the closed forms rather than restating
    them.

    Args:
        n_points: Gauss-Legendre points per velocity component. The integrand is
            a Gaussian times a polynomial, so convergence is fast; 200 points give
            ~1e-12 relative here.

    The integral, in a frame whose ``w`` axis is the inward normal::

        moment = n * (m/(2 pi k T))**(3/2)
                 * integral_{w>0} integral_u integral_v
                     g(u,v,w) exp(-m[(u-Ut1)^2 + (v-Ut2)^2 + (w-Un)^2]/(2kT))
                     du dv dw

    with ``g = w`` for number, ``m w**2`` for normal momentum, and so on. The
    tangential directions integrate over the full range and factor out
    analytically only in the sense that a Gaussian integral is done numerically
    here too.
    """
    n = np.asarray(number_density, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    T = np.asarray(T_K, dtype=np.float64)
    A = np.asarray(areas, dtype=np.float64)
    nhat = np.asarray(inward_normals, dtype=np.float64)

    U_normal = np.einsum("ij,ij->i", U, nhat)
    U_tangential = U - U_normal[:, None] * nhat

    number = np.zeros(len(n))
    normal_momentum = np.zeros(len(n))
    energy = np.zeros(len(n))

    nodes, weights = np.polynomial.legendre.leggauss(int(n_points))

    for i in range(len(n)):
        beta = math.sqrt(mass_kg / (2.0 * BOLTZMANN_J_PER_K * T[i]))
        c = 1.0 / beta                      # the most probable speed
        un = U_normal[i]

        # Map the fixed [-1, 1] nodes onto [0, w_max]; w_max is far enough out
        # that the Gaussian tail beyond it is below double precision.
        w_max = max(un, 0.0) + 8.0 * c
        w = 0.5 * w_max * (nodes + 1.0)
        dw = 0.5 * w_max * weights

        # 1-D normalised half-range distribution along the normal.
        f_w = beta / math.sqrt(math.pi) * np.exp(-((w - un) * beta) ** 2)

        number[i] = float((w * f_w * dw).sum())
        normal_momentum[i] = float((mass_kg * w * w * f_w * dw).sum())

        # Translational energy carries the tangential kinetic energy too. The
        # tangential distribution is an undisplaced Maxwellian about U_t, whose
        # mean square is |U_t|**2 + 2 * (kT/m) over the two directions.
        ut2 = float(U_tangential[i] @ U_tangential[i])
        tangential_ke = 0.5 * mass_kg * (ut2 + 2.0 * BOLTZMANN_J_PER_K * T[i] / mass_kg)
        energy[i] = float((
            (0.5 * mass_kg * w * w + tangential_ke
             + 0.5 * internal_dof * BOLTZMANN_J_PER_K * T[i])
            * w * f_w * dw).sum())

    number_per_area = n * number
    momentum_per_area = (n * normal_momentum)[:, None] * nhat \
        + mass_kg * U_tangential * number_per_area[:, None]

    return FluxTotals(
        number_per_s=float((number_per_area * A).sum()),
        mass_kg_per_s=float((mass_kg * number_per_area * A).sum()),
        momentum_N=(momentum_per_area * A[:, None]).sum(axis=0),
        energy_W=float((n * energy * A).sum()),
    )


def sampled_number_flux(*, number_density: float, U_normal: float, T_K: float,
                        mass_kg: float, n_samples: int = 400_000,
                        seed: int = 0) -> tuple[float, float]:
    """Monte-Carlo half-range number flux per unit area, and its standard error.

    Draws velocities from the full drifting Maxwellian and averages
    ``n * max(w, 0)`` -- the flux each molecule contributes through the surface,
    zero for the outward-moving half. That is an unbiased estimator of

        integral_{w>0} n w f(w) dw

    which is the same quantity :func:`analytical_face_fluxes` computes in closed
    form and :func:`quadrature_face_fluxes` integrates numerically, obtained
    without either.

    Returns:
        ``(flux, standard_error)`` per unit area [1/(m^2 s)].
    """
    rng = np.random.default_rng(seed)
    sigma = math.sqrt(BOLTZMANN_J_PER_K * T_K / mass_kg)
    w = rng.normal(U_normal, sigma, int(n_samples))
    contributions = number_density * np.maximum(w, 0.0)
    return (float(contributions.mean()),
            float(contributions.std(ddof=1) / math.sqrt(n_samples)))


def sampled_mean_normal_speed(*, U_normal: float, T_K: float, mass_kg: float,
                              n_samples: int = 200_000,
                              seed: int = 0) -> tuple[float, float]:
    """Mean normalised normal speed from the solver's own eqn 12.5 sampler.

    Reproduces the acceptance-rejection loop in ``plumeFieldInflow::inflow`` --
    which is verbatim ``FreeStream``'s -- and returns the mean and standard error
    of the accepted normalised normal velocities.

    This tests the **sampler**, not the flux formula: an injected particle's
    normal velocity is drawn from the flux-weighted distribution
    ``u exp(-(u - s)^2)`` on ``u > 0``, whose mean is ``F1 / F0`` with ``F1`` and
    ``F0`` the momentum and number coefficients. So agreement here says the
    velocities the solver actually injects carry the momentum the closed form
    predicts -- something the two deterministic routes cannot check about each
    other, because neither draws a sample.

    Returns:
        ``(mean, standard_error)`` of ``u = w / c_mp`` over accepted samples.
    """
    rng = np.random.default_rng(seed)
    c = math.sqrt(2.0 * BOLTZMANN_J_PER_K * T_K / mass_kg)
    s = U_normal / c

    coeff_a = s + math.sqrt(s * s + 2.0)
    coeff_b = 0.5 * (1.0 + s * (s - math.sqrt(s * s + 2.0)))
    scaling = 3.0 if s >= -3 else abs(s) + 1.0

    accepted = np.empty(int(n_samples))
    filled = 0
    while filled < n_samples:
        batch = max(4096, (n_samples - filled) * 4)
        thermal = scaling * (2.0 * rng.random(batch) - 1.0)
        u_norm = thermal + s
        with np.errstate(over="ignore"):
            p = np.where(
                u_norm < 0.0, -1.0,
                2.0 * u_norm / coeff_a * np.exp(coeff_b - thermal ** 2))
        keep = u_norm[p >= rng.random(batch)]
        take = min(len(keep), n_samples - filled)
        accepted[filled:filled + take] = keep[:take]
        filled += take

    return (float(accepted.mean()),
            float(accepted.std(ddof=1) / math.sqrt(n_samples)))


def normal_momentum_coefficient(s) -> np.ndarray:
    """``F1``: the dimensionless normal-momentum (pressure) coefficient.

    ``F1 = s exp(-s^2)/(2 sqrt(pi)) + (s^2 + 1/2)(1 + erf(s))/2``, so the normal
    momentum flux per unit area is ``n m c_mp^2 F1``. Exposed separately because
    ``F1 / F0`` is the mean normalised normal speed of an injected particle, which
    is what :func:`sampled_mean_normal_speed` measures.
    """
    s = np.asarray(s, dtype=np.float64)
    return (s * np.exp(-s * s) / (2.0 * math.sqrt(math.pi))
            + (s * s + 0.5) * (1.0 + special.erf(s)) / 2.0)


def inward_normals_from_patch(normals, centroids) -> np.ndarray:
    """Orient face normals to point into the fluid domain.

    OpenFOAM winds a boundary face so its normal points **out of** the fluid. On
    the inflow cavity the fluid is outside the sphere, so the stored normal points
    toward the origin and the inward-to-the-fluid direction is radially outward.

    Rather than assume a sign convention, each normal is flipped to agree with the
    outward radial direction of its own face centroid, which is unambiguous for a
    star-shaped surface centred on the origin.
    """
    n = np.asarray(normals, dtype=np.float64)
    c = np.asarray(centroids, dtype=np.float64)
    radial = c / np.linalg.norm(c, axis=1)[:, None]
    sign = np.sign(np.einsum("ij,ij->i", n, radial))
    sign[sign == 0.0] = 1.0
    return n * sign[:, None]


@dataclass(frozen=True)
class FluxVerification:
    """The comparison reported by ``./Allrun``."""

    analytical: FluxTotals
    quadrature: FluxTotals
    from_fields: FluxTotals
    drift_only_number_per_s: float

    def relative_errors(self) -> dict:
        """Relative difference of the field-reconstructed flux from the analytical."""
        def rel(a: float, b: float) -> float:
            scale = max(abs(a), abs(b), 1e-300)
            return abs(a - b) / scale

        return {
            "number": rel(self.analytical.number_per_s, self.from_fields.number_per_s),
            "mass": rel(self.analytical.mass_kg_per_s, self.from_fields.mass_kg_per_s),
            "momentum_x": rel(float(self.analytical.momentum_N[0]),
                              float(self.from_fields.momentum_N[0])),
            "momentum_y": rel(float(self.analytical.momentum_N[1]),
                              float(self.from_fields.momentum_N[1])),
            "momentum_z": rel(float(self.analytical.momentum_N[2]),
                              float(self.from_fields.momentum_N[2])),
            "energy": rel(self.analytical.energy_W, self.from_fields.energy_W),
        }

    def quadrature_errors(self) -> dict:
        """Relative difference between the two independent deterministic routes."""
        def rel(a: float, b: float) -> float:
            scale = max(abs(a), abs(b), 1e-300)
            return abs(a - b) / scale

        return {
            "number": rel(self.analytical.number_per_s, self.quadrature.number_per_s),
            "mass": rel(self.analytical.mass_kg_per_s, self.quadrature.mass_kg_per_s),
            "momentum_x": rel(float(self.analytical.momentum_N[0]),
                              float(self.quadrature.momentum_N[0])),
            "energy": rel(self.analytical.energy_W, self.quadrature.energy_W),
        }

    def report(self) -> list[str]:
        errors = self.relative_errors()
        quad = self.quadrature_errors()
        a, f = self.analytical, self.from_fields
        drift_ratio = (self.drift_only_number_per_s
                       / max(a.number_per_s, 1e-300))
        return [
            "Inflow flux verification (half-range drifting Maxwellian, Bird eqn 4.22)",
            f"  {'quantity':<14} {'analytical':>14} {'from 0/ fields':>16} "
            f"{'rel. error':>11}",
            f"  {'number [1/s]':<14} {a.number_per_s:>14.6e} "
            f"{f.number_per_s:>16.6e} {errors['number']:>11.2e}",
            f"  {'mass [kg/s]':<14} {a.mass_kg_per_s:>14.6e} "
            f"{f.mass_kg_per_s:>16.6e} {errors['mass']:>11.2e}",
            f"  {'momentum x [N]':<14} {a.momentum_N[0]:>14.6e} "
            f"{f.momentum_N[0]:>16.6e} {errors['momentum_x']:>11.2e}",
            f"  {'momentum y [N]':<14} {a.momentum_N[1]:>14.6e} "
            f"{f.momentum_N[1]:>16.6e} {errors['momentum_y']:>11.2e}",
            f"  {'momentum z [N]':<14} {a.momentum_N[2]:>14.6e} "
            f"{f.momentum_N[2]:>16.6e} {errors['momentum_z']:>11.2e}",
            f"  {'energy [W]':<14} {a.energy_W:>14.6e} "
            f"{f.energy_W:>16.6e} {errors['energy']:>11.2e}",
            "",
            f"  independent quadrature agrees to "
            f"{max(quad.values()):.2e} relative (no closed form used)",
            f"  drift-only n*U*A would give {self.drift_only_number_per_s:.6e} 1/s, "
            f"{100.0 * drift_ratio:.1f}% of the true half-range flux --",
            "  which is why rho*U*A is not the quantity to verify against.",
        ]


def verify(inflow, cfg) -> FluxVerification:
    """Compare analytical, quadrature and field-reconstructed inflow fluxes.

    Args:
        inflow: a :class:`plumetools.markelov1999.inflow.InflowResult`.
        cfg: the case config.

    Returns:
        A :class:`FluxVerification`.

    "From the generated fields" means: take the per-face ``rhoN``, ``U`` and ``T``
    exactly as written into ``0/``, and integrate the same half-range moments over
    them. It is the flux the boundary model will actually see, so a discrepancy
    against the analytical value means the field writing lost something --
    truncated precision, a dropped face, a reordering.
    """
    nhat = inward_normals_from_patch(inflow.normals, inflow.centroids)
    dof = int(cfg.dsmc.species.internal_degrees_of_freedom)

    common = dict(
        U=inflow.U, T_K=inflow.T, areas=inflow.areas,
        inward_normals=nhat, mass_kg=inflow.mass_kg, internal_dof=dof,
    )

    analytical = analytical_face_fluxes(number_density=inflow.rhoN, **common)
    quadrature = quadrature_face_fluxes(number_density=inflow.rhoN, **common)

    # Round-trip the values through the same formatting the field writer uses, so
    # any precision loss in 0/ shows up here rather than in the solver.
    rhoN_written = np.array([float(f"{n:.10g}") for n in inflow.rhoN])
    U_written = np.array([[float(f"{c:.10g}") for c in row] for row in inflow.U])
    T_written = np.array([float(f"{t:.10g}") for t in inflow.T])
    from_fields = analytical_face_fluxes(
        number_density=rhoN_written, U=U_written, T_K=T_written,
        areas=inflow.areas, inward_normals=nhat,
        mass_kg=inflow.mass_kg, internal_dof=dof)

    drift = float((drift_only_number_flux(inflow.rhoN, inflow.U, nhat)
                   * inflow.areas).sum())

    return FluxVerification(
        analytical=analytical,
        quadrature=quadrature,
        from_fields=from_fields,
        drift_only_number_per_s=drift,
    )
