r"""Argon VHS properties, and the conversions Cai's case matrix is built from.

Three conversions define the whole study:

.. math::

    \lambda_0 = \mathrm{Kn}\, D
    \qquad
    n_0 = n_0(\lambda_0)
    \qquad
    U_0 = S_0 \sqrt{2 R T_0}

The first is Cai's own statement -- he uses the nozzle **diameter** as the
characteristic length. The second and third are worked out below.

The mean free path
------------------
For a variable-hard-sphere gas the collision cross-section depends on the
relative speed, :math:`d = d_{ref}(c_{r,ref}/c_r)^{\omega - 1/2}`. Carrying that
through the equilibrium collision rate,

.. math::

    \nu = n \langle \sigma c_r \rangle
        = 4\sqrt{\pi}\, n\, d_{ref}^2 \sqrt{k/m}\;
          T_{ref}^{\,\omega-1/2}\, T^{\,1-\omega},

and dividing the mean thermal speed :math:`\bar c = \sqrt{8kT/\pi m}` by it, the
:math:`\Gamma(5/2-\omega)` factors that appear in both cancel and what is left is

.. math::

    \lambda_{VHS} = \frac{1}{\sqrt{2}\,\pi d_{ref}^2 n}
                    \left(\frac{T}{T_{ref}}\right)^{\omega - 1/2}.
    \tag{G1}

Two things follow, and both matter here:

* at :math:`T = T_{ref}` the VHS mean free path **is** the hard-sphere one,
  :math:`\lambda = 1/(\sqrt2 \pi d^2 n)`, which is the expression
  :mod:`plumetools.markelov1999.resolution` already uses;
* away from :math:`T_{ref}` the two differ by :math:`(T/T_{ref})^{\omega-1/2}`,
  which for argon at :math:`T_0 = 300` K against :math:`T_{ref} = 273` K is
  **1.0229** -- 2.3%, small but not round-off.

Cai does not state which convention GRASP used, so the choice is recorded as an
``[ASSUMPTION]``: this package defaults to the full VHS form (G1), and the
hard-sphere form is reachable by setting ``gas.mean_free_path_convention:
hard_sphere``. Inverting (G1) gives the number density each Knudsen case needs:

.. math::

    n_0 = \frac{1}{\sqrt{2}\,\pi d_{ref}^2 \lambda_0}
          \left(\frac{T_0}{T_{ref}}\right)^{\omega - 1/2}.
    \tag{G2}

The speed ratio
---------------
Cai's :math:`S_0 = U_0/\sqrt{2RT_0}` with :math:`R = k_B/m` for the species, so
:math:`U_0 = S_0\sqrt{2RT_0}` exactly -- no gas-dynamic relation and no
:math:`\gamma` is involved. :math:`\sqrt{2RT_0}` is the most probable thermal
speed, and :math:`\beta_0 = 1/(2RT_0)` is Cai's normalising constant, so
:math:`U\sqrt{\beta_0}` is simply "the local speed ratio".
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from plumetools.cai2012.constants import AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K

#: Mean-free-path conventions this module implements.
#:
#: ``vhs``          equation (G1), temperature-corrected. The default.
#: ``hard_sphere``  ``1 / (sqrt(2) pi d^2 n)``; (G1) evaluated at ``T = T_ref``.
MEAN_FREE_PATH_CONVENTIONS = ("vhs", "hard_sphere")


@dataclass(frozen=True)
class VhsSpecies:
    """One monatomic VHS species.

    Attributes:
        name: the ``typeId`` written into ``dsmcProperties``.
        mass_kg: molecular mass [kg].
        diameter_m: VHS reference diameter [m], defined *at* ``t_ref_K``.
        omega: VHS viscosity-temperature exponent.
        t_ref_K: the reference temperature ``diameter_m`` belongs to [K].
        internal_degrees_of_freedom: 0 for a monatomic gas.

    ``diameter_m`` and ``t_ref_K`` are a **pair**. Changing one without the other
    silently changes the gas: the cross-section is :math:`\\pi d_{ref}^2` at
    :math:`T_{ref}` and scales from there, so 4.17e-10 m at 273 K and 4.17e-10 m
    at 300 K are different species.
    """

    name: str = "Ar"
    mass_kg: float = 6.63e-26
    diameter_m: float = 4.17e-10
    omega: float = 0.74
    t_ref_K: float = 273.0
    internal_degrees_of_freedom: int = 0

    def __post_init__(self) -> None:
        for field_name in ("mass_kg", "diameter_m", "t_ref_K"):
            value = getattr(self, field_name)
            if not (isinstance(value, (int, float)) and value > 0.0
                    and math.isfinite(value)):
                raise ValueError(
                    f"{self.name}: {field_name} must be finite and positive, "
                    f"got {value!r}")
        if not 0.5 <= float(self.omega) <= 1.0:
            raise ValueError(
                f"{self.name}: omega is {self.omega}; the VHS exponent lies in "
                f"[0.5, 1.0] (0.5 is hard-sphere, 1.0 is Maxwell molecules)")
        if int(self.internal_degrees_of_freedom) < 0:
            raise ValueError(
                f"{self.name}: internal_degrees_of_freedom must be >= 0")

    @property
    def gas_constant_j_per_kg_k(self) -> float:
        """Specific gas constant ``R = k_B / m`` [J/(kg K)].

        This is the ``R`` in Cai's :math:`S_0 = U_0/\\sqrt{2RT_0}`, per unit
        mass of the species -- not the molar gas constant.
        """
        return BOLTZMANN_J_PER_K / self.mass_kg

    @property
    def gamma(self) -> float:
        """Ratio of specific heats, from the degrees of freedom.

        Reported, not configured: a monatomic gas with no internal modes has
        ``gamma = 5/3`` and there is no second place to state it. Nothing in the
        Cai solution uses it -- the collisionless solution is built from the
        exit distribution, not from an isentropic relation -- but
        ``dsmcProperties`` states the DoF count and this is the same fact.
        """
        dof = 3 + int(self.internal_degrees_of_freedom)
        return (dof + 2) / dof

    @property
    def molar_mass_g_per_mol(self) -> float:
        """Molar mass [g/mol] implied by ``mass_kg``. Reported, never an input."""
        return self.mass_kg * AVOGADRO_PER_MOL * 1000.0

    def most_probable_speed(self, temperature_K: float) -> float:
        """``sqrt(2 R T)`` [m/s] -- Cai's velocity normalisation.

        Also ``1/sqrt(beta0)`` when evaluated at ``T0``.
        """
        return math.sqrt(2.0 * self.gas_constant_j_per_kg_k * _positive(
            temperature_K, "temperature_K"))

    def mean_thermal_speed(self, temperature_K: float) -> float:
        """``sqrt(8 R T / pi)`` [m/s], the Maxwellian mean speed."""
        return math.sqrt(
            8.0 * self.gas_constant_j_per_kg_k
            * _positive(temperature_K, "temperature_K") / math.pi)

    def beta(self, temperature_K: float) -> float:
        """Cai's ``beta = 1 / (2 R T)`` [s^2/m^2].

        An inverse squared speed, so ``sqrt(beta)`` times a velocity is
        dimensionless -- which is why Cai plots ``U1 sqrt(beta0)``.
        """
        return 1.0 / (2.0 * self.gas_constant_j_per_kg_k
                      * _positive(temperature_K, "temperature_K"))


#: Argon as this repository already models it.
#:
#: ``mass`` and ``diameter`` are the values ``cases/3d-inflow`` runs with;
#: ``omega`` is that case's 0.74. Bird's tabulated argon is ``omega = 0.81`` at
#: the same 4.17e-10 m / 273 K reference, so 0.74 is an **existing repository
#: value, not a paper value or a literature value** -- see
#: ``docs/cai2012-case.md``. It is kept because §3 of the case specification asks
#: to start from the repository's argon rather than introduce new molecular
#: parameters, and because the quantity it affects here (the Kn -> n0 conversion,
#: through ``(T0/Tref)**(omega-1/2)``) moves by 0.6% between the two.
ARGON = VhsSpecies(
    name="Ar",
    mass_kg=6.63e-26,
    diameter_m=4.17e-10,
    omega=0.74,
    t_ref_K=273.0,
    internal_degrees_of_freedom=0,
)


def _positive(value: float, what: str) -> float:
    value = float(value)
    if not (value > 0.0 and math.isfinite(value)):
        raise ValueError(f"{what} must be finite and positive, got {value}")
    return value


def vhs_temperature_factor(species: VhsSpecies, temperature_K: float) -> float:
    """``(T / T_ref)**(omega - 1/2)`` -- the VHS correction in (G1).

    1.0 exactly when ``T == T_ref``, which is why the hard-sphere expression is
    not a different model but the same one evaluated at the reference state.
    """
    return (_positive(temperature_K, "temperature_K") / species.t_ref_K) ** (
        float(species.omega) - 0.5)


def mean_free_path(number_density_per_m3: float, species: VhsSpecies,
                   temperature_K: float, *, convention: str = "vhs") -> float:
    """Equilibrium mean free path [m], equation (G1).

    Args:
        number_density_per_m3: number density [1/m^3].
        species: the VHS species.
        temperature_K: gas temperature [K]. Ignored by ``hard_sphere``.
        convention: ``"vhs"`` or ``"hard_sphere"``.

    Returns:
        The mean free path [m]. ``inf`` for a zero density, which is the correct
        limit and keeps a vacuum cell from raising in a diagnostic loop.

    Raises:
        ValueError: for a negative density or an unknown convention.
    """
    if convention not in MEAN_FREE_PATH_CONVENTIONS:
        raise ValueError(
            f"unknown mean-free-path convention {convention!r}; "
            f"valid: {list(MEAN_FREE_PATH_CONVENTIONS)}")
    n = float(number_density_per_m3)
    if n < 0.0 or not math.isfinite(n):
        raise ValueError(f"number density must be finite and >= 0, got {n}")
    if n == 0.0:
        return float("inf")

    hard_sphere = 1.0 / (math.sqrt(2.0) * math.pi * species.diameter_m ** 2 * n)
    if convention == "hard_sphere":
        return hard_sphere
    return hard_sphere * vhs_temperature_factor(species, temperature_K)


def number_density_from_mean_free_path(mean_free_path_m: float,
                                       species: VhsSpecies,
                                       temperature_K: float, *,
                                       convention: str = "vhs") -> float:
    """Number density [1/m^3] from a mean free path -- equation (G2).

    The exact inverse of :func:`mean_free_path`; the round trip is tested.
    """
    lam = _positive(mean_free_path_m, "mean_free_path_m")
    if convention not in MEAN_FREE_PATH_CONVENTIONS:
        raise ValueError(
            f"unknown mean-free-path convention {convention!r}; "
            f"valid: {list(MEAN_FREE_PATH_CONVENTIONS)}")
    hard_sphere = 1.0 / (math.sqrt(2.0) * math.pi * species.diameter_m ** 2 * lam)
    if convention == "hard_sphere":
        return hard_sphere
    return hard_sphere * vhs_temperature_factor(species, temperature_K)


def mean_free_path_from_kn(knudsen: float, characteristic_length_m: float) -> float:
    """``lambda0 = Kn * L`` [m].

    ``L`` is the nozzle **diameter** for this study, which Cai states explicitly.
    Using the radius instead would halve every mean free path and so double every
    number density -- a factor of two hiding in a definition, which is why the
    characteristic length is a named parameter and not a constant.
    """
    return _positive(knudsen, "knudsen") * _positive(
        characteristic_length_m, "characteristic_length_m")


def number_density_from_kn(knudsen: float, characteristic_length_m: float,
                           species: VhsSpecies, temperature_K: float, *,
                           convention: str = "vhs") -> float:
    """Exit number density [1/m^3] for one Knudsen case.

    ``Kn -> lambda0 = Kn L -> n0`` by (G2). This is the function the whole case
    matrix rests on: no case anywhere in this family sets a number density by
    hand, so a density can only ever be wrong by being derived from the wrong Kn.
    """
    lam = mean_free_path_from_kn(knudsen, characteristic_length_m)
    return number_density_from_mean_free_path(
        lam, species, temperature_K, convention=convention)


def knudsen_from_number_density(number_density_per_m3: float,
                                characteristic_length_m: float,
                                species: VhsSpecies, temperature_K: float, *,
                                convention: str = "vhs") -> float:
    """The Knudsen number a given density implies. The inverse audit path."""
    lam = mean_free_path(number_density_per_m3, species, temperature_K,
                         convention=convention)
    return lam / _positive(characteristic_length_m, "characteristic_length_m")


def speed_from_speed_ratio(speed_ratio: float, species: VhsSpecies,
                           temperature_K: float) -> float:
    """``U0 = S0 sqrt(2 R T0)`` [m/s].

    Cai's definition exactly. A negative speed ratio is accepted -- it means a
    drift into the domain from the other side -- but zero is not rejected either,
    since ``S0 = 0`` is the pure-effusion limit the analytical solution is checked
    against in the unit tests.
    """
    return float(speed_ratio) * species.most_probable_speed(temperature_K)


def speed_ratio_from_speed(speed_m_per_s: float, species: VhsSpecies,
                           temperature_K: float) -> float:
    """``S = U / sqrt(2 R T)``. The inverse of :func:`speed_from_speed_ratio`."""
    return float(speed_m_per_s) / species.most_probable_speed(temperature_K)


def maxwellian_number_flux(number_density_per_m3: float, speed_ratio: float,
                           species: VhsSpecies, temperature_K: float) -> float:
    """One-way number flux through a plane, for a drifting Maxwellian [1/(m^2 s)].

    Bird eq. 4.22, with ``s`` the speed ratio **normal to the plane**:

    .. math::

        \\Gamma = \\frac{n}{2\\sqrt{\\pi}\\,\\beta}
                  \\left[e^{-s^2} + \\sqrt{\\pi}\\, s\\,(1+\\mathrm{erf}\\,s)\\right],
        \\qquad \\beta = 1/\\sqrt{2RT}.

    This is what ``plumeFieldInflow`` accumulates per face per timestep, so it is
    the theoretical value the OpenFOAM flux test compares the measured injection
    rate against. At ``s = 0`` it reduces to the effusion flux ``n c_bar / 4``,
    which is asserted in the unit tests.
    """
    n = float(number_density_per_m3)
    if n < 0.0:
        raise ValueError(f"number density must be >= 0, got {n}")
    s = float(speed_ratio)
    most_probable = species.most_probable_speed(temperature_K)
    return (n * most_probable / (2.0 * math.sqrt(math.pi))) * (
        math.exp(-s * s) + math.sqrt(math.pi) * s * (1.0 + math.erf(s)))


def collision_time(number_density_per_m3: float, species: VhsSpecies,
                   temperature_K: float, *, convention: str = "vhs") -> float:
    """Mean collision time [s], ``lambda / c_bar``.

    **[ASSUMPTION].** Cai states a reference collision time ``t0`` and requires
    ``dt/t0 = 1``, but does not define ``t0``. This is the most common reading --
    the time to travel one mean free path at the mean thermal speed, which is
    also ``1/nu`` for the collision rate ``nu = c_bar/lambda`` that (G1) was
    derived from. Any other convention differs by an O(1) factor, so the value is
    reported alongside the time step actually used rather than being adopted as
    if it were a paper number. See ``docs/cai2012-case.md``.
    """
    lam = mean_free_path(number_density_per_m3, species, temperature_K,
                         convention=convention)
    return lam / species.mean_thermal_speed(temperature_K)


def describe(species: VhsSpecies, temperature_K: float) -> list[str]:
    """Report lines for the gas model, printed before a case is meshed or run."""
    return [
        f"  species                  {species.name}  "
        f"(m = {species.mass_kg:.4e} kg, M = {species.molar_mass_g_per_mol:.4f} g/mol)",
        f"  VHS                      d_ref = {species.diameter_m:.4e} m, "
        f"omega = {species.omega:g}, T_ref = {species.t_ref_K:g} K",
        f"  R = k_B/m                {species.gas_constant_j_per_kg_k:.4f} J/(kg K)",
        f"  sqrt(2 R T0)             {species.most_probable_speed(temperature_K):.3f} m/s"
        f"   (T0 = {float(temperature_K):g} K)",
        f"  beta0 = 1/(2 R T0)       {species.beta(temperature_K):.6e} s^2/m^2",
        f"  VHS factor (T0/T_ref)^(omega-1/2)   "
        f"{vhs_temperature_factor(species, temperature_K):.6f}",
    ]
