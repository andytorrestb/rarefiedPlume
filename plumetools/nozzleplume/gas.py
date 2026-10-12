r"""One VHS species, and the kinetic-theory relations an exit state is built from.

Generalised from :mod:`plumetools.cai2012.gas`, which is argon and monatomic by
construction. Nothing here assumes a species: the mass, the diameter, the
viscosity exponent and the internal degrees of freedom are all inputs.

The mean free path
------------------
For a variable-hard-sphere gas the equilibrium mean free path is

.. math::

    \lambda = \frac{1}{\sqrt{2}\,\pi d_{ref}^2\, n}
              \left(\frac{T}{T_{ref}}\right)^{\omega - 1/2},
    \tag{G1}

derived in :mod:`plumetools.cai2012.gas`. At :math:`\omega = 1/2` the correction
is identically one, the reference temperature drops out, and (G1) is the
hard-sphere mean free path with collision cross-section
:math:`\sigma = \pi d^2` -- which is the model a data sheet is describing when it
quotes a single "molecular diameter" and a single "collision cross-section" and
no reference temperature. That is why ``omega`` defaults to 0.5 here rather than
to a species' tabulated value: with it, ``dsmcFoam``'s ``VariableHardSphere``
collides molecules at exactly the cross-section quoted.

Fluxes through a plane
----------------------
A drifting Maxwellian at speed ratio :math:`s = u_n/\sqrt{2RT}` normal to a
plane carries, in the molecules crossing it **one way**,

.. math::

    \Gamma &= \frac{n c_m}{2\sqrt\pi}
              \left[e^{-s^2} + \sqrt\pi\, s\,(1+\mathrm{erf}\,s)\right] \\
    p_i    &= \frac{\rho c_m^2}{2\sqrt\pi}
              \left[s\,e^{-s^2}
                    + \sqrt\pi\left(\tfrac12 + s^2\right)(1+\mathrm{erf}\,s)\right] \\
    q_i    &= \frac{\rho c_m^3}{4\sqrt\pi}
              \left[(s^2+2)\,e^{-s^2}
                    + \sqrt\pi\, s\left(s^2+\tfrac52\right)(1+\mathrm{erf}\,s)\right]

with :math:`c_m = \sqrt{2RT}` (Bird 1994, eqs. 4.22, 4.25, 4.27). The first is
what ``plumeFieldInflow`` accumulates per face per step, so it is the injection
rate. The second, over the exit area, is the momentum the jet carries -- its
thrust, and the load a plate that stopped the whole jet would feel before any
re-emission. The third is the translational power. All three reduce to the
continuum values :math:`nU`, :math:`\rho(U^2 + RT)` and
:math:`\rho U(U^2/2 + \tfrac52 RT)` once :math:`s \gg 1`, and to the effusion
values at :math:`s = 0`; both limits are asserted in the unit tests.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

from plumetools.nozzleplume.constants import AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K

#: Mean-free-path conventions this module implements.
#:
#: ``vhs``          equation (G1), temperature-corrected. The default.
#: ``hard_sphere``  ``1 / (sqrt(2) pi d^2 n)``; (G1) evaluated at ``T = T_ref``.
#:
#: The two are the same number when ``omega = 0.5``.
MEAN_FREE_PATH_CONVENTIONS = ("vhs", "hard_sphere")

#: What a species name may look like. It becomes an OpenFOAM ``word`` -- the
#: ``typeId`` -- and part of a file name, ``0/boundaryNumberDensity_<name>``, so
#: it is held to an identifier rather than to everything a ``word`` tolerates.
SPECIES_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


@dataclass(frozen=True)
class VhsSpecies:
    """One VHS species.

    Attributes:
        name: the ``typeId`` written into ``dsmcProperties``.
        mass_kg: molecular mass [kg].
        diameter_m: VHS reference diameter [m], defined *at* ``t_ref_K``.
        omega: VHS viscosity-temperature exponent, in ``[0.5, 1.0]``. 0.5 is a
            hard sphere.
        t_ref_K: the reference temperature ``diameter_m`` belongs to [K].
            Irrelevant at ``omega = 0.5``, where the cross-section does not
            depend on temperature.
        internal_degrees_of_freedom: rotational (and any other fully excited)
            modes. 0 for a monatomic gas, 2 for a linear molecule.

    ``diameter_m`` and ``t_ref_K`` are a **pair** whenever ``omega != 0.5``:
    changing one without the other silently changes the gas.
    """

    name: str
    mass_kg: float
    diameter_m: float
    omega: float = 0.5
    t_ref_K: float = 273.0
    internal_degrees_of_freedom: int = 0

    def __post_init__(self) -> None:
        if not SPECIES_NAME.match(str(self.name)):
            raise ValueError(
                f"species name {self.name!r} is not an identifier. It becomes the "
                f"dsmcFoam typeId and part of the file name "
                f"0/boundaryNumberDensity_<name>, so use letters, digits and "
                f"underscores, starting with a letter.")
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
        """Specific gas constant ``R = k_B / m`` [J/(kg K)]."""
        return BOLTZMANN_J_PER_K / self.mass_kg

    @property
    def gamma(self) -> float:
        """Ratio of specific heats, from the degrees of freedom.

        ``(f + 2) / f`` with ``f = 3 + internal_degrees_of_freedom``. Reported,
        not configured: the DSMC gas has exactly the modes ``dsmcProperties``
        gives it, so this is the ``gamma`` the *simulated* gas has -- which is
        not necessarily the one the nozzle flow was computed with, and the
        difference is worth seeing printed.
        """
        dof = 3 + int(self.internal_degrees_of_freedom)
        return (dof + 2) / dof

    @property
    def molar_mass_g_per_mol(self) -> float:
        """Molar mass [g/mol] implied by ``mass_kg``. Reported, never an input."""
        return self.mass_kg * AVOGADRO_PER_MOL * 1000.0

    @property
    def collision_cross_section_m2(self) -> float:
        """``pi d_ref^2`` [m^2], the total cross-section at ``T_ref``.

        At ``omega = 0.5`` this is the cross-section at every temperature.
        """
        return math.pi * self.diameter_m ** 2

    def most_probable_speed(self, temperature_K: float) -> float:
        """``sqrt(2 R T)`` [m/s] -- the speed a speed ratio is normalised by."""
        return math.sqrt(2.0 * self.gas_constant_j_per_kg_k * _positive(
            temperature_K, "temperature_K"))

    def mean_thermal_speed(self, temperature_K: float) -> float:
        """``sqrt(8 R T / pi)`` [m/s], the Maxwellian mean speed."""
        return math.sqrt(
            8.0 * self.gas_constant_j_per_kg_k
            * _positive(temperature_K, "temperature_K") / math.pi)

    def thermal_sigma(self, temperature_K: float) -> float:
        """``sqrt(R T)`` [m/s], the standard deviation of one velocity component."""
        return math.sqrt(self.gas_constant_j_per_kg_k * _positive(
            temperature_K, "temperature_K"))

    def sound_speed(self, temperature_K: float) -> float:
        """``sqrt(gamma R T)`` [m/s], with :attr:`gamma` from the mode count."""
        return math.sqrt(self.gamma * self.gas_constant_j_per_kg_k * _positive(
            temperature_K, "temperature_K"))


def _positive(value: float, what: str) -> float:
    value = float(value)
    if not (value > 0.0 and math.isfinite(value)):
        raise ValueError(f"{what} must be finite and positive, got {value}")
    return value


def _convention(convention: str) -> str:
    if convention not in MEAN_FREE_PATH_CONVENTIONS:
        raise ValueError(
            f"unknown mean-free-path convention {convention!r}; "
            f"valid: {list(MEAN_FREE_PATH_CONVENTIONS)}")
    return convention


def vhs_temperature_factor(species: VhsSpecies, temperature_K: float) -> float:
    """``(T / T_ref)**(omega - 1/2)`` -- the VHS correction in (G1).

    Exactly 1.0 when ``omega == 0.5``, whatever the temperatures are.
    """
    return (_positive(temperature_K, "temperature_K") / species.t_ref_K) ** (
        float(species.omega) - 0.5)


def mean_free_path(number_density_per_m3: float, species: VhsSpecies,
                   temperature_K: float, *, convention: str = "vhs") -> float:
    """Equilibrium mean free path [m], equation (G1).

    Returns ``inf`` for a zero density, which is the correct limit and keeps a
    vacuum cell from raising in a diagnostic loop.

    Raises:
        ValueError: for a negative density or an unknown convention.
    """
    _convention(convention)
    n = float(number_density_per_m3)
    if n < 0.0 or not math.isfinite(n):
        raise ValueError(f"number density must be finite and >= 0, got {n}")
    if n == 0.0:
        return float("inf")

    hard_sphere = 1.0 / (math.sqrt(2.0) * species.collision_cross_section_m2 * n)
    if convention == "hard_sphere":
        return hard_sphere
    return hard_sphere * vhs_temperature_factor(species, temperature_K)


def number_density_from_mean_free_path(mean_free_path_m: float,
                                       species: VhsSpecies,
                                       temperature_K: float, *,
                                       convention: str = "vhs") -> float:
    """Number density [1/m^3] from a mean free path. The inverse of (G1)."""
    _convention(convention)
    lam = _positive(mean_free_path_m, "mean_free_path_m")
    hard_sphere = 1.0 / (math.sqrt(2.0) * species.collision_cross_section_m2 * lam)
    if convention == "hard_sphere":
        return hard_sphere
    return hard_sphere * vhs_temperature_factor(species, temperature_K)


def collision_time(number_density_per_m3: float, species: VhsSpecies,
                   temperature_K: float, *, convention: str = "vhs") -> float:
    """Mean collision time [s], ``lambda / c_bar``.

    The time a DSMC step has to stay under for moving and colliding to be
    separable. ``inf`` in a vacuum.
    """
    lam = mean_free_path(number_density_per_m3, species, temperature_K,
                         convention=convention)
    return lam / species.mean_thermal_speed(temperature_K)


def _flux_terms(speed_ratio: float) -> tuple[float, float]:
    """``(exp(-s^2), sqrt(pi) (1 + erf s))`` -- shared by the three fluxes."""
    s = float(speed_ratio)
    return math.exp(-s * s), math.sqrt(math.pi) * (1.0 + math.erf(s))


def maxwellian_number_flux(number_density_per_m3: float, speed_ratio: float,
                           species: VhsSpecies, temperature_K: float) -> float:
    """One-way number flux through a plane [1/(m^2 s)]. Bird eq. 4.22.

    ``speed_ratio`` is the component **normal to the plane**. This is what
    ``plumeFieldInflow`` injects, so it -- not ``n U`` -- is the injection rate:
    the two differ by 0.04% at ``s = 2`` and by everything at ``s = 0``.
    """
    n = float(number_density_per_m3)
    if n < 0.0:
        raise ValueError(f"number density must be >= 0, got {n}")
    s = float(speed_ratio)
    gaussian, error = _flux_terms(s)
    c_m = species.most_probable_speed(temperature_K)
    return n * c_m / (2.0 * math.sqrt(math.pi)) * (gaussian + s * error)


def maxwellian_momentum_flux(number_density_per_m3: float, speed_ratio: float,
                             species: VhsSpecies, temperature_K: float) -> float:
    """One-way flux of plane-normal momentum [Pa]. Bird eq. 4.25.

    The pressure the crossing molecules would exert on a surface that absorbed
    them. Over a nozzle exit it is the jet's thrust per unit area,
    ``rho (U^2 + R T)`` once ``s >> 1``.
    """
    n = float(number_density_per_m3)
    if n < 0.0:
        raise ValueError(f"number density must be >= 0, got {n}")
    s = float(speed_ratio)
    gaussian, error = _flux_terms(s)
    c_m = species.most_probable_speed(temperature_K)
    return (n * species.mass_kg * c_m ** 2 / (2.0 * math.sqrt(math.pi))
            * (s * gaussian + (0.5 + s * s) * error))


def maxwellian_energy_flux(number_density_per_m3: float, speed_ratio: float,
                           species: VhsSpecies, temperature_K: float) -> float:
    """One-way energy flux through a plane [W/m^2]. Bird eq. 4.27, plus internal.

    Translational energy from the distribution, and the equipartition internal
    energy ``(zeta/2) k T`` carried by each crossing molecule -- the same
    ``equipartitionInternalEnergy`` the inflow model gives a parcel on average.
    """
    n = float(number_density_per_m3)
    if n < 0.0:
        raise ValueError(f"number density must be >= 0, got {n}")
    s = float(speed_ratio)
    gaussian, error = _flux_terms(s)
    c_m = species.most_probable_speed(temperature_K)
    translational = (n * species.mass_kg * c_m ** 3 / (4.0 * math.sqrt(math.pi))
                     * ((s * s + 2.0) * gaussian + s * (s * s + 2.5) * error))
    internal = (maxwellian_number_flux(n, s, species, temperature_K)
                * 0.5 * int(species.internal_degrees_of_freedom)
                * BOLTZMANN_J_PER_K * float(temperature_K))
    return translational + internal


def describe(species: VhsSpecies) -> list[str]:
    """Report lines for the gas model, printed before a case is meshed or run."""
    hard = float(species.omega) == 0.5
    return [
        "Gas model (one VHS species)",
        f"  species                  {species.name}",
        f"  molecular mass           {species.mass_kg:.5e} kg"
        f"   (M = {species.molar_mass_g_per_mol:.4f} g/mol)",
        f"  molecular diameter       {species.diameter_m:.4e} m",
        f"  collision cross-section  {species.collision_cross_section_m2:.5e} m^2"
        f"   = pi d^2   [DERIVED]",
        f"  omega                    {species.omega:g}"
        + ("   (hard sphere: the cross-section does not vary with T)" if hard
           else f"   (T_ref = {species.t_ref_K:g} K)"),
        f"  internal DoF             {int(species.internal_degrees_of_freedom)}"
        f"   -> gamma = {species.gamma:.4f} for the SIMULATED gas",
        f"  R = k_B/m                {species.gas_constant_j_per_kg_k:.3f} J/(kg K)",
    ]
