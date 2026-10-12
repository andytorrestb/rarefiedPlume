r"""The exit state, and everything it implies.

Three numbers define the inflow -- a number density, a bulk speed and a
temperature, uniform across the exit disk -- and this module is where they
become a single object, however ``case.yaml`` chose to state them:

.. code-block:: text

    density   number_density_per_m3          n
              mass_density_kg_per_m3         n = rho / m
              knudsen                        n from lambda = Kn D      (G1)

    speed     speed_m_per_s                  U
              speed_ratio                    U = S sqrt(2 R T)

A nozzle code or a data sheet states the first form of each; a similarity study
such as Cai & Wang's states the last. Both arrive at the same
:class:`ExitState`, and :meth:`ExitState.describe` prints every other form next
to the one that was given, so a value typed from a data sheet can be checked
against the rest of that sheet's row.

The derived quantities are the ones worth checking a new case against before
any of it is meshed: the **Knudsen number** says how far from resolvable the exit
is, and the **mass flow** and **thrust** are numbers an engine has a name plate
for.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from plumetools.nozzleplume import gas as gas_module
from plumetools.nozzleplume.constants import (
    BOLTZMANN_J_PER_K,
    STANDARD_GRAVITY_M_PER_S2,
)
from plumetools.nozzleplume.gas import VhsSpecies


@dataclass(frozen=True)
class ExitState:
    """The uniform drifting Maxwellian on the nozzle exit disk.

    Attributes:
        species: the VHS species.
        temperature_K: exit static temperature [K].
        speed_m_per_s: bulk speed ``U``, along ``+x`` [m/s].
        number_density_per_m3: ``n`` [1/m^3].
        diameter_m: the exit diameter ``D`` [m] -- the length the Knudsen number
            is built on.
        convention: the mean-free-path convention.
        density_given_as, speed_given_as: which ``case.yaml`` key supplied each,
            so a report can say which numbers were inputs.
    """

    species: VhsSpecies
    temperature_K: float
    speed_m_per_s: float
    number_density_per_m3: float
    diameter_m: float
    convention: str = "vhs"
    density_given_as: str = "number_density_per_m3"
    speed_given_as: str = "speed_m_per_s"

    # --- the state itself ----------------------------------------------------
    @property
    def mass_density_kg_per_m3(self) -> float:
        """``rho = n m`` [kg/m^3]."""
        return self.number_density_per_m3 * self.species.mass_kg

    @property
    def static_pressure_pa(self) -> float:
        """``p = n k T`` [Pa]."""
        return self.number_density_per_m3 * BOLTZMANN_J_PER_K * self.temperature_K

    @property
    def dynamic_pressure_pa(self) -> float:
        """``rho U^2 / 2`` [Pa]."""
        return 0.5 * self.mass_density_kg_per_m3 * self.speed_m_per_s ** 2

    @property
    def area_m2(self) -> float:
        """``pi D^2 / 4`` [m^2], the full exit area."""
        return 0.25 * math.pi * self.diameter_m ** 2

    # --- speeds --------------------------------------------------------------
    @property
    def most_probable_speed_m_per_s(self) -> float:
        """``sqrt(2 R T)`` [m/s]."""
        return self.species.most_probable_speed(self.temperature_K)

    @property
    def mean_thermal_speed_m_per_s(self) -> float:
        """``sqrt(8 R T / pi)`` [m/s]."""
        return self.species.mean_thermal_speed(self.temperature_K)

    @property
    def thermal_sigma_m_per_s(self) -> float:
        """Per-component thermal standard deviation ``sqrt(R T)`` [m/s]."""
        return self.species.thermal_sigma(self.temperature_K)

    @property
    def characteristic_speed_m_per_s(self) -> float:
        """``U + 3 sigma`` [m/s] -- what the time step must not outrun.

        Three standard deviations rather than the mean: the mean would let the
        fast tail cross several cells per step unnoticed, and it is the tail that
        breaks the collision sampling.
        """
        return self.speed_m_per_s + 3.0 * self.thermal_sigma_m_per_s

    @property
    def speed_ratio(self) -> float:
        """``S = U / sqrt(2 R T)``."""
        return self.speed_m_per_s / self.most_probable_speed_m_per_s

    @property
    def mach_number(self) -> float:
        """``U / sqrt(gamma R T)``, with ``gamma`` of the **simulated** gas."""
        return self.speed_m_per_s / self.species.sound_speed(self.temperature_K)

    # --- rarefaction ---------------------------------------------------------
    @property
    def mean_free_path_m(self) -> float:
        """Exit mean free path ``lambda`` [m]."""
        return gas_module.mean_free_path(
            self.number_density_per_m3, self.species, self.temperature_K,
            convention=self.convention)

    @property
    def knudsen(self) -> float:
        """``Kn = lambda / D``, on the exit diameter."""
        return self.mean_free_path_m / self.diameter_m

    @property
    def collision_time_s(self) -> float:
        """Exit mean collision time ``lambda / c_bar`` [s]."""
        return self.mean_free_path_m / self.mean_thermal_speed_m_per_s

    # --- fluxes through the exit plane --------------------------------------
    @property
    def number_flux_per_m2_s(self) -> float:
        """One-way number flux through the exit plane [1/(m^2 s)]."""
        return gas_module.maxwellian_number_flux(
            self.number_density_per_m3, self.speed_ratio, self.species,
            self.temperature_K)

    @property
    def momentum_flux_pa(self) -> float:
        """One-way axial momentum flux through the exit plane [Pa]."""
        return gas_module.maxwellian_momentum_flux(
            self.number_density_per_m3, self.speed_ratio, self.species,
            self.temperature_K)

    @property
    def energy_flux_w_per_m2(self) -> float:
        """One-way energy flux through the exit plane [W/m^2]."""
        return gas_module.maxwellian_energy_flux(
            self.number_density_per_m3, self.speed_ratio, self.species,
            self.temperature_K)

    def injection_rate_per_s(self, area_m2: float | None = None) -> float:
        """Molecules per second entering through an exit of ``area_m2``.

        Pass the **meshed** patch area to get what the solver actually injects:
        the staircase patch is not the disk, and the difference is the mass-flow
        error of the mesh. With no argument, the exact exit area.
        """
        return self.number_flux_per_m2_s * float(
            self.area_m2 if area_m2 is None else area_m2)

    def mass_flow_kg_per_s(self, area_m2: float | None = None) -> float:
        """Mass flow through the exit [kg/s]."""
        return self.injection_rate_per_s(area_m2) * self.species.mass_kg

    def thrust_n(self, area_m2: float | None = None) -> float:
        """Axial momentum leaving the exit per second [N].

        ``rho (U^2 + R T) A`` at a supersonic exit: momentum thrust plus
        pressure thrust into vacuum.
        """
        return self.momentum_flux_pa * float(
            self.area_m2 if area_m2 is None else area_m2)

    def power_w(self, area_m2: float | None = None) -> float:
        """Energy leaving the exit per second [W]."""
        return self.energy_flux_w_per_m2 * float(
            self.area_m2 if area_m2 is None else area_m2)

    @property
    def specific_impulse_s(self) -> float:
        """``thrust / (mass flow * g0)`` [s]. Independent of the area."""
        return self.momentum_flux_pa / (
            self.number_flux_per_m2_s * self.species.mass_kg
            * STANDARD_GRAVITY_M_PER_S2)

    # --- reporting -----------------------------------------------------------
    def as_dict(self) -> dict:
        """Machine-readable, for the case summary."""
        return {
            "species": self.species.name,
            "diameter_m": float(self.diameter_m),
            "area_m2": float(self.area_m2),
            "temperature_K": float(self.temperature_K),
            "speed_m_per_s": float(self.speed_m_per_s),
            "number_density_per_m3": float(self.number_density_per_m3),
            "mass_density_kg_per_m3": float(self.mass_density_kg_per_m3),
            "static_pressure_pa": float(self.static_pressure_pa),
            "dynamic_pressure_pa": float(self.dynamic_pressure_pa),
            "speed_ratio": float(self.speed_ratio),
            "mach_number": float(self.mach_number),
            "gamma_simulated": float(self.species.gamma),
            "mean_free_path_m": float(self.mean_free_path_m),
            "knudsen": float(self.knudsen),
            "collision_time_s": float(self.collision_time_s),
            "number_flux_per_m2_s": float(self.number_flux_per_m2_s),
            "mass_flow_kg_per_s": float(self.mass_flow_kg_per_s()),
            "thrust_n": float(self.thrust_n()),
            "specific_impulse_s": float(self.specific_impulse_s),
            "power_w": float(self.power_w()),
            "mean_free_path_convention": self.convention,
            "density_given_as": self.density_given_as,
            "speed_given_as": self.speed_given_as,
        }

    def describe(self) -> list[str]:
        """Report lines, printed before meshing and before running.

        Each value is tagged ``[GIVEN]`` if ``case.yaml`` supplied it and
        ``[DERIVED]`` otherwise, so the report reads as one row of a data sheet
        with the inputs marked.
        """
        def tag(*keys: str) -> str:
            given = {self.density_given_as, self.speed_given_as}
            return "[GIVEN]" if given & set(keys) else "[DERIVED]"

        return [
            "Exit state (uniform drifting Maxwellian on the exit disk)",
            f"  exit diameter D          {self.diameter_m:.6g} m   [GIVEN]",
            f"  exit area                {self.area_m2:.6f} m^2"
            f"   = pi D^2 / 4   [DERIVED]",
            f"  temperature T            {self.temperature_K:g} K   [GIVEN]",
            f"  bulk speed U             {self.speed_m_per_s:.6g} m/s"
            f"   {tag('speed_m_per_s')}",
            f"  speed ratio S            {self.speed_ratio:.4f}"
            f"   = U / sqrt(2 R T)   {tag('speed_ratio')}",
            f"  Mach number              {self.mach_number:.3f}"
            f"   (gamma = {self.species.gamma:.4f}, the simulated gas)",
            f"  number density n         {self.number_density_per_m3:.4e} 1/m^3"
            f"   {tag('number_density_per_m3')}",
            f"  mass density rho         {self.mass_density_kg_per_m3:.4e} kg/m^3"
            f"   = n m   {tag('mass_density_kg_per_m3')}",
            f"  static pressure          {self.static_pressure_pa:.4g} Pa"
            f"   = n k T",
            f"  mean free path           {self.mean_free_path_m:.4e} m"
            f"   ({self.convention} convention)",
            f"  Knudsen number           {self.knudsen:.4e}"
            f"   = lambda / D   {tag('knudsen')}",
            f"  mean collision time      {self.collision_time_s:.4e} s",
            f"  number flux              {self.number_flux_per_m2_s:.4e} 1/(m^2 s)",
            f"  mass flow                {self.mass_flow_kg_per_s():.5g} kg/s"
            f"   through the exact exit area",
            f"  thrust                   {self.thrust_n():.5g} N"
            f"   = rho (U^2 + R T) A;  Isp {self.specific_impulse_s:.1f} s",
            f"  jet power                {self.power_w():.4e} W",
        ]


def from_config(cfg) -> ExitState:
    """Build the exit state from a validated ``case.yaml``.

    :func:`plumetools.nozzleplume.config.validate` has already insisted on
    exactly one density key and exactly one speed key, so each branch below is
    the only one that can apply.
    """
    species = cfg.species()
    temperature = float(cfg.exit.temperature_K)
    diameter = float(cfg.nozzle.exit_diameter_m)
    convention = cfg.gas.mean_free_path_convention

    if cfg.exit.speed_m_per_s is not None:
        speed, speed_given_as = float(cfg.exit.speed_m_per_s), "speed_m_per_s"
    else:
        speed = (float(cfg.exit.speed_ratio)
                 * species.most_probable_speed(temperature))
        speed_given_as = "speed_ratio"

    if cfg.exit.number_density_per_m3 is not None:
        density = float(cfg.exit.number_density_per_m3)
        density_given_as = "number_density_per_m3"
    elif cfg.exit.mass_density_kg_per_m3 is not None:
        density = float(cfg.exit.mass_density_kg_per_m3) / species.mass_kg
        density_given_as = "mass_density_kg_per_m3"
    else:
        density = gas_module.number_density_from_mean_free_path(
            float(cfg.exit.knudsen) * diameter, species, temperature,
            convention=convention)
        density_given_as = "knudsen"

    return ExitState(
        species=species,
        temperature_K=temperature,
        speed_m_per_s=speed,
        number_density_per_m3=density,
        diameter_m=diameter,
        convention=convention,
        density_given_as=density_given_as,
        speed_given_as=speed_given_as,
    )
