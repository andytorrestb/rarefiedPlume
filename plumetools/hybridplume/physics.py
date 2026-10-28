r"""The single source of every derived physical quantity in the hybrid plume.

Nothing else in :mod:`plumetools.hybridplume` derives a density, a speed or a
transport coefficient: it asks a :class:`SourceState` or a :class:`Gas`. That is
what keeps the CFD precursor, the DSMC inflow, the collisionless reference and
the diagnostics describing the *same* gas.

Coordinates
-----------
The source disk occupies ``x = 0, y^2 + z^2 <= Re^2`` and the plume expands
into ``+x``. All geometry is full three-dimensional Cartesian; there is no
axisymmetric weighting anywhere in this package.

The gas: what is inherited, what is changed, and why
----------------------------------------------------
=====================  ==============  ======================================
quantity               value           provenance
=====================  ==============  ======================================
molecular mass ``m``   6.63e-26 kg     repository (``cai2012.gas.ARGON``);
                                       39.948 u is 6.6335e-26 kg
VHS ``d_ref``          4.17e-10 m      repository = Bird (1994) Table A1
``T_ref``              273 K           repository = Bird (1994) Table A1
VHS ``omega``          **0.81**        Bird (1994) Table A1 argon -- **not**
                                       the repository's 0.74
``gamma``              5/3             monatomic, no internal modes
``Pr``                 2/3             monatomic Eucken value
=====================  ==============  ======================================

The Cai-family packages use ``omega = 0.74``, which ``cai2012/gas.py`` itself
records as the OpenFOAM N2 tutorial exponent rather than an argon value. There
it enters only through a 0.4 % change in ``n0``. Here it sets the CFD
viscosity over a temperature range of more than a decade, so the hybrid family
uses Bird's argon pair ``(d_ref, omega) = (4.17e-10 m, 0.81)`` and records the
departure. ``gas.omega: 0.74`` in the configuration restores the Cai-family
value; nothing else changes.

Transport consistent with the DSMC collision model
--------------------------------------------------
OpenFOAM's ``VariableHardSphere::sigmaTcR`` (v2512) uses Bird's VHS
cross-section

.. math::

    \sigma_T = \pi d_{ref}^2
        \left(\frac{2kT_{ref}}{m_r c_r^2}\right)^{\omega-1/2}
        \Big/ \Gamma(5/2-\omega).

The first Chapman--Enskog approximation for that cross-section is Bird (1994)
eq. 4.62, which relates ``d_ref`` and the viscosity at ``T_ref``:

.. math::

    \mu_{ref} = \frac{15\sqrt{m k T_{ref}/\pi}}
                     {2(5-2\omega)(7-2\omega)\,d_{ref}^2},
    \qquad \mu(T) = \mu_{ref}(T/T_{ref})^\omega .

This is checked against Bird's own table: argon at (4.17e-10 m, 0.81, 273 K)
gives 2.117e-5 Pa s, the tabulated value. The conductivity is
``kappa = mu cp / Pr`` with ``Pr = 2/3`` (Eucken, exact for the first
Chapman--Enskog approximation of a monatomic gas).

Mean-free-path convention
-------------------------
The diagnostic convention, named ``vhs``, is Bird eq. 4.76 written with the
viscosity of eq. 4.62 substituted, which collapses to

.. math::

    \lambda(n,T) = \frac{(T/T_{ref})^{\omega-1/2}}{\sqrt2\,\pi d_{ref}^2\,n}.

It is a *diagnostic* length, not the free flight of every molecule: in an
expanding plume the distribution is not a single Maxwellian. Two other
conventions appear in the literature and differ by a constant factor; both
are reported by :meth:`Gas.mean_free_path_conventions` so a threshold quoted
under one can be converted:

* ``hard_sphere``: ``1/(sqrt(2) pi d^2 n)``, the ``vhs`` value at ``T_ref``;
* ``viscosity_simple``: ``(mu/p) sqrt(pi R T / 2)``, larger than ``vhs`` by
  ``15 pi / (2 (5-2w)(7-2w))`` = 1.296 for ``omega = 0.81``.

The source state
----------------
From ``D``, ``S0``, ``Te``, ``Kne`` and the gas::

    ae   = sqrt(2 R Te)              thermal speed scale (NOT the sound speed)
    Ue   = S0 ae                     bulk velocity, along +x
    ce   = sqrt(gamma R Te)          sound speed
    Me   = Ue/ce = S0 sqrt(2/gamma)
    lambda_e = Kne D
    ne   = (Te/Tref)^(omega-1/2) / (sqrt(2) pi d^2 lambda_e)
    rho_e = m ne,   pe = ne k Te
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np

#: Boltzmann constant [J/K] (SI 2019, exact).
K_BOLTZMANN = 1.380649e-23
#: Avogadro constant [1/mol] (SI 2019, exact).
AVOGADRO = 6.02214076e23

MEAN_FREE_PATH_CONVENTIONS = ("vhs", "hard_sphere", "viscosity_simple")


def _positive(value, what: str) -> float:
    value = float(value)
    if not (math.isfinite(value) and value > 0.0):
        raise ValueError(f"{what} must be finite and positive, got {value!r}")
    return value


@dataclass(frozen=True)
class Gas:
    """A single monatomic VHS species, and the continuum model consistent with it.

    Attribute names match what :func:`plumetools.rarefaction.core.mean_free_path`
    reads (``diameter_m``, ``omega``, ``t_ref_K``), so the generic diagnostics
    can be pointed at it unchanged.
    """

    name: str = "Ar"
    mass_kg: float = 6.63e-26
    diameter_m: float = 4.17e-10
    t_ref_K: float = 273.0
    omega: float = 0.81
    prandtl: float = 2.0 / 3.0
    provenance: dict = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        for name in ("mass_kg", "diameter_m", "t_ref_K", "prandtl"):
            _positive(getattr(self, name), f"gas.{name}")
        if not 0.5 <= float(self.omega) <= 1.0:
            raise ValueError(
                f"gas.omega = {self.omega}: the VHS exponent lies in [0.5, 1]")

    # -- thermodynamics: calorically perfect monatomic gas -------------------
    @property
    def gas_constant(self) -> float:
        """Specific gas constant ``R = k_B/m`` [J/(kg K)]."""
        return K_BOLTZMANN / self.mass_kg

    @property
    def gamma(self) -> float:
        return 5.0 / 3.0

    @property
    def cv(self) -> float:
        return 1.5 * self.gas_constant

    @property
    def cp(self) -> float:
        return 2.5 * self.gas_constant

    @property
    def molar_mass_g_per_mol(self) -> float:
        return self.mass_kg * AVOGADRO * 1.0e3

    # -- transport consistent with VHS ---------------------------------------
    @property
    def mu_ref(self) -> float:
        """Viscosity at ``T_ref`` [Pa s], Bird (1994) eq. 4.62 for VHS."""
        w = float(self.omega)
        return (15.0 * math.sqrt(self.mass_kg * K_BOLTZMANN * self.t_ref_K / math.pi)
                / (2.0 * (5.0 - 2.0 * w) * (7.0 - 2.0 * w) * self.diameter_m ** 2))

    def viscosity(self, temperature):
        """``mu(T) = mu_ref (T/T_ref)^omega`` [Pa s]."""
        t = np.asarray(temperature, dtype=np.float64)
        return self.mu_ref * (t / self.t_ref_K) ** self.omega

    def conductivity(self, temperature):
        """``kappa = mu cp / Pr`` [W/(m K)]."""
        return self.viscosity(temperature) * self.cp / self.prandtl

    # -- kinetic scales ------------------------------------------------------
    def thermal_speed(self, temperature):
        """``a = sqrt(2 R T)`` [m/s] -- the most probable speed."""
        return np.sqrt(2.0 * self.gas_constant * np.asarray(temperature, dtype=np.float64))

    def sound_speed(self, temperature):
        """``c = sqrt(gamma R T)`` [m/s]."""
        return np.sqrt(self.gamma * self.gas_constant * np.asarray(temperature, dtype=np.float64))

    def mean_speed(self, temperature):
        """``c_bar = sqrt(8 R T / pi)`` [m/s]."""
        return np.sqrt(8.0 * self.gas_constant * np.asarray(temperature, dtype=np.float64) / math.pi)

    def mean_free_path(self, number_density, temperature, convention: str = "vhs"):
        """Mean free path [m] under a named convention; NaN where undefined."""
        n = np.asarray(number_density, dtype=np.float64)
        t = np.asarray(temperature, dtype=np.float64)
        ok = (n > 0) & np.isfinite(n) & (t > 0) & np.isfinite(t)
        with np.errstate(divide="ignore", invalid="ignore"):
            hard = 1.0 / (math.sqrt(2.0) * math.pi * self.diameter_m ** 2 * n)
            if convention == "hard_sphere":
                lam = hard
            elif convention == "vhs":
                lam = hard * (t / self.t_ref_K) ** (self.omega - 0.5)
            elif convention == "viscosity_simple":
                p = n * K_BOLTZMANN * t
                lam = self.viscosity(t) / p * np.sqrt(math.pi * self.gas_constant * t / 2.0)
            else:
                raise ValueError(f"unknown mean-free-path convention {convention!r}; "
                                 f"known: {MEAN_FREE_PATH_CONVENTIONS}")
        return np.where(ok, lam, np.nan)

    def number_density_for_mean_free_path(self, mean_free_path, temperature,
                                          convention: str = "vhs") -> float:
        """Invert :meth:`mean_free_path` for ``n`` [1/m^3]."""
        lam = _positive(mean_free_path, "mean free path")
        t = _positive(temperature, "temperature")
        unit = float(self.mean_free_path(1.0, t, convention))
        return unit / lam

    def collision_time(self, number_density, temperature, convention: str = "vhs"):
        """Nominal ``tau_c = lambda / c_bar`` [s]."""
        return (self.mean_free_path(number_density, temperature, convention)
                / self.mean_speed(temperature))

    def mean_free_path_conventions(self, temperature) -> dict:
        """Ratios of the other conventions to ``vhs`` at one temperature."""
        t = _positive(temperature, "temperature")
        base = float(self.mean_free_path(1.0, t, "vhs"))
        return {c: float(self.mean_free_path(1.0, t, c)) / base
                for c in MEAN_FREE_PATH_CONVENTIONS}

    def describe(self) -> dict:
        out = {k: v for k, v in asdict(self).items() if k != "provenance"}
        out.update(gas_constant_J_per_kgK=self.gas_constant, gamma=self.gamma,
                   cv=self.cv, cp=self.cp, mu_ref_Pa_s=self.mu_ref,
                   molar_mass_g_per_mol=self.molar_mass_g_per_mol,
                   provenance=dict(self.provenance))
        return out


#: Bird (1994), Table A1, argon. ``mass_kg`` is the repository's value.
BIRD_ARGON = Gas(
    name="Ar", mass_kg=6.63e-26, diameter_m=4.17e-10, t_ref_K=273.0, omega=0.81,
    prandtl=2.0 / 3.0,
    provenance={
        "mass_kg": "repository (plumetools.cai2012.gas.ARGON); 39.948 u = 6.6335e-26 kg",
        "diameter_m": "Bird (1994) Table A1 argon VHS; same as the repository",
        "t_ref_K": "Bird (1994) Table A1; same as the repository",
        "omega": "Bird (1994) Table A1 argon; the Cai family uses 0.74 (N2 tutorial value)",
        "prandtl": "monatomic Eucken value 2/3",
    })


@dataclass(frozen=True)
class SourceState:
    """The uniform drifting-Maxwellian circular source, and everything derived from it.

    Inputs are the four Cai quantities plus the gas. Every other number is a
    property, so no derived value can be stored inconsistently.
    """

    gas: Gas
    diameter_m: float = 0.2
    speed_ratio: float = 2.0
    temperature_K: float = 300.0
    knudsen: float = 0.01
    mean_free_path_convention: str = "vhs"

    def __post_init__(self) -> None:
        _positive(self.diameter_m, "source.diameter_m")
        _positive(self.temperature_K, "source.temperature_K")
        _positive(self.knudsen, "source.knudsen")
        if not (math.isfinite(self.speed_ratio) and self.speed_ratio >= 0.0):
            raise ValueError(f"source.speed_ratio must be >= 0, got {self.speed_ratio}")
        if self.mean_free_path_convention not in MEAN_FREE_PATH_CONVENTIONS:
            raise ValueError(f"unknown mean_free_path_convention "
                             f"{self.mean_free_path_convention!r}")

    @property
    def radius_m(self) -> float:
        return 0.5 * self.diameter_m

    @property
    def area_m2(self) -> float:
        """Exact disk area ``pi D^2 / 4`` [m^2]."""
        return math.pi * self.radius_m ** 2

    @property
    def thermal_speed(self) -> float:
        return float(self.gas.thermal_speed(self.temperature_K))

    @property
    def velocity(self) -> float:
        return self.speed_ratio * self.thermal_speed

    @property
    def velocity_vector(self) -> np.ndarray:
        return np.array([self.velocity, 0.0, 0.0])

    @property
    def sound_speed(self) -> float:
        return float(self.gas.sound_speed(self.temperature_K))

    @property
    def mach(self) -> float:
        return self.speed_ratio * math.sqrt(2.0 / self.gas.gamma)

    @property
    def mean_free_path(self) -> float:
        return self.knudsen * self.diameter_m

    @property
    def number_density(self) -> float:
        return self.gas.number_density_for_mean_free_path(
            self.mean_free_path, self.temperature_K, self.mean_free_path_convention)

    @property
    def mass_density(self) -> float:
        return self.gas.mass_kg * self.number_density

    @property
    def pressure(self) -> float:
        return self.number_density * K_BOLTZMANN * self.temperature_K

    @property
    def viscosity(self) -> float:
        return float(self.gas.viscosity(self.temperature_K))

    @property
    def reynolds(self) -> float:
        return self.mass_density * self.velocity * self.diameter_m / self.viscosity

    @property
    def collision_time(self) -> float:
        return float(self.gas.collision_time(self.number_density, self.temperature_K,
                                             self.mean_free_path_convention))

    @property
    def stagnation_temperature(self) -> float:
        return self.temperature_K + self.velocity ** 2 / (2.0 * self.gas.cp)

    @property
    def continuum_mass_flow(self) -> float:
        """``rho_e Ue A`` [kg/s] -- NOT the kinetic source's gross injection."""
        return self.mass_density * self.velocity * self.area_m2

    def with_knudsen(self, knudsen: float) -> "SourceState":
        return SourceState(self.gas, self.diameter_m, self.speed_ratio,
                           self.temperature_K, knudsen, self.mean_free_path_convention)

    def describe(self) -> dict:
        return {
            "inputs": {
                "diameter_m": self.diameter_m, "speed_ratio": self.speed_ratio,
                "temperature_K": self.temperature_K, "knudsen": self.knudsen,
                "mean_free_path_convention": self.mean_free_path_convention,
            },
            "derived": {
                "radius_m": self.radius_m, "area_m2": self.area_m2,
                "thermal_speed_m_s": self.thermal_speed,
                "velocity_m_s": self.velocity, "sound_speed_m_s": self.sound_speed,
                "mach": self.mach, "mean_free_path_m": self.mean_free_path,
                "number_density_per_m3": self.number_density,
                "mass_density_kg_m3": self.mass_density,
                "pressure_Pa": self.pressure, "viscosity_Pa_s": self.viscosity,
                "reynolds_D": self.reynolds,
                "collision_time_s": self.collision_time,
                "stagnation_temperature_K": self.stagnation_temperature,
                "continuum_mass_flow_kg_s": self.continuum_mass_flow,
            },
            "gas": self.gas.describe(),
            "mean_free_path_convention_ratios": self.gas.mean_free_path_conventions(
                self.temperature_K),
        }
