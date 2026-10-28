r"""The Cai exit continued as a Markelov angular source, to a downstream sphere.

At ``Kn_e = 1e-3`` and ``1e-4`` the dense core next to Cai's exit cannot be
resolved by a 3-D DSMC mesh (``docs/cai2012-kn-extension.md`` §5). The
remedy this module supports is to remove that core from the DSMC domain and
inject, on a sphere downstream, the state a continuum expansion of the same
exit would have reached there. Everything here is pure: no I/O, no OpenFOAM.

Model identifier: ``cai2012ext_markelov_streamtube``.

What is the paper's, and what is this module's
----------------------------------------------
``[PAPER]`` -- Lumpkin, Stewart & Markelov, AIAA 99-3455, printed p. 3: the
angular law, its normalisation and the limiting speed. Implemented once, in
:mod:`plumetools.markelov1999.sourceflow`, and **called** from here, never
re-typed::

    theta_L = (pi/2) [sqrt((g+1)/(g-1)) - 1]                 sf.limiting_angle
    f(theta) = cos(pi theta / (2 theta_L))**((g+0.41)/(g-1))  sf.angular
    I = int_0^theta_L f sin(theta) dtheta                     sf.normalization_integral
    A_M = 0.5 sqrt((g-1)/(g+1)) / I                           sf.normalization_coefficient
    V_inf = sqrt(2 cp T_t)                                    sf.limiting_velocity

``[PAPER]`` -- Stewart & Lumpkin (2012), Eq. (2): the spherical area-Mach,
temperature and density relations; Eq. (1): the Bird parameter; Eqs. (3)-(4):
its simplified radial form and the hard-sphere inversion. :func:`area_mach`,
:func:`bird_parameter`, :func:`stewart_hard_sphere_density`.

``[ASSUMPTION]`` -- this module's, and nobody's published result. They join a
point-source model to a finite circular exit, which neither paper does:

=====  =======================================================================
CP-1   The source sits at the exit centre, ``x_s = x_e``, on the exit axis.
CP-2   The conserved targets are the exit's GROSS outward molecular fluxes
       (half-range drifting Maxwellian, :func:`exit_targets`), not its
       continuum net flux. The two differ by 0.02% at S0 = 2; both reported.
CP-3   Mass flux per steradian has Markelov's shape at every radius:
       ``j(theta) = mdot f(theta) / (2 pi I)``. The tail beyond the core is not
       redistributed into it.
CP-4   Independent radial streamtubes, each an isentropic quasi-1-D expansion
       through area ``r^2``: ``rho U r^2 = j``, ``cp T + U^2/2 = h_target``.
       No transverse pressure gradient, no finite-exit structure.
CP-5   The closure anchor: ``U(r_a) = U_e`` and ``rho(r_a, 0) = rho_e``. It
       defines an equivalent sphere ``r_a``, not a geometric nozzle exit.
=====  =======================================================================

The resulting fields, with ``T_t = h_target/cp``, ``T_a = T_t - U_e^2/(2 cp)``,
``M_a = U_e/sqrt(g R T_a)`` and ``r_a = sqrt(mdot / (2 pi I rho_e U_e))``::

    F(M) = (1/M) [(2/(g+1)) (1 + (g-1) M^2/2)]**((g+1)/(2(g-1)))
    F(M(r)) = F(M_a) (r/r_a)^2                    supersonic branch only
    T(r) = T_t / (1 + (g-1) M^2/2),  U(r) = M sqrt(g R T),  U_vec = U r_hat
    rho(r, theta) = j(theta) / (U r^2),  n = rho/m

Along a ray ``n/n(r_a) = (T/T_a)**(1/(g-1))``: it IS the isentropic expansion.

The 1999 terminal model, compared
---------------------------------
Mapped through the same mass flow, Markelov's printed density is the
``U -> V_inf`` limit of this closure: ``rho_M = j/(V_inf r^2)``. Only the
product ``p_t r_o^2`` is fixed by the mass flow, and :func:`markelov_terminal`
evaluates the frozen :func:`~plumetools.markelov1999.sourceflow.mass_density`
with it, so the comparison runs through the published implementation. The
terminal model's own velocity ``V_inf`` and temperature ``T_t`` together carry
twice the available enthalpy; that is why it is compared and never injected.

Bird's parameter
----------------
``B = -(lambda/(rho c_bar)) U . grad(rho)``, signed, positive in an expansion.
Along a ray of this closure ``d ln U/dr = 2/(r (M^2 - 1))`` exactly, so::

    B_general       = (2 lambda U / (r c_bar)) M^2 / (M^2 - 1)     default
    B_radial_approx = 2 lambda U / (r c_bar)                       Stewart Eq. (3)

``lambda`` is the case's own convention (:mod:`plumetools.cai2012.gas`, VHS
with the local temperature by default). The two definitions differ by
``M^2/(M^2-1)`` -- 26% at the anchor -- and every surface records which one
chose it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import integrate, optimize, special

from plumetools.cai2012 import inflow as cai_inflow
from plumetools.cai2012.geometry import from_config as cai_geometry
# The configuration -- interface.yaml's schema, validation and reading -- is in
# interface_config; its names stay importable from here.
from plumetools.cai2012ext.interface_config import (  # noqa: F401
    BIRD_DEFINITIONS,
    PURPOSES,
    RUN_MODES,
    SECTIONS,
    SETTINGS_FILE,
    SURFACES,
    TAIL_TREATMENTS,
    InterfaceConfig,
    InterfaceConfigError,
    load_interface,
    parse_interface,
    read_settings,
)
from plumetools.markelov1999 import sourceflow as sf

#: The model identifier recorded in every manifest row and report.
MODEL_NAME = "cai2012ext_markelov_streamtube"

_SQRT_PI = math.sqrt(math.pi)

#: Selection outcomes.
STATUS_ADMISSIBLE = "ADMISSIBLE"
STATUS_NO_INTERFACE = "NO_ADMISSIBLE_INTERFACE"
STATUS_INVALID_STATE = "INVALID_STATE"
STATUS_INSUFFICIENT_DOMAIN = "INSUFFICIENT_DOMAIN"
STATUS_INCOMPLETE = "CONFIGURATION_INCOMPLETE"


# --------------------------------------------------------------------------- #
# B -- one-way molecular fluxes of a drifting Maxwellian
# --------------------------------------------------------------------------- #

def one_way_moments(u_n, temperature_K, gas_constant: float):
    r"""``(J1, J3)``: the positive normal-velocity moments per unit density.

    With ``a = sqrt(2 R T)``, ``s = u_n/a``, ``H = erfc(-s)/2`` and
    ``E = exp(-s^2)/(2 sqrt(pi))``::

        J1 = u_n H + a E                                 [m/s]
        J3 = (u_n^3 + 1.5 u_n a^2) H + a (u_n^2 + a^2) E   [m^3/s^3]

    so the gross number flux is ``n J1`` and the gross translational energy
    flux ``(m n/2)[J3 + (|U_t|^2 + 2 R T) J1]``.

    For ``s < 0`` both are differences of nearly equal terms. They are
    evaluated with the scaled ``erfcx(-s) = exp(s^2) erfc(-s)`` and the common
    factor ``exp(-s^2)`` taken out, so the relative error stays near machine
    precision for every ``s`` at which the result is representable (below
    ``s ~ -27`` it underflows to an exact zero).
    """
    u = np.asarray(u_n, dtype=np.float64)
    t = np.asarray(temperature_K, dtype=np.float64)
    a = np.sqrt(2.0 * float(gas_constant) * t)
    s = u / a
    j1 = np.empty(np.broadcast(u, t).shape)
    j3 = np.empty_like(j1)
    s, a, u = np.broadcast_arrays(s, a, u)

    positive = s >= 0.0
    if np.any(positive):
        sp, ap, up = s[positive], a[positive], u[positive]
        h = 0.5 * special.erfc(-sp)
        e = np.exp(-sp * sp) / (2.0 * _SQRT_PI)
        j1[positive] = up * h + ap * e
        j3[positive] = (up ** 3 + 1.5 * up * ap * ap) * h + ap * (up * up + ap * ap) * e
    negative = ~positive
    if np.any(negative):
        sn, an = s[negative], a[negative]
        scale = np.exp(-sn * sn)
        x = special.erfcx(-sn)
        g1 = sn * x + 1.0 / _SQRT_PI
        g3 = (sn ** 3 + 1.5 * sn) * x + (sn * sn + 1.0) / _SQRT_PI
        j1[negative] = 0.5 * an * scale * g1
        j3[negative] = 0.5 * an ** 3 * scale * g3
    return j1, j3


@dataclass(frozen=True)
class ExitTargets:
    """The physical exit's fluxes: what the equivalent source must carry.

    Attributes:
        area_m2: ``A_e = pi D^2/4``.
        number_flux_per_m2_s: ``Gamma_plus = n_e J1(U_e, T_e)``.
        energy_flux_W_per_m2: ``Q_plus``.
        mdot_kg_per_s: ``A_e m Gamma_plus`` (the GROSS target, CP-2).
        edot_W: ``A_e Q_plus``.
        h_J_per_kg: ``edot/mdot``.
        mdot_net_kg_per_s: ``rho_e U_e A_e``, the continuum net flow.
        h_net_J_per_kg: ``cp T_e + U_e^2/2``.
    """

    area_m2: float
    number_flux_per_m2_s: float
    energy_flux_W_per_m2: float
    mdot_kg_per_s: float
    edot_W: float
    h_J_per_kg: float
    mdot_net_kg_per_s: float
    h_net_J_per_kg: float

    def as_dict(self) -> dict:
        return {
            "area_m2": self.area_m2,
            "gross_number_flux_per_m2_s": self.number_flux_per_m2_s,
            "gross_energy_flux_W_per_m2": self.energy_flux_W_per_m2,
            "mdot_target_kg_per_s": self.mdot_kg_per_s,
            "edot_target_W": self.edot_W,
            "h_target_J_per_kg": self.h_J_per_kg,
            "mdot_continuum_net_kg_per_s": self.mdot_net_kg_per_s,
            "h_continuum_net_J_per_kg": self.h_net_J_per_kg,
            "gross_over_net_mass": self.mdot_kg_per_s / self.mdot_net_kg_per_s,
            "gross_over_net_enthalpy": self.h_J_per_kg / self.h_net_J_per_kg,
            "basis": "gross outward half-range drifting-Maxwellian fluxes "
                     "through the exit disk (assumption CP-2); the continuum "
                     "net values are reported beside them",
        }


def exit_targets(*, number_density, velocity, temperature, mass_kg,
                 gas_constant, cp, area_m2) -> ExitTargets:
    """Gross and net fluxes of the uniform axial exit (``U_t = 0``)."""
    j1, j3 = one_way_moments(velocity, temperature, gas_constant)
    j1, j3 = float(j1), float(j3)
    gamma_plus = float(number_density) * j1
    q_plus = 0.5 * mass_kg * float(number_density) * (
        j3 + 2.0 * gas_constant * float(temperature) * j1)
    mdot = area_m2 * mass_kg * gamma_plus
    edot = area_m2 * q_plus
    return ExitTargets(
        area_m2=float(area_m2),
        number_flux_per_m2_s=gamma_plus,
        energy_flux_W_per_m2=q_plus,
        mdot_kg_per_s=mdot,
        edot_W=edot,
        h_J_per_kg=edot / mdot,
        mdot_net_kg_per_s=mass_kg * float(number_density) * float(velocity) * area_m2,
        h_net_J_per_kg=cp * float(temperature) + 0.5 * float(velocity) ** 2,
    )


# --------------------------------------------------------------------------- #
# D -- area-Mach
# --------------------------------------------------------------------------- #

def area_mach(mach, gamma: float):
    """``F(M)``: the spherical (or quasi-1-D) area ratio ``A/A*``."""
    m = np.asarray(mach, dtype=np.float64)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (1.0 / m) * ((2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * m * m)
                        ) ** exponent


def supersonic_mach(area_ratio, gamma: float, *, lower: float = 1.0):
    """Solve ``F(M) = area_ratio`` on the supersonic branch ``M > lower``.

    Vectorised bisection in ``ln M``: ``F`` is monotone for ``M > 1``, so a
    bracket that holds at both ends holds at every step. The upper end starts
    at ``2 lower`` and doubles until it encloses every root. Never returns the
    subsonic root; an area ratio below ``F(lower)`` has no root here and
    comes back NaN.
    """
    target = np.asarray(area_ratio, dtype=np.float64)
    lo_value = max(float(lower), 1.0)
    low = np.full(target.shape, math.log(lo_value))
    high = np.full(target.shape, math.log(2.0 * lo_value + 1.0))
    for _ in range(200):
        short = area_mach(np.exp(high), gamma) < target
        if not np.any(short):
            break
        high = np.where(short, high + math.log(2.0), high)
    else:
        raise ValueError("supersonic_mach: no upper bracket below M ~ 1e60")
    for _ in range(200):
        mid = 0.5 * (low + high)
        below = area_mach(np.exp(mid), gamma) < target
        low = np.where(below, mid, low)
        high = np.where(below, high, mid)
        if np.all(high - low < 1.0e-15):
            break
    mach = np.exp(0.5 * (low + high))
    invalid = target < area_mach(lo_value, gamma) * (1.0 - 1.0e-14)
    return np.where(invalid, np.nan, mach)


# --------------------------------------------------------------------------- #
# C + D -- the equivalent source
# --------------------------------------------------------------------------- #

def mean_free_path(number_density, temperature_K, species, convention: str):
    """Vectorised :func:`plumetools.cai2012.gas.mean_free_path` [m].

    The same expression -- VHS with the local temperature, or hard sphere --
    on arrays; :mod:`tests.unit.test_cai2012ext_interface` pins it against
    the scalar original.
    """
    n = np.asarray(number_density, dtype=np.float64)
    hard = 1.0 / (math.sqrt(2.0) * math.pi * species.diameter_m ** 2 * n)
    if convention == "hard_sphere":
        return hard
    if convention != "vhs":
        raise ValueError(f"unknown mean-free-path convention {convention!r}")
    t = np.asarray(temperature_K, dtype=np.float64)
    return hard * (t / species.t_ref_K) ** (float(species.omega) - 0.5)


@dataclass(frozen=True)
class EquivalentSource:
    """The Cai exit mapped onto a Markelov angular source (CP-1 .. CP-5).

    Built by :func:`equivalent_source`. Every attribute is derived; the
    methods evaluate the finite-temperature closure.
    """

    species: object
    convention: str
    gamma: float
    gas_constant: float
    cp: float
    mass_kg: float
    diameter_m: float
    knudsen: float
    n_e: float
    T_e: float
    U_e: float
    rho_e: float
    targets: ExitTargets
    theta_L: float
    k_theta: float
    I: float
    A_M: float
    exponent_offset: float
    T_t: float
    T_a: float
    U_a: float
    M_a: float
    r_a: float
    r_star: float
    V_inf: float
    p_t_r_o2: float
    source_point: tuple = (0.0, 0.0, 0.0)
    exit_origin: tuple = (0.0, 0.0, 0.0)
    axis: tuple = (1.0, 0.0, 0.0)

    # --- the angular law, through the frozen 1999 implementation -------------
    def angular(self, theta):
        """Markelov's ``f(theta)``, :func:`plumetools.markelov1999.sourceflow.angular`."""
        return sf.angular(self.gamma, theta, self.exponent_offset)

    def j(self, theta):
        """Mass flow per steradian ``mdot f/(2 pi I)`` [kg/(s sr)] (CP-3)."""
        return self.targets.mdot_kg_per_s * self.angular(theta) / (2.0 * math.pi * self.I)

    def core_fraction(self, theta_c: float) -> float:
        """``int_0^theta_c f sin / I`` -- normalised by the FULL plume's ``I``."""
        theta_c = min(float(theta_c), self.theta_L)
        if theta_c <= 0.0:
            return 0.0
        value, _ = integrate.quad(
            lambda t: math.sin(t) * float(self.angular(t)), 0.0, theta_c,
            limit=200, epsabs=0.0, epsrel=1.0e-12)
        return float(value) / self.I

    # --- the radial closure --------------------------------------------------
    def mach(self, radius):
        """``M(r)`` on the supersonic branch through the anchor [-]."""
        r = np.asarray(radius, dtype=np.float64)
        return supersonic_mach(area_mach(self.M_a, self.gamma) * (r / self.r_a) ** 2,
                               self.gamma)

    def state(self, radius, theta) -> dict:
        """``M, T, U, rho, n, p, lambda, c_bar`` at ``(r, theta)``; arrays."""
        r = np.asarray(radius, dtype=np.float64)
        th = np.asarray(theta, dtype=np.float64)
        mach = self.mach(r)
        temperature = self.T_t / (1.0 + 0.5 * (self.gamma - 1.0) * mach * mach)
        speed = mach * np.sqrt(self.gamma * self.gas_constant * temperature)
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = self.j(th) / (speed * r * r)
            n = rho / self.mass_kg
            lam = mean_free_path(n, temperature, self.species, self.convention)
        c_bar = np.sqrt(8.0 * self.gas_constant * temperature / math.pi)
        return {"M": mach, "T": temperature, "U": speed, "rho": rho, "n": n,
                "p": rho * self.gas_constant * temperature,
                "lambda": lam, "c_bar": c_bar}

    def bird(self, radius, theta, definition: str = "general"):
        """Bird's parameter at ``(r, theta)`` [-], signed (positive: expansion)."""
        st = self.state(radius, theta)
        return bird_parameter(st, np.asarray(radius, dtype=np.float64), definition)

    def as_dict(self) -> dict:
        return {
            "model": MODEL_NAME,
            "gamma": self.gamma, "cp_J_per_kg_K": self.cp,
            "gas_constant_J_per_kg_K": self.gas_constant,
            "mean_free_path_convention": self.convention,
            "exit": {"Kn": self.knudsen, "D_m": self.diameter_m,
                     "n_e_per_m3": self.n_e, "T_e_K": self.T_e,
                     "U_e_m_per_s": self.U_e, "rho_e_kg_per_m3": self.rho_e},
            "targets": self.targets.as_dict(),
            "markelov": {
                "theta_L_rad": self.theta_L, "theta_L_deg": math.degrees(self.theta_L),
                "k_theta": self.k_theta, "I": self.I, "A_M": self.A_M,
                "exponent_offset": self.exponent_offset,
                "p_t_r_o2_Pa_m2": self.p_t_r_o2,
                "p_t_r_o2_note": "only this product is fixed by the mass flow "
                                 "once T_t is; no unique reservoir pressure or "
                                 "orifice radius is claimed",
                "implementation": "plumetools.markelov1999.sourceflow (frozen)",
            },
            "closure": {
                "T_t_K": self.T_t, "T_a_K": self.T_a, "U_a_m_per_s": self.U_a,
                "M_a": self.M_a, "r_a_m": self.r_a, "r_a_over_D": self.r_a / self.diameter_m,
                "r_star_m": self.r_star, "V_inf_m_per_s": self.V_inf,
                "assumptions": list(ASSUMPTIONS),
            },
            "source_point_m": list(self.source_point),
            "exit_origin_m": list(self.exit_origin),
            "axis": list(self.axis),
        }


#: The coupling assumptions, by name, for every report and manifest row.
ASSUMPTIONS = (
    "CP-1 source at the exit centre (x_s = x_e) on the exit axis",
    "CP-2 conserved targets are the exit's gross outward molecular fluxes",
    "CP-3 mass flow per steradian has Markelov's f(theta) shape, full-plume I",
    "CP-4 independent isentropic radial streamtubes (rho U r^2 = j, "
    "cp T + U^2/2 = h_target)",
    "CP-5 anchor U(r_a) = U_e and rho(r_a, 0) = rho_e",
)


def equivalent_source(cfg, *, exponent_offset: float = 0.41) -> EquivalentSource:
    """Map one Cai ``case.yaml``'s exit onto the equivalent source.

    Raises:
        ValueError: if the anchor is not a supersonic, positive-temperature
            state (``T_a <= 0`` or ``M_a <= 1``), or the geometry is not the
            full, unhalved domain the mapping assumes.
    """
    if cfg.geometry.symmetry_mode != "none":
        raise ValueError(
            f"geometry.symmetry_mode is {cfg.geometry.symmetry_mode!r}; the "
            f"interface mapping is written for the full 3-D domain")
    exit_state = cai_inflow.from_config(cfg)
    geom = cai_geometry(cfg)
    species = exit_state.species
    gamma = float(species.gamma)
    big_r = float(species.gas_constant_j_per_kg_k)
    cp = gamma * big_r / (gamma - 1.0)
    m = float(species.mass_kg)
    n_e = float(exit_state.number_density_per_m3)
    t_e = float(exit_state.T0_K)
    u_e = float(exit_state.velocity_m_per_s)
    rho_e = m * n_e

    targets = exit_targets(number_density=n_e, velocity=u_e, temperature=t_e,
                           mass_kg=m, gas_constant=big_r, cp=cp,
                           area_m2=math.pi * geom.radius_m ** 2)

    theta_l = sf.limiting_angle(gamma)
    k_theta = (gamma + exponent_offset) / (gamma - 1.0)
    big_i = sf.normalization_integral(gamma, exponent_offset)
    a_m = sf.normalization_coefficient(gamma, exponent_offset)

    t_t = targets.h_J_per_kg / cp
    t_a = (targets.h_J_per_kg - 0.5 * u_e * u_e) / cp
    if not t_a > 0.0:
        raise ValueError(
            f"the anchor temperature T_a = {t_a:g} K is not positive: the exit "
            f"carries no thermal enthalpy beyond its drift, and no supersonic "
            f"anchor exists")
    m_a = u_e / math.sqrt(gamma * big_r * t_a)
    if not m_a > 1.0:
        raise ValueError(
            f"the anchor Mach number M_a = {m_a:g} is not supersonic; one-way "
            f"injection of an outward continuum state needs M > 1")
    r_a = math.sqrt(targets.mdot_kg_per_s / (2.0 * math.pi * big_i * rho_e * u_e))
    r_star = r_a / math.sqrt(float(area_mach(m_a, gamma)))
    v_inf = sf.limiting_velocity(gamma, t_t, m)
    product = (targets.mdot_kg_per_s * v_inf
               / (4.0 * math.pi * a_m * big_i * sf.sonic_density_ratio(gamma)))

    return EquivalentSource(
        species=species, convention=exit_state.convention, gamma=gamma,
        gas_constant=big_r, cp=cp, mass_kg=m, diameter_m=float(geom.diameter_m),
        knudsen=float(exit_state.knudsen), n_e=n_e, T_e=t_e, U_e=u_e, rho_e=rho_e,
        targets=targets, theta_L=theta_l, k_theta=k_theta, I=big_i, A_M=a_m,
        exponent_offset=float(exponent_offset), T_t=t_t, T_a=t_a, U_a=u_e,
        M_a=m_a, r_a=r_a, r_star=r_star, V_inf=v_inf, p_t_r_o2=product,
    )


def markelov_terminal(source: EquivalentSource, radius, theta) -> dict:
    """The 1999 terminal model at ``(r, theta)``, mapped through ``p_t r_o^2``.

    Evaluated by the frozen :func:`plumetools.markelov1999.sourceflow.mass_density`
    with ``r_o = r_star`` and ``p_t = p_t r_o^2 / r_star^2``; any other
    ``r_o`` with the same product gives the same density. Returned for
    COMPARISON only: its velocity ``V_inf`` and temperature ``T_t`` are not a
    consistent local state (together they hold twice ``h_target``).
    """
    r = np.asarray(radius, dtype=np.float64)
    rho = sf.mass_density(
        theta, r, gamma=source.gamma, T0_K=source.T_t, mass_kg=source.mass_kg,
        p0_pa=source.p_t_r_o2 / source.r_star ** 2,
        orifice_radius_m=source.r_star, exponent_offset=source.exponent_offset)
    return {"rho": rho, "n": rho / source.mass_kg,
            "U": np.full(np.shape(rho), source.V_inf),
            "T": np.full(np.shape(rho), source.T_t)}


def compare_with_terminal_model(source: EquivalentSource, radii_over_D,
                                theta_deg=0.0) -> list[dict]:
    """The closure against the 1999 terminal model, ray by ray.

    Returns one row per radius: both densities, both speeds, the closure's
    temperature, and the specific total energy each pair implies. The
    closure holds ``h_target`` by construction; the terminal pair
    ``(V_inf, T_t)`` holds ``cp T_t + V_inf^2/2 = 2 h_target``.
    """
    radii = np.asarray(radii_over_D, dtype=np.float64) * source.diameter_m
    theta = np.full(radii.shape, math.radians(float(theta_deg)))
    closure = source.state(radii, theta)
    terminal = markelov_terminal(source, radii, theta)
    rows = []
    for i, r in enumerate(radii):
        rows.append({
            "r_over_D": float(r / source.diameter_m),
            "theta_deg": float(theta_deg),
            "M": float(closure["M"][i]),
            "rho_closure": float(closure["rho"][i]),
            "rho_terminal": float(terminal["rho"][i]),
            "rho_closure_over_terminal": float(closure["rho"][i] / terminal["rho"][i]),
            "U_closure": float(closure["U"][i]),
            "U_terminal": float(source.V_inf),
            "T_closure": float(closure["T"][i]),
            "T_terminal_as_injected": float(source.T_t),
            "h_closure_over_target": float(
                (source.cp * closure["T"][i] + 0.5 * closure["U"][i] ** 2)
                / source.targets.h_J_per_kg),
            "h_terminal_over_target": float(
                (source.cp * source.T_t + 0.5 * source.V_inf ** 2)
                / source.targets.h_J_per_kg),
        })
    return rows


def angular_shape_comparison(source: EquivalentSource, *, radius_over_D: float = 10.0,
                             angles_deg=(5, 10, 15, 20, 30, 45, 60, 75, 90),
                             n_theta: int = 360) -> dict:
    """Cumulative mass fraction inside ``theta``: Markelov's law against Cai's.

    Markelov's ``f`` belongs to a SONIC orifice (``theta_L = 90 deg`` at
    ``gamma = 5/3``); Cai's exit is at ``M_a ~ 2.2``, whose continuum
    Prandtl-Meyer turning limit is smaller. Cai's collisionless solution
    (:func:`plumetools.cai2012.analytical.moments`, Kn -> infinity, at
    ``radius_over_D``) is the other bound. The real plume at a given Kn is
    measured by DSMC, not by either; this says how far apart the two bounds
    are. CP-3 is not adjusted to close the gap.
    """
    from plumetools.cai2012 import analytical
    from plumetools.cai2012ext import comparison

    d = source.diameter_m
    big_r = float(radius_over_D) * d
    edges = np.linspace(0.0, 0.5 * math.pi, n_theta + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    x, lateral = big_r * np.cos(mid), big_r * np.sin(mid)
    exact = analytical.moments(x, lateral, np.zeros_like(x), 0.5 * d,
                               source.U_e / math.sqrt(2.0 * source.gas_constant * source.T_e))
    u = exact["U_sqrt_beta0"]
    u_r = np.cos(mid) * u[:, 0] + np.sin(mid) * u[:, 1]
    flux = comparison.mass_flux_per_steradian(
        exact["n_over_n0"], u_r, big_r, 0.5 * d,
        source.U_e / math.sqrt(2.0 * source.gas_constant * source.T_e))
    shell = flux * 2.0 * math.pi * (np.cos(edges[:-1]) - np.cos(edges[1:]))
    cumulative = np.concatenate([[0.0], np.cumsum(shell)]) / np.sum(shell)
    g = source.gamma
    nu = lambda m: (math.sqrt((g + 1) / (g - 1))
                    * math.atan(math.sqrt((g - 1) / (g + 1) * (m * m - 1)))
                    - math.atan(math.sqrt(m * m - 1)))
    rows = []
    for angle in angles_deg:
        t = math.radians(float(angle))
        rows.append({"theta_deg": float(angle),
                     "markelov_fraction": source.core_fraction(t),
                     "cai_collisionless_fraction": float(np.interp(t, edges, cumulative))})
    axis = []
    s0 = source.U_e / math.sqrt(2.0 * source.gas_constant * source.T_e)
    for r_over_d in (2.0, 5.0, float(radius_over_D)):
        r = r_over_d * d
        closure = float(source.state(np.array([r]), np.array([0.0]))["n"][0]) / source.n_e
        cai = float(analytical.centerline_density_ratio(np.array([r]), 0.5 * d, s0)[0])
        axis.append({"r_over_D": r_over_d, "closure_n_over_n_e": closure,
                     "cai_collisionless_n_over_n0": cai,
                     "closure_over_collisionless": closure / cai})
    return {
        "radius_over_D": float(radius_over_D),
        "rows": rows,
        "axis_density": axis,
        "markelov_limiting_angle_deg": math.degrees(source.theta_L),
        "prandtl_meyer_limit_from_anchor_deg": math.degrees(
            0.5 * math.pi * (math.sqrt((g + 1) / (g - 1)) - 1.0) - nu(source.M_a)),
        "note": "Markelov's f is the sonic-orifice law; CP-3 imposes it on a "
                "supersonic exit. Not tuned: the overlap study measures what "
                "it costs.",
    }


# --------------------------------------------------------------------------- #
# E -- Bird's parameter
# --------------------------------------------------------------------------- #

def bird_parameter(state: dict, radius, definition: str = "general"):
    """``B`` from a closure state at ``radius`` [-].

    ``general``: ``(2 lambda U/(r c_bar)) M^2/(M^2-1)`` -- the exact value of
    ``-(lambda/(rho c_bar)) U . grad(rho)`` on a radial streamtube.
    ``radial_approx``: Stewart's Eq. (3), without the acceleration factor.
    NaN at and below ``M = 1``, where the spherical derivative is singular.
    """
    if definition not in BIRD_DEFINITIONS:
        raise ValueError(f"unknown Bird definition {definition!r}; "
                         f"valid: {list(BIRD_DEFINITIONS)}")
    r = np.asarray(radius, dtype=np.float64)
    m2 = state["M"] ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        simple = 2.0 * state["lambda"] * state["U"] / (r * state["c_bar"])
        if definition == "radial_approx":
            return simple
        return np.where(m2 > 1.0, simple * m2 / (m2 - 1.0), np.nan)


def stewart_spherical(mach, gamma: float) -> dict:
    """Stewart & Lumpkin Eq. (2): ``T/T*, n/n*, r/r*`` against ``M``.

    The spherical source WITHOUT an angular law -- the independent
    regression the closure is checked against.
    """
    m = np.asarray(mach, dtype=np.float64)
    base = (gamma + 1.0) / (2.0 + (gamma - 1.0) * m * m)
    return {"T_over_T_star": base,
            "n_over_n_star": base ** (1.0 / (gamma - 1.0)),
            "r_over_r_star": np.sqrt(area_mach(m, gamma))}


def stewart_hard_sphere_density(mach, radius_m, bird: float, gamma: float,
                                diameter_m: float):
    """Stewart & Lumpkin Eq. (4): ``n = 0.5 sqrt(g/pi) M / (B d^2 r)``.

    The hard-sphere inversion of :func:`bird_parameter`'s ``radial_approx``
    form. Hard sphere only -- the VHS production convention has a
    temperature factor this expression does not.
    """
    return (0.5 * math.sqrt(gamma / math.pi) * np.asarray(mach, dtype=np.float64)
            / (float(bird) * diameter_m ** 2 * np.asarray(radius_m, dtype=np.float64)))


# --------------------------------------------------------------------------- #
# E -- choosing the sphere
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class SphereSelection:
    """The outcome of :func:`select_sphere`.

    Attributes:
        status: one of the ``STATUS_*`` values.
        reason: why, in a sentence a person can act on.
        radius_m: the chosen ``R_interface``, ``None`` unless admissible.
        limited_by: ``bird`` (the B root), ``domain`` (the box or the
            diagnostic sphere stopped it first) or ``None``.
        r_min_m, r_max_m: the searched interval.
        max_core_bird_at_r_min: ``max B`` over the core at the smallest
            admissible radius -- the lower bound on any threshold that could
            admit a sphere here.
        axis_bird_at_r_min: ``B(r_min, 0)``; independent of the core angle.
        non_monotonic: ``max B`` fell back below the target beyond the first
            crossing. Only the connected interval from ``r_min`` is used.
        scan: ``(r, max_core_B, valid)`` samples, thinned, for the report.
    """

    status: str
    reason: str
    radius_m: float | None
    limited_by: str | None
    r_min_m: float | None
    r_max_m: float | None
    bird_threshold: float
    bird_definition: str
    core_half_angle_deg: float | None
    core_fraction: float | None
    max_core_bird_at_r_min: float | None = None
    axis_bird_at_r_min: float | None = None
    max_core_bird_at_radius: float | None = None
    axis_bird_at_radius: float | None = None
    non_monotonic: bool = False
    scan: list = field(default_factory=list)

    @property
    def admissible(self) -> bool:
        return self.status == STATUS_ADMISSIBLE

    def as_dict(self, diameter_m: float, r_a: float) -> dict:
        def over(value, scale):
            return None if value is None else float(value) / scale
        return {
            "status": self.status,
            "reason": self.reason,
            "R_interface_m": self.radius_m,
            "R_interface_over_D": over(self.radius_m, diameter_m),
            "R_interface_over_r_a": over(self.radius_m, r_a),
            "limited_by": self.limited_by,
            "searched_r_over_D": [over(self.r_min_m, diameter_m),
                                  over(self.r_max_m, diameter_m)],
            "bird_threshold": self.bird_threshold,
            "bird_definition": self.bird_definition,
            "bird_threshold_note": "a study choice, not a universal "
                                   "continuum-breakdown criterion",
            "core_half_angle_deg": self.core_half_angle_deg,
            "core_mass_fraction": self.core_fraction,
            "excluded_tail_mass_fraction": None if self.core_fraction is None
            else 1.0 - self.core_fraction,
            "max_core_bird_at_r_min": self.max_core_bird_at_r_min,
            "axis_bird_at_r_min": self.axis_bird_at_r_min,
            "max_core_bird_at_R": self.max_core_bird_at_radius,
            "axis_bird_at_R": self.axis_bird_at_radius,
            "non_monotonic": bool(self.non_monotonic),
            "scan": self.scan,
        }


def domain_radius_limit(cfg, settings: InterfaceConfig) -> tuple[float, str]:
    """The largest sphere radius the domain and the diagnostic sphere allow.

    The hemisphere sits on ``x = 0``; it must stay ``domain_clearance_over_D``
    clear of ``x_max`` and the four lateral faces, and
    ``diagnostic_clearance_over_D`` inside the angular-comparison sphere.
    Returns ``(radius_m, what_limits_it)``.
    """
    geom = cai_geometry(cfg)
    d = geom.diameter_m
    box = min(geom.x_max_m, geom.y_max_m, -geom.y_min_m, geom.z_max_m,
              -geom.z_min_m) - settings.domain_clearance_over_D * d
    diagnostic = (settings.diagnostic_radius_over_D
                  - settings.diagnostic_clearance_over_D) * d
    if diagnostic < box:
        return diagnostic, "diagnostic sphere"
    return box, "domain box"


def max_core_bird(source: EquivalentSource, radius_m: float, theta_c: float,
                  definition: str, n_theta: int) -> tuple[float, bool]:
    """``(max B over 0 <= theta <= theta_c, every state valid)`` on one sphere.

    Valid means finite, ``T > 0``, ``M > 1`` (supersonic outward flow for
    one-way coupling) and ``n > 0`` across the whole core -- the edge included,
    which is where ``f`` is smallest and ``B`` largest.
    """
    theta = np.linspace(0.0, theta_c, int(n_theta))
    radius = np.full(theta.shape, float(radius_m))
    st = source.state(radius, theta)
    b = bird_parameter(st, radius, definition)
    valid = bool(np.all(np.isfinite(b)) and np.all(st["T"] > 0.0)
                 and np.all(st["M"] > 1.0) and np.all(st["n"] > 0.0)
                 and np.all(np.isfinite(st["n"])))
    return (float(np.nanmax(b)) if np.any(np.isfinite(b)) else float("nan")), valid


def select_sphere(source: EquivalentSource, settings: InterfaceConfig, *,
                  r_max_m: float, r_max_reason: str = "domain",
                  provisional_r_min_over_D: float | None = None) -> SphereSelection:
    """Choose ``R_interface``: the largest admissible connected radius.

    Args:
        source: the equivalent source.
        settings: the case's interface settings.
        r_max_m: the domain/diagnostic limit (:func:`domain_radius_limit`).
        r_max_reason: what sets ``r_max_m``, for the report.
        provisional_r_min_over_D: a development-smoke stand-in when
            ``source_validity_min_radius_over_D`` is null. Never used for a
            scientific purpose.

    The search runs outward from ``r_min = max(r_a, source validity radius)``.
    If ``max B`` over the core already exceeds the target there the answer is
    :data:`STATUS_NO_INTERFACE` -- moving downstream would only put the
    surface in gas further from equilibrium. Otherwise ``radius_scan_points``
    radii are sampled to ``r_max`` to find the first crossing (or an invalid
    state), which brackets a root of ``max B - B_target``; the samples beyond
    it are kept to detect a non-monotonic profile. A profile that never
    crosses is limited by the domain, and says so.
    """
    d = source.diameter_m
    theta_c_deg = settings.core_half_angle_deg
    common = dict(bird_threshold=float(settings.bird_threshold),
                  bird_definition=settings.bird_definition,
                  core_half_angle_deg=theta_c_deg)

    # Angle-independent first. theta = 0 belongs to every core, so if B on the
    # axis alone exceeds the target at every radius that could be searched --
    # from the anchor, or the validity radius if one is set, to the domain
    # limit -- no core angle and no larger validity radius can admit a sphere.
    floor_r = source.r_a
    if settings.source_validity_min_radius_over_D is not None:
        floor_r = max(floor_r, settings.source_validity_min_radius_over_D * d)
    if r_max_m > floor_r:
        grid = np.geomspace(floor_r, r_max_m, settings.radius_scan_points)
        axis = source.bird(grid, np.zeros_like(grid), settings.bird_definition)
        least = float(np.nanmin(axis))
        if least > float(settings.bird_threshold):
            at = float(grid[int(np.nanargmin(axis))])
            return SphereSelection(
                STATUS_NO_INTERFACE,
                f"B on the axis alone is at least {least:.4g} (at "
                f"{at / d:.3f} D) everywhere from {floor_r / d:.3f} D "
                f"(= {floor_r / source.r_a:.3f} r_a) to the {r_max_reason} limit "
                f"{r_max_m / d:.3f} D, above the target "
                f"{settings.bird_threshold:g}. theta = 0 is in every core, so no "
                f"core angle or larger validity radius can admit a sphere; the "
                f"exit state itself is already past this threshold.",
                None, None, floor_r, r_max_m, core_fraction=None,
                axis_bird_at_r_min=float(axis[0]), **common)

    if theta_c_deg is None:
        return SphereSelection(STATUS_INCOMPLETE, settings.missing_decisions()[0],
                               None, None, None, None, core_fraction=None, **common)
    theta_c = math.radians(theta_c_deg)
    if theta_c > source.theta_L:
        return SphereSelection(
            STATUS_INVALID_STATE,
            f"the core angle {theta_c_deg:g} deg exceeds Markelov's limiting "
            f"angle {math.degrees(source.theta_L):g} deg, where f = 0",
            None, None, None, None, core_fraction=None, **common)
    fraction = source.core_fraction(theta_c)
    common["core_fraction"] = fraction

    validity = settings.source_validity_min_radius_over_D
    if validity is None:
        if settings.scientific or provisional_r_min_over_D is None:
            return SphereSelection(
                STATUS_INCOMPLETE,
                "source_validity_min_radius_over_D is null; a scientific run "
                "needs an explicit, justified value",
                None, None, None, None, **common)
        validity = provisional_r_min_over_D
    r_min = max(source.r_a, float(validity) * d)
    if r_max_m <= r_min:
        return SphereSelection(
            STATUS_INSUFFICIENT_DOMAIN,
            f"the smallest admissible radius {r_min / d:.3f} D does not fit: "
            f"the {r_max_reason} limits the sphere to {r_max_m / d:.3f} D",
            None, None, r_min, r_max_m, **common)

    definition = settings.bird_definition
    target = float(settings.bird_threshold)
    n_theta = settings.theta_scan_points
    b_min, valid_min = max_core_bird(source, r_min, theta_c, definition, n_theta)
    axis_min = float(source.bird(np.array([r_min]), np.array([0.0]), definition)[0])
    common.update(max_core_bird_at_r_min=b_min, axis_bird_at_r_min=axis_min)
    if not valid_min:
        return SphereSelection(
            STATUS_INVALID_STATE,
            f"the closure state at r_min = {r_min / d:.3f} D is not finite, "
            f"positive and supersonic across the core",
            None, None, r_min, r_max_m, **common)
    if b_min > target:
        return SphereSelection(
            STATUS_NO_INTERFACE,
            f"max B over the {theta_c_deg:g} deg core is already {b_min:.4g} "
            f"(axis {axis_min:.4g}) at the smallest admissible radius "
            f"{r_min / d:.3f} D (= {r_min / source.r_a:.3f} r_a), above the "
            f"target {target:g}. B grows outward along every ray, so no sphere "
            f"from r_min outward satisfies it; a larger radius would be further "
            f"from equilibrium, not closer.",
            None, None, r_min, r_max_m, **common)

    radii = np.geomspace(r_min, r_max_m, settings.radius_scan_points)
    samples = []
    for r in radii:
        b, ok = max_core_bird(source, r, theta_c, definition, n_theta)
        samples.append((float(r), b, ok))

    crossing = next((i for i, (_, b, ok) in enumerate(samples)
                     if not ok or b > target), None)
    scan = [{"r_over_D": r / d, "max_core_B": b, "valid": ok}
            for r, b, ok in samples[::max(1, len(samples) // 25)]]

    if crossing is None:
        radius, limited = r_max_m, r_max_reason
        non_monotonic = False
    else:
        lo_r, _, _ = samples[crossing - 1]
        hi_r, hi_b, hi_ok = samples[crossing]
        if not hi_ok:
            # An invalid state ends the admissible interval; no root to refine.
            radius, limited = lo_r, "invalid state beyond"
        else:
            radius = optimize.brentq(
                lambda r: max_core_bird(source, r, theta_c, definition,
                                        n_theta)[0] - target,
                lo_r, hi_r, xtol=1.0e-12 * hi_r, rtol=1.0e-12, maxiter=200)
            limited = "bird"
        non_monotonic = any(ok and b <= target
                            for _, b, ok in samples[crossing + 1:])

    b_at, _ = max_core_bird(source, radius, theta_c, definition, n_theta)
    axis_at = float(source.bird(np.array([radius]), np.array([0.0]), definition)[0])
    reason = (f"max B over the {theta_c_deg:g} deg core reaches {target:g} at "
              f"{radius / d:.4f} D" if limited == "bird" else
              f"max B stays below {target:g} out to {radius / d:.4f} D, where "
              f"the {limited} stops the search; the sphere is placed there")
    return SphereSelection(
        STATUS_ADMISSIBLE, reason, float(radius), limited, r_min, r_max_m,
        max_core_bird_at_radius=b_at, axis_bird_at_radius=axis_at,
        non_monotonic=non_monotonic, scan=scan, **common)


def radius_table(source: EquivalentSource, settings: InterfaceConfig, *,
                 angles_deg, thresholds, r_min_over_D: float,
                 r_max_m: float) -> list[dict]:
    """``R_interface`` for each (core angle, threshold), for the decision.

    Purely informational -- what each choice WOULD give -- so that the angle
    and the threshold can be chosen before any DSMC result exists, and not
    against one.
    """
    rows = []
    for angle in angles_deg:
        for threshold in thresholds:
            trial = InterfaceConfig(
                enabled=True, purpose="development-smoke", run_mode="exploratory",
                bird_definition=settings.bird_definition,
                bird_threshold=float(threshold), core_half_angle_deg=float(angle),
                source_validity_min_radius_over_D=float(r_min_over_D),
                radius_scan_points=min(settings.radius_scan_points, 200),
                theta_scan_points=settings.theta_scan_points,
                domain_clearance_over_D=settings.domain_clearance_over_D,
                diagnostic_radius_over_D=settings.diagnostic_radius_over_D,
                diagnostic_clearance_over_D=settings.diagnostic_clearance_over_D)
            pick = select_sphere(source, trial, r_max_m=r_max_m)
            rows.append({
                "core_half_angle_deg": float(angle),
                "bird_threshold": float(threshold),
                "status": pick.status,
                "R_interface_over_D": None if pick.radius_m is None
                else pick.radius_m / source.diameter_m,
                "limited_by": pick.limited_by,
                "core_mass_fraction": pick.core_fraction,
                "max_core_bird_at_r_min": pick.max_core_bird_at_r_min,
            })
    return rows


# --------------------------------------------------------------------------- #
# F -- the inlet faces
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class FaceInflow:
    """Per-face injection state, every array in patch-face order.

    Attributes:
        centres, inward_normals, areas: the faces; ``inward_normals`` point
            INTO the DSMC domain (the negated OpenFOAM boundary normal).
        radius, theta: from the source point.
        in_core: ``theta <= theta_c`` -- the injecting faces.
        number_density: ``n`` on core faces, exactly 0 elsewhere.
        velocity, temperature: the closure state on EVERY face (positive
            ``T`` everywhere: ``plumeFieldInflow`` aborts on a zero).
        u_normal: ``U . n_in``.
    """

    centres: np.ndarray
    inward_normals: np.ndarray
    areas: np.ndarray
    radius: np.ndarray
    theta: np.ndarray
    in_core: np.ndarray
    number_density: np.ndarray
    velocity: np.ndarray
    temperature: np.ndarray
    u_normal: np.ndarray
    mach: np.ndarray

    @property
    def n_faces(self) -> int:
        return int(len(self.areas))


def face_inflow(source: EquivalentSource, settings: InterfaceConfig, *,
                centres, outward_normals, areas) -> FaceInflow:
    """Evaluate the closure on actual mesh faces, in their order.

    Args:
        centres: ``(n, 3)`` face centres [m].
        outward_normals: ``(n, 3)`` OpenFOAM boundary-face unit normals, which
            point OUT of the fluid -- toward the source, on this cavity.
        areas: ``(n,)`` face areas [m^2].

    Raises:
        ValueError: if any face is not oriented with the outward flow
            (``n_in . r_hat <= 0``), sits at or inside the sonic radius, or has
            a non-finite or non-supersonic state. Each means the patch is not
            the surface the model was chosen for.
    """
    c = np.asarray(centres, dtype=np.float64) - np.asarray(source.source_point)
    n_in = -np.asarray(outward_normals, dtype=np.float64)
    a = np.asarray(areas, dtype=np.float64)
    r = np.linalg.norm(c, axis=1)
    if np.any(r <= source.r_star):
        raise ValueError(
            f"{int(np.count_nonzero(r <= source.r_star))} inlet face(s) lie at or "
            f"inside the sonic radius {source.r_star:.4g} m")
    r_hat = c / r[:, None]
    theta = sf.off_axis_angle(c, plume_axis=source.axis)
    facing = np.einsum("ij,ij->i", n_in, r_hat)
    if np.any(facing <= 0.0):
        raise ValueError(
            f"{int(np.count_nonzero(facing <= 0.0))} inlet face(s) have an "
            f"inward normal pointing back at the source; the patch orientation "
            f"is not that of a cavity around it")

    st = source.state(r, theta)
    bad = ~(np.isfinite(st["T"]) & (st["T"] > 0.0) & (st["M"] > 1.0)
            & np.isfinite(st["U"]))
    if np.any(bad):
        raise ValueError(f"{int(np.count_nonzero(bad))} inlet face(s) have a "
                         f"non-finite, non-positive or subsonic closure state")
    velocity = st["U"][:, None] * r_hat
    core = theta <= settings.core_half_angle_rad + 1.0e-15
    if settings.tails == "inject_exploratory":
        density = np.where(theta < source.theta_L, st["n"], 0.0)
    else:
        density = np.where(core, st["n"], 0.0)
    return FaceInflow(
        centres=np.asarray(centres, dtype=np.float64), inward_normals=n_in,
        areas=a, radius=r, theta=theta, in_core=core,
        number_density=np.asarray(density, dtype=np.float64),
        velocity=velocity, temperature=np.asarray(st["T"], dtype=np.float64),
        u_normal=np.einsum("ij,ij->i", velocity, n_in), mach=st["M"])


def flux_balance(source: EquivalentSource, inflow: FaceInflow,
                 core_fraction: float) -> dict:
    """What the inlet faces inject, against what the exit delivers.

    Gross (what ``plumeFieldInflow`` injects, Bird eqn 4.22 per face) and net
    continuum (``rho u_n`` and ``rho u_n (cp T + |U|^2/2)``) mass and energy,
    each against the core share of the exit's gross target. The difference
    between the two is the return half a full Maxwellian would carry, which
    the injection omits; DSMC particles crossing back out are deleted. Nothing
    is rescaled to close either balance.
    """
    m = source.mass_kg
    n = inflow.number_density
    u_n = inflow.u_normal
    temp = inflow.temperature
    u_t2 = (np.einsum("ij,ij->i", inflow.velocity, inflow.velocity) - u_n * u_n)
    j1, j3 = one_way_moments(u_n, temp, source.gas_constant)
    back1, _ = one_way_moments(-u_n, temp, source.gas_constant)
    area = inflow.areas
    gross_mass = float(np.sum(area * m * n * j1))
    gross_energy = float(np.sum(
        area * 0.5 * m * n * (j3 + (u_t2 + 2.0 * source.gas_constant * temp) * j1)))
    rho = m * n
    speed2 = np.einsum("ij,ij->i", inflow.velocity, inflow.velocity)
    net_mass = float(np.sum(area * rho * u_n))
    net_energy = float(np.sum(area * rho * u_n * (source.cp * temp + 0.5 * speed2)))
    return_mass = float(np.sum(area * m * n * back1))
    mdot_core = core_fraction * source.targets.mdot_kg_per_s
    edot_core = core_fraction * source.targets.edot_W

    def rel(value, target):
        return float(value / target - 1.0) if target else None

    return {
        "core_fraction": float(core_fraction),
        "mdot_core_target_kg_per_s": mdot_core,
        "edot_core_target_W": edot_core,
        "gross_injected_mass_kg_per_s": gross_mass,
        "gross_injected_energy_W": gross_energy,
        "continuum_net_mass_kg_per_s": net_mass,
        "continuum_net_energy_W": net_energy,
        "maxwellian_return_mass_kg_per_s": return_mass,
        "gross_mass_residual": rel(gross_mass, mdot_core),
        "gross_energy_residual": rel(gross_energy, edot_core),
        "net_mass_residual": rel(net_mass, mdot_core),
        "net_energy_residual": rel(net_energy, edot_core),
        "gross_minus_net_mass_fraction": rel(gross_mass, net_mass),
        "injecting_faces": int(np.count_nonzero(n > 0.0)),
        "faces": inflow.n_faces,
        "injecting_area_m2": float(np.sum(area[n > 0.0])),
        "note": "gross = the half-range flux plumeFieldInflow injects; net = "
                "the continuum surrogate. Their difference is the return half "
                "the injection omits. The returning DSMC flux is measured after "
                "a run, not assumed; nothing here is rescaled.",
    }


def ideal_sphere_faces(radius_m: float, theta_max: float, n_theta: int = 360,
                       n_phi: int = 720):
    """An exact spherical cap ``0 <= theta <= theta_max`` as small flat faces.

    ``(centres, outward_normals, areas)`` with OpenFOAM's orientation (normals
    toward the source). Used to quote the flux balance before a mesh exists;
    the meshed patch is checked against the same targets afterwards.
    """
    edges_t = np.linspace(0.0, float(theta_max), n_theta + 1)
    edges_p = np.linspace(0.0, 2.0 * math.pi, n_phi + 1)
    tc = 0.5 * (edges_t[:-1] + edges_t[1:])
    pc = 0.5 * (edges_p[:-1] + edges_p[1:])
    tt, pp = np.meshgrid(tc, pc, indexing="ij")
    band = (np.cos(edges_t[:-1]) - np.cos(edges_t[1:]))[:, None]
    areas = (radius_m ** 2 * band * (2.0 * math.pi / n_phi)) * np.ones_like(tt)
    r_hat = np.stack([np.cos(tt), np.sin(tt) * np.cos(pp), np.sin(tt) * np.sin(pp)],
                     axis=-1).reshape(-1, 3)
    return radius_m * r_hat, -r_hat, areas.reshape(-1)


# --------------------------------------------------------------------------- #
# G -- local time scales and the reference cell
# --------------------------------------------------------------------------- #

def reference_state(source: EquivalentSource, radius_m: float,
                    cell_size_m: float) -> dict:
    """The declared reference cell: on the axis, just outside the sphere.

    Its centre is half a core cell beyond ``R_interface`` at ``theta = 0`` --
    the densest retained cell, where the weight puts 40 parcels and where the
    collision time is shortest. Reported by identity; occupancy elsewhere is
    measured after the run, not implied.
    """
    r = float(radius_m) + 0.5 * float(cell_size_m)
    st = source.state(np.array([r]), np.array([0.0]))
    sigma = math.sqrt(source.gas_constant * float(st["T"][0]))
    return {
        "identity": "core cell on the +x axis whose centre is half a core cell "
                    "beyond R_interface",
        "centre_m": [r, 0.0, 0.0],
        "number_density_per_m3": float(st["n"][0]),
        "temperature_K": float(st["T"][0]),
        "speed_m_per_s": float(st["U"][0]),
        "mach": float(st["M"][0]),
        "mean_free_path_m": float(st["lambda"][0]),
        "mean_thermal_speed_m_per_s": float(st["c_bar"][0]),
        "collision_time_s": float(st["lambda"][0] / st["c_bar"][0]),
        "characteristic_speed_m_per_s": float(st["U"][0]) + 3.0 * sigma,
    }
