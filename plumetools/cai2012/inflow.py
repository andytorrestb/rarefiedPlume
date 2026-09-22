r"""The derived exit state, and the run settings that follow from it.

Nothing in this family types a number density, a velocity, a particle weight or a
time step into a configuration file. All four are derived here, from the three
values Cai actually states -- :math:`\mathrm{Kn}`, :math:`S_0` and :math:`D` --
plus the one documented assumption, :math:`T_0`:

.. code-block:: text

    Kn, D          -> lambda0 = Kn D            gas.mean_free_path_from_kn
    lambda0, T0    -> n0                        gas.number_density_from_kn   (G2)
    S0, T0         -> U0 = S0 sqrt(2 R T0)      gas.speed_from_speed_ratio
    n0, U0, T0     -> the injected flux         gas.maxwellian_number_flux

    n0, cell       -> particle weight           RunSettings
    cell, U0, T0   -> deltaT                    RunSettings

That is the whole chain, and it is why a case labelled ``Kn0p1`` cannot run
another case's density: the density does not exist anywhere until the Knudsen
number produces it.

Cai's resolution statements
---------------------------
Cai fixes :math:`\Delta x/\lambda_0 = 1` and :math:`\Delta t/t_0 = 1`, both
referred to the **Kn = 0.01 exit properties** whatever case is running, and runs
at least :math:`10^4 t_0` before sampling. Those reference quantities are
computed here and reported for every case -- but the time step is **not** set to
:math:`t_0`.

The reason is arithmetic. With argon at :math:`T_0 = 300` K,
:math:`\lambda_{ref} = 2` mm and :math:`t_0 = \lambda_{ref}/\bar c = 5.0\;\mu s`,
a molecule at :math:`U_0 + 3\sigma \approx 1460` m/s covers 7.3 mm in one
:math:`t_0` -- 3.6 cells of the size Cai's own criterion asks for. A DSMC step
that lets molecules skip cells without being offered a collision partner in them
is not a valid step, and :mod:`plumetools.cai2012.checks` fails the case on it.

Since Cai does not define :math:`t_0` precisely enough to reconstruct his step,
:class:`RunSettings` derives :math:`\Delta t` from an explicit Courant target and
reports Cai's value and the ratio alongside it. That is the honest position: the
criterion is represented, not claimed to be reproduced.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from plumetools.cai2012 import gas as gas_module
from plumetools.cai2012.gas import VhsSpecies


@dataclass(frozen=True)
class ExitState:
    """The uniform drifting Maxwellian on the nozzle exit disk.

    Every attribute except ``knudsen``, ``speed_ratio``, ``T0_K`` and the
    species is ``[DERIVED]``.

    Attributes:
        species: the VHS species.
        knudsen: ``Kn``. **[PAPER]**
        speed_ratio: ``S0``. **[PAPER]**
        T0_K: exit temperature [K]. **[ASSUMPTION]**
        characteristic_length_m: the length ``Kn`` is built on [m].
        mean_free_path_m: ``lambda0 = Kn L``.
        number_density_per_m3: ``n0``.
        velocity_m_per_s: ``U0``, along ``+x``.
        beta0: ``1 / (2 R T0)`` [s^2/m^2].
        convention: which mean-free-path convention produced ``n0``.
    """

    species: VhsSpecies
    knudsen: float
    speed_ratio: float
    T0_K: float
    characteristic_length_m: float
    mean_free_path_m: float
    number_density_per_m3: float
    velocity_m_per_s: float
    beta0: float
    convention: str = "vhs"

    @property
    def most_probable_speed_m_per_s(self) -> float:
        """``sqrt(2 R T0)`` [m/s]; ``1/sqrt(beta0)``."""
        return self.species.most_probable_speed(self.T0_K)

    @property
    def mean_thermal_speed_m_per_s(self) -> float:
        """``sqrt(8 R T0 / pi)`` [m/s]."""
        return self.species.mean_thermal_speed(self.T0_K)

    @property
    def thermal_sigma_m_per_s(self) -> float:
        """Per-component thermal standard deviation ``sqrt(R T0)`` [m/s]."""
        return math.sqrt(self.species.gas_constant_j_per_kg_k * self.T0_K)

    @property
    def characteristic_speed_m_per_s(self) -> float:
        """``U0 + 3 sigma`` [m/s] -- what the time step must not outrun.

        Three standard deviations rather than the mean: the mean would let the
        fast tail cross several cells per step unnoticed, and it is the tail that
        breaks the collision sampling.
        """
        return self.velocity_m_per_s + 3.0 * self.thermal_sigma_m_per_s

    @property
    def number_flux_per_m2_s(self) -> float:
        """One-way number flux through the exit plane [1/(m^2 s)]. Bird eq. 4.22.

        Not ``n0 U0``: that omits the thermal spread, which for ``S0 = 2`` is a
        0.04% difference and for ``S0 = 0`` is the whole quantity.
        """
        return gas_module.maxwellian_number_flux(
            self.number_density_per_m3, self.speed_ratio, self.species, self.T0_K)

    def injection_rate_per_s(self, area_m2: float) -> float:
        """Molecules per second entering through an exit of ``area_m2``.

        The theoretical value the OpenFOAM flux test measures against. Passing
        the **meshed** patch area rather than ``pi R0^2`` is deliberate: the
        staircase patch is what the solver actually injects through.
        """
        return self.number_flux_per_m2_s * float(area_m2)

    def collision_time_s(self) -> float:
        """Mean collision time at the exit state [s]. ``lambda0 / c_bar``."""
        return self.mean_free_path_m / self.mean_thermal_speed_m_per_s

    def as_dict(self) -> dict:
        """Machine-readable, for the manifest and the case summary."""
        return {
            "species": self.species.name,
            "knudsen": float(self.knudsen),
            "speed_ratio": float(self.speed_ratio),
            "T0_K": float(self.T0_K),
            "characteristic_length_m": float(self.characteristic_length_m),
            "mean_free_path_m": float(self.mean_free_path_m),
            "number_density_per_m3": float(self.number_density_per_m3),
            "velocity_m_per_s": float(self.velocity_m_per_s),
            "beta0_s2_per_m2": float(self.beta0),
            "most_probable_speed_m_per_s": float(self.most_probable_speed_m_per_s),
            "number_flux_per_m2_s": float(self.number_flux_per_m2_s),
            "collision_time_s": float(self.collision_time_s()),
            "mean_free_path_convention": self.convention,
        }

    def describe(self) -> list[str]:
        """Report lines, printed before meshing and before running."""
        return [
            "Cai 2012 exit state (DERIVED from Kn, S0, T0)",
            f"  Kn                       {self.knudsen:g}"
            f"   (on L = {self.characteristic_length_m:g} m)",
            f"  lambda0 = Kn L           {self.mean_free_path_m:.6e} m",
            f"  n0                       {self.number_density_per_m3:.6e} 1/m^3"
            f"   ({self.convention} convention)",
            f"  S0                       {self.speed_ratio:g}",
            f"  U0 = S0 sqrt(2 R T0)     {self.velocity_m_per_s:.4f} m/s",
            f"  T0                       {self.T0_K:g} K   [ASSUMPTION]",
            f"  number flux              {self.number_flux_per_m2_s:.6e} 1/(m^2 s)",
            f"  mean collision time      {self.collision_time_s():.6e} s",
        ]


def from_config(cfg) -> ExitState:
    """Derive the exit state from a validated ``case.yaml``."""
    species = cfg.species()
    length = cfg.characteristic_length_m
    convention = cfg.gas.mean_free_path_convention
    t0 = float(cfg.exit.T0_K)

    return ExitState(
        species=species,
        knudsen=float(cfg.exit.knudsen),
        speed_ratio=float(cfg.exit.speed_ratio),
        T0_K=t0,
        characteristic_length_m=length,
        mean_free_path_m=gas_module.mean_free_path_from_kn(
            cfg.exit.knudsen, length),
        number_density_per_m3=gas_module.number_density_from_kn(
            cfg.exit.knudsen, length, species, t0, convention=convention),
        velocity_m_per_s=gas_module.speed_from_speed_ratio(
            cfg.exit.speed_ratio, species, t0),
        beta0=species.beta(t0),
        convention=convention,
    )


def reference_state(cfg) -> ExitState:
    """The ``Kn = resolution.reference_knudsen`` exit state.

    Cai refers ``dx/lambda0`` and ``dt/t0`` to the **Kn = 0.01** exit properties
    for every case, not to the case's own. This builds that state so the
    reference length and time can be reported next to whatever the case actually
    uses.
    """
    species = cfg.species()
    length = cfg.characteristic_length_m
    convention = cfg.gas.mean_free_path_convention
    t0 = float(cfg.exit.T0_K)
    kn = float(cfg.resolution.reference_knudsen)

    return ExitState(
        species=species,
        knudsen=kn,
        speed_ratio=float(cfg.exit.speed_ratio),
        T0_K=t0,
        characteristic_length_m=length,
        mean_free_path_m=gas_module.mean_free_path_from_kn(kn, length),
        number_density_per_m3=gas_module.number_density_from_kn(
            kn, length, species, t0, convention=convention),
        velocity_m_per_s=gas_module.speed_from_speed_ratio(
            cfg.exit.speed_ratio, species, t0),
        beta0=species.beta(t0),
        convention=convention,
    )


# --------------------------------------------------------------------------- #
# run settings
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class RunSettings:
    """Time step, particle weight and run duration, all derived.

    Attributes:
        delta_t_s: the time step written into ``controlDict``.
        delta_t_source: ``"derived"`` or ``"case.yaml"``.
        courant: fraction of the smallest cell a molecule at
            ``U0 + 3 sigma`` crosses in one step.
        n_equivalent_particles: molecules per simulated parcel, **after** the
            numerical-particle multiplier has been applied.
        baseline_n_equivalent_particles: the weight before that division -- the
            value a 1x case runs at, and the value every other level is defined
            against.
        numerical_particle_multiplier: the factor the parcel population was
            raised by. ``nEquivalentParticles = baseline / multiplier``.
        weight_source: ``"derived"`` or ``"case.yaml"``.
        exit_particles_per_cell: occupancy in the exit cell, by construction the
            resolution target when the weight is derived.
        average_start_s: when ``fieldAverage`` starts accumulating. **Not**
            moved by ``run_time_multiplier``: the transient is a physical
            statement about when the plume is established, and a longer run is
            for collecting more samples after it, not a longer startup.
        end_time_s: end of the run. Always a whole number of time steps, and a
            whole number of write intervals, so there is a write exactly at the
            end -- see :func:`derive_run_settings`.
        baseline_end_time_s: the unscaled end time, before
            ``run_time_multiplier``. Not snapped to the write schedule; it is
            the number the multiplier is defined against.
        run_time_multiplier: the requested scale on the end time.
        sampling_time_s: ``end_time_s - average_start_s`` -- the window
            ``fieldAverage`` actually accumulates over.
        write_interval_steps: ``controlDict``'s ``writeInterval``, in **steps**.
            ``writeControl timeStep`` rather than ``runTime`` because the step
            index is global and survives a restart, where the ``runTime``
            schedule restarts from the resume point.
        write_interval_s: the same interval in seconds, for reporting.
        baseline_write_interval_steps: the interval before
            ``output_frequency_multiplier``, in steps.
        output_frequency_multiplier: the requested increase in write frequency.
        achieved_output_frequency_multiplier: what the integer step interval
            actually delivers -- ``baseline / write_interval_steps``. Reported
            rather than assumed, because ``writeControl timeStep`` cannot
            express a fractional interval.
        n_steps_total: steps from 0 to ``end_time_s``.
        transient_basis: which basis set ``average_start_s``.
        transient_cai_s: ``transient_collision_times * t0_ref`` -- Cai's own
            transient, always reported whether or not it was used.
        transient_transits_s: the domain-transit alternative.
        domain_transit_s: domain length over ``U0``.
        reference_mean_free_path_m: ``lambda_ref`` at ``Kn = 0.01``.
        reference_cell_size_m: Cai's ``dx = target_cell_over_mfp * lambda_ref``.
        reference_collision_time_s: ``t0`` at the reference state.
        cai_delta_t_s: the step Cai's ``dt/t0 = 1`` implies.
        cai_courant: the Courant number **that** step would give here.
    """

    delta_t_s: float
    delta_t_source: str
    courant: float
    n_equivalent_particles: float
    weight_source: str
    exit_particles_per_cell: float
    average_start_s: float
    end_time_s: float
    write_interval_s: float
    write_interval_steps: int
    n_steps_total: int
    baseline_n_equivalent_particles: float
    numerical_particle_multiplier: float
    baseline_end_time_s: float
    run_time_multiplier: float
    baseline_write_interval_steps: int
    output_frequency_multiplier: float
    transient_basis: str
    transient_cai_s: float
    transient_transits_s: float
    domain_transit_s: float
    reference_mean_free_path_m: float
    reference_cell_size_m: float
    reference_collision_time_s: float
    cai_delta_t_s: float
    cai_courant: float

    @property
    def n_steps(self) -> int:
        """Number of time steps the run takes. Exact: the end time is snapped
        to a whole number of steps, and to a whole number of write intervals."""
        return int(self.n_steps_total)

    @property
    def n_writes(self) -> int:
        """Number of time directories the run will produce, the last at endTime."""
        return int(self.n_steps_total // self.write_interval_steps)

    @property
    def sampling_time_s(self) -> float:
        """The window ``fieldAverage`` accumulates over [s]."""
        return self.end_time_s - self.average_start_s

    @property
    def n_sampled_writes(self) -> int:
        """Time directories written at or after ``average_start_s``.

        The ones that carry a mean field. A write before the averaging start is
        a transient snapshot and has no ``rhoNMean`` in it.
        """
        first = int(math.ceil(self.average_start_s / self.write_interval_s))
        return max(0, self.n_writes - first + 1)

    @property
    def achieved_output_frequency_multiplier(self) -> float:
        """How much more often this writes than the baseline schedule.

        ``writeControl timeStep`` takes an integer number of steps, so the
        requested multiplier is met to within that rounding. This is what was
        actually achieved, and it is what the manifest records.
        """
        return self.baseline_write_interval_steps / self.write_interval_steps

    def as_dict(self) -> dict:
        return {
            "delta_t_s": self.delta_t_s,
            "delta_t_source": self.delta_t_source,
            "courant": self.courant,
            "n_equivalent_particles": self.n_equivalent_particles,
            "weight_source": self.weight_source,
            "exit_particles_per_cell": self.exit_particles_per_cell,
            "average_start_s": self.average_start_s,
            "end_time_s": self.end_time_s,
            "write_interval_s": self.write_interval_s,
            "write_interval_steps": self.write_interval_steps,
            "n_steps": self.n_steps,
            "n_writes": self.n_writes,
            "n_sampled_writes": self.n_sampled_writes,
            "sampling_time_s": self.sampling_time_s,
            "baseline_n_equivalent_particles":
                self.baseline_n_equivalent_particles,
            "numerical_particle_multiplier": self.numerical_particle_multiplier,
            "baseline_end_time_s": self.baseline_end_time_s,
            "run_time_multiplier": self.run_time_multiplier,
            "baseline_write_interval_steps": self.baseline_write_interval_steps,
            "output_frequency_multiplier": self.output_frequency_multiplier,
            "achieved_output_frequency_multiplier":
                self.achieved_output_frequency_multiplier,
            "transient_basis": self.transient_basis,
            "transient_cai_s": self.transient_cai_s,
            "transient_transits_s": self.transient_transits_s,
            "domain_transit_s": self.domain_transit_s,
            "reference_mean_free_path_m": self.reference_mean_free_path_m,
            "reference_cell_size_m": self.reference_cell_size_m,
            "reference_collision_time_s": self.reference_collision_time_s,
            "cai_delta_t_s": self.cai_delta_t_s,
            "cai_courant": self.cai_courant,
        }

    def describe(self) -> list[str]:
        lines = [
            "Cai resolution references (Kn = 0.01 exit properties) [PAPER]",
            f"  lambda_reference         {self.reference_mean_free_path_m:.6e} m",
            f"  cell_size_reference      {self.reference_cell_size_m:.6e} m"
            f"   (Cai: dx / lambda0 = 1)",
            f"  collision_time_reference {self.reference_collision_time_s:.6e} s",
            f"  Cai deltaT (dt/t0 = 1)   {self.cai_delta_t_s:.6e} s"
            f"   -> Courant {self.cai_courant:.2f} on this mesh",
            "",
            "Run settings (DERIVED)",
            f"  deltaT                   {self.delta_t_s:.6e} s"
            f"   ({self.delta_t_source}) -> Courant {self.courant:.3f}",
            f"  particle weight          {self.n_equivalent_particles:.6e}"
            f"   ({self.weight_source})",
            f"  exit cell occupancy      {self.exit_particles_per_cell:.1f}"
            f" particles/cell",
            f"  domain transit           {self.domain_transit_s:.6e} s",
            f"  transient (Cai, 10^4 t0) {self.transient_cai_s:.6e} s",
            f"  transient (transits)     {self.transient_transits_s:.6e} s",
            f"  sampling starts at       {self.average_start_s:.6e} s"
            f"   (basis: {self.transient_basis})",
            f"  end time                 {self.end_time_s:.6e} s"
            f"   = {self.n_steps} steps",
            f"  sampling window          {self.sampling_time_s:.6e} s"
            f"   ({self.n_sampled_writes} of {self.n_writes} writes)",
            f"  writes                   every {self.write_interval_steps} steps"
            f" ({self.write_interval_s:.6e} s), {self.n_writes} in total,"
            f" the last AT endTime",
        ]
        if self.numerical_particle_multiplier != 1.0:
            lines.append(
                f"  numerical particles      {self.numerical_particle_multiplier:g}x"
                f" the baseline: weight"
                f" {self.baseline_n_equivalent_particles:.6e} /"
                f" {self.numerical_particle_multiplier:g}")
        if self.run_time_multiplier != 1.0:
            lines.append(
                f"  run time                 {self.run_time_multiplier:g}x the"
                f" baseline {self.baseline_end_time_s:.6e} s; the transient is"
                f" unchanged, so all of it is sampling")
        if self.output_frequency_multiplier != 1.0:
            lines.append(
                f"  output frequency         {self.output_frequency_multiplier:g}x"
                f" requested, {self.achieved_output_frequency_multiplier:.3f}x"
                f" achieved ({self.baseline_write_interval_steps} ->"
                f" {self.write_interval_steps} steps; writeControl timeStep takes"
                f" an integer)")
        if self.transient_basis != "cai":
            lines.append(
                f"  NOTE: sampling starts at {self.average_start_s / self.reference_collision_time_s:.0f}"
                f" reference collision times, not the 10 000 Cai states. "
                f"dsmc.transient_basis is '{self.transient_basis}'.")
        return lines


def derive_run_settings(cfg, exit_state: ExitState, geom, *,
                        min_cell_size_m: float,
                        exit_cell_volume_m3: float) -> RunSettings:
    """Derive every solver time and the particle weight.

    Args:
        cfg: the case config.
        exit_state: this case's exit state.
        geom: the :class:`~plumetools.cai2012.geometry.CaiGeometry`.
        min_cell_size_m: the smallest cell in the mesh [m]; sets the Courant
            number and hence the time step.
        exit_cell_volume_m3: volume of a cell at the exit [m^3]; sets the
            particle weight.

    Returns:
        A :class:`RunSettings`.

    Raises:
        ValueError: for a non-positive cell size or volume.
    """
    if not (min_cell_size_m > 0.0 and exit_cell_volume_m3 > 0.0):
        raise ValueError(
            f"min_cell_size_m ({min_cell_size_m}) and exit_cell_volume_m3 "
            f"({exit_cell_volume_m3}) must both be positive")

    reference = reference_state(cfg)
    reference_collision_time = reference.collision_time_s()
    reference_cell = float(cfg.mesh.target_cell_over_mfp) * reference.mean_free_path_m

    speed = exit_state.characteristic_speed_m_per_s

    # --- time step -----------------------------------------------------------
    if cfg.dsmc.delta_t_s is not None:
        delta_t = float(cfg.dsmc.delta_t_s)
        delta_t_source = "case.yaml"
    else:
        delta_t = float(cfg.dsmc.courant_target) * min_cell_size_m / speed
        delta_t_source = "derived"
    courant = speed * delta_t / min_cell_size_m
    cai_courant = speed * reference_collision_time / min_cell_size_m

    # --- particle weight -----------------------------------------------------
    #
    # The BASELINE weight first, exactly as it has always been derived, and then
    # the numerical-particle multiplier applied to it. One parcel stands for
    # nEquivalentParticles molecules, so the parcel count goes as the RECIPROCAL
    # of the weight and raising the population by F means dividing by F.
    #
    # The multiplier is applied whichever way the baseline arose. A case that
    # pins dsmc.n_equivalent_particles is still entitled to a 2x variant, and
    # defining the levels against whatever the 1x case runs at is the only way
    # the sweep stays synchronised with the study it extends.
    if cfg.dsmc.n_equivalent_particles is not None:
        baseline_weight = float(cfg.dsmc.n_equivalent_particles)
        weight_source = "case.yaml"
    else:
        baseline_weight = (exit_state.number_density_per_m3 * exit_cell_volume_m3
                           / float(cfg.resolution.target_particles_per_cell))
        weight_source = "derived"

    particle_multiplier = float(cfg.dsmc.numerical_particle_multiplier)
    if not particle_multiplier > 0.0:
        raise ValueError(
            f"dsmc.numerical_particle_multiplier is {particle_multiplier}; it "
            f"divides the particle weight, so it must be positive. 1.0 is the "
            f"baseline.")
    weight = baseline_weight / particle_multiplier
    occupancy = (exit_state.number_density_per_m3 * exit_cell_volume_m3 / weight)

    # --- durations -----------------------------------------------------------
    transit = (geom.x_max_m - geom.x_min_m) / exit_state.velocity_m_per_s
    transient_cai = float(cfg.dsmc.transient_collision_times) * reference_collision_time
    transient_transits = float(cfg.dsmc.transient_domain_transits) * transit
    basis = cfg.dsmc.transient_basis
    transient = transient_cai if basis == "cai" else transient_transits

    if cfg.dsmc.average_start_s is not None:
        average_start = float(cfg.dsmc.average_start_s)
    else:
        average_start = transient

    if cfg.dsmc.end_time_s is not None:
        baseline_end_time = float(cfg.dsmc.end_time_s)
    else:
        baseline_end_time = (average_start
                             + float(cfg.dsmc.sampling_domain_transits) * transit)

    # The run-time multiplier scales the END, not the transient. average_start
    # is a physical statement -- when the plume is established -- and moving it
    # with the run length would buy a longer startup instead of more samples.
    # Every second the multiplier adds therefore lands inside the averaging
    # window: at 3x here the run is three times as long and the sampled window
    # is (3 E - A) / (E - A) times as long, which is more than three.
    run_time_multiplier = float(cfg.dsmc.run_time_multiplier)
    if not run_time_multiplier > 0.0:
        raise ValueError(
            f"dsmc.run_time_multiplier is {run_time_multiplier}; it scales the "
            f"end time, so it must be positive. 1.0 is the baseline.")
    end_time = run_time_multiplier * baseline_end_time

    if not end_time > average_start:
        raise ValueError(
            f"the run ends at {end_time:g} s but sampling would start at "
            f"{average_start:g} s, so nothing would be averaged. Raise "
            f"dsmc.end_time_s / dsmc.sampling_domain_transits / "
            f"dsmc.run_time_multiplier, or shorten the transient.")

    # --- the write schedule, counted in STEPS -------------------------------
    #
    # `writeControl runTime` measures its schedule from the **start time of the
    # current run**: ``Time::operator++`` computes the write index as
    # ``(value - startTime)/writeInterval``. On a RESUME that start time is the
    # resume point, so the schedule begins again there -- and a resumed leg
    # shorter than one write interval writes nothing at all. Measured: resuming
    # Kn100 at t = 8.487e-3 s with a 1.980e-3 s interval ran the remaining 1030
    # steps to endTime and produced no time directory.
    #
    # `writeControl timeStep` counts the **global** step index, which each time
    # directory stores in ``uniform/time`` and which continues across restarts.
    # Rounding the run up to a whole number of write intervals then puts a write
    # exactly at endTime, wherever the run happened to be resumed from.
    n_steps_raw = max(1, int(math.ceil(end_time / delta_t - 1.0e-9)))

    # The BASELINE write interval is derived from the BASELINE run, not from the
    # lengthened one. Deriving it from the scaled end time would make the write
    # interval a function of dsmc.run_time_multiplier, and then "four times more
    # often than before" would not be four times anything anybody can point at.
    baseline_n_steps_raw = max(1, int(math.ceil(
        baseline_end_time / delta_t - 1.0e-9)))
    if cfg.dsmc.write_interval_s is not None:
        baseline_write_steps = max(
            1, int(round(float(cfg.dsmc.write_interval_s) / delta_t)))
    else:
        # At least two writes inside the sampling window: one is enough to
        # post-process, two show whether the average has settled.
        target_writes = max(2, int(math.ceil(
            2.0 * baseline_end_time / (baseline_end_time - average_start))))
        baseline_write_steps = max(
            1, int(math.ceil(baseline_n_steps_raw / target_writes)))

    # writeControl is timeStep -- see the comment above -- so the interval is a
    # whole number of steps and the requested frequency is met to within that
    # rounding. NEAREST integer, and never below 1: a fractional interval cannot
    # be written, and 0 would mean writing nothing.
    output_multiplier = float(cfg.dsmc.output_frequency_multiplier)
    if not output_multiplier > 0.0:
        raise ValueError(
            f"dsmc.output_frequency_multiplier is {output_multiplier}; it "
            f"divides the write interval, so it must be positive. 1.0 is the "
            f"baseline.")
    write_steps = max(1, int(round(baseline_write_steps / output_multiplier)))

    n_writes = max(1, int(math.ceil(n_steps_raw / write_steps)))
    n_steps = write_steps * n_writes
    end_time = n_steps * delta_t
    write_interval = write_steps * delta_t

    return RunSettings(
        write_interval_steps=write_steps,
        n_steps_total=n_steps,
        baseline_n_equivalent_particles=baseline_weight,
        numerical_particle_multiplier=particle_multiplier,
        baseline_end_time_s=baseline_end_time,
        run_time_multiplier=run_time_multiplier,
        baseline_write_interval_steps=baseline_write_steps,
        output_frequency_multiplier=output_multiplier,
        delta_t_s=delta_t,
        delta_t_source=delta_t_source,
        courant=courant,
        n_equivalent_particles=weight,
        weight_source=weight_source,
        exit_particles_per_cell=occupancy,
        average_start_s=average_start,
        end_time_s=end_time,
        write_interval_s=write_interval,
        transient_basis=basis,
        transient_cai_s=transient_cai,
        transient_transits_s=transient_transits,
        domain_transit_s=transit,
        reference_mean_free_path_m=reference.mean_free_path_m,
        reference_cell_size_m=reference_cell,
        reference_collision_time_s=reference_collision_time,
        cai_delta_t_s=reference_collision_time,
        cai_courant=cai_courant,
    )
