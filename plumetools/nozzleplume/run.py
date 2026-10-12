r"""Time step, particle weight and run duration.

Generalised from :func:`plumetools.cai2012.inflow.derive_run_settings`, without
the sweep multipliers that family layers on top -- a single case has nothing to
sweep -- and with one emphasis reversed.

The particle weight is an input here
------------------------------------
Standard ``dsmcFoam`` has **one** particle weight for the whole domain
(``DSMCCloud::nParticle_`` is a single scalar). The Cai family derives it so the
exit cell -- the densest in the domain -- holds a target number of parcels, and
never needs to do anything else, because its densest case is still rarefied.

A real nozzle exit is not. At an exit Knudsen number of order :math:`10^{-4}`
the exit cell that resolves the nozzle holds ~:math:`10^{17}` molecules, and the
weight that puts twenty parcels in it puts tens of millions in the domain. So
``dsmc.n_equivalent_particles`` is written **exactly as given** when it is set,
and what this module does instead of deriving it is say what it costs:

.. code-block:: text

    exit occupancy   =  n0 * V_exit / weight           parcels in an exit cell
    injection        =  flux * A_nozzle * dt / weight  parcels per step
    population       ~  injection rate * residence     parcels in the domain

``null`` still derives the weight from ``dsmc.target_particles_per_exit_cell``.

The parcel estimate is rough, and says so
-----------------------------------------
The free jet is simple: gas takes ``L / U`` to reach whatever stops it, so the
jet holds that many seconds of the injection rate. The gas held up at a plate is
not: it leaves the wall thermalised to the wall temperature -- an order of
magnitude slower than it arrived -- and has to travel to an edge to escape. The
estimate charges it half the plate's shorter side at the wall's mean thermal
speed. That is an order-of-magnitude figure for choosing a weight before a run,
not a prediction; the solver log's parcel count is the measurement.

The write schedule is counted in steps
--------------------------------------
``writeControl timeStep``, with the end time rounded up to a whole number of
write intervals, for the reason :mod:`plumetools.cai2012.inflow` documents:
``runTime`` measures its schedule from the start of the *current* run, so a
resumed leg shorter than one interval writes nothing. The step index is global
and survives a restart, and a run this expensive will be restarted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class RunSettings:
    """Time step, particle weight and run duration.

    Attributes:
        delta_t_s: the time step written into ``controlDict``.
        delta_t_source: ``"derived"`` or ``"case.yaml"``.
        courant: fraction of the smallest cell a molecule at ``U + 3 sigma``
            crosses in one step.
        characteristic_speed_m_per_s: that speed.
        collision_time_s: the exit mean collision time.
        n_equivalent_particles: real molecules per simulated parcel.
        weight_source: ``"derived"`` or ``"case.yaml"``.
        exit_particles_per_cell: parcels in an exit cell at that weight.
        transit_s: domain length over the exit speed.
        average_start_s: when ``fieldAverage`` starts accumulating.
        end_time_s: end of the run -- a whole number of steps and of write
            intervals, so there is a write exactly at the end.
        write_interval_steps: ``controlDict``'s ``writeInterval``, in steps.
        n_steps_total: steps from 0 to ``end_time_s``.
        injection_rate_per_s: molecules per second through the meshed nozzle
            patch, in the modelled domain.
        flow_time_s: how long the jet takes to reach the plate, or the far
            boundary without one.
        plate_residence_s: rough time gas is held up at the plate before it
            escapes round an edge; 0 with no plate.
    """

    delta_t_s: float
    delta_t_source: str
    courant: float
    characteristic_speed_m_per_s: float
    collision_time_s: float
    n_equivalent_particles: float
    weight_source: str
    exit_particles_per_cell: float
    transit_s: float
    average_start_s: float
    end_time_s: float
    write_interval_steps: int
    n_steps_total: int
    injection_rate_per_s: float
    flow_time_s: float
    plate_residence_s: float

    # --- time ----------------------------------------------------------------
    @property
    def n_steps(self) -> int:
        """Time steps the run takes. Exact: the end time is a whole number."""
        return int(self.n_steps_total)

    @property
    def write_interval_s(self) -> float:
        """The write interval in seconds, for reporting."""
        return self.write_interval_steps * self.delta_t_s

    @property
    def n_writes(self) -> int:
        """Time directories the run produces, the last at ``endTime``."""
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
        first = int(math.ceil(self.average_start_s / self.write_interval_s - 1.0e-9))
        return max(0, self.n_writes - max(1, first) + 1)

    @property
    def time_step_over_collision_time(self) -> float:
        """``deltaT`` over the exit mean collision time. DSMC wants under 1."""
        return self.delta_t_s / self.collision_time_s

    # --- parcels -------------------------------------------------------------
    @property
    def parcels_per_step(self) -> float:
        """Parcels injected per time step, on average."""
        return self.injection_rate_per_s * self.delta_t_s / self.n_equivalent_particles

    @property
    def jet_parcels(self) -> float:
        """Parcels in the free jet at steady state: injection times flight time."""
        return self.injection_rate_per_s * self.flow_time_s / self.n_equivalent_particles

    @property
    def plate_parcels(self) -> float:
        """Rough count of parcels held up at the plate. See the module docstring."""
        return (self.injection_rate_per_s * self.plate_residence_s
                / self.n_equivalent_particles)

    @property
    def estimated_parcels(self) -> float:
        """Order-of-magnitude steady-state parcel count in the domain."""
        return self.jet_parcels + self.plate_parcels

    @property
    def settling_time_s(self) -> float:
        """Rough time for the flow to become steady [s]: flight plus hold-up."""
        return self.flow_time_s + self.plate_residence_s

    def as_dict(self) -> dict:
        return {
            "delta_t_s": self.delta_t_s,
            "delta_t_source": self.delta_t_source,
            "courant": self.courant,
            "characteristic_speed_m_per_s": self.characteristic_speed_m_per_s,
            "exit_collision_time_s": self.collision_time_s,
            "time_step_over_collision_time": self.time_step_over_collision_time,
            "n_equivalent_particles": self.n_equivalent_particles,
            "weight_source": self.weight_source,
            "exit_particles_per_cell": self.exit_particles_per_cell,
            "domain_transit_s": self.transit_s,
            "average_start_s": self.average_start_s,
            "end_time_s": self.end_time_s,
            "sampling_time_s": self.sampling_time_s,
            "write_interval_s": self.write_interval_s,
            "write_interval_steps": self.write_interval_steps,
            "n_steps": self.n_steps,
            "n_writes": self.n_writes,
            "n_sampled_writes": self.n_sampled_writes,
            "injection_rate_per_s": self.injection_rate_per_s,
            "parcels_per_step": self.parcels_per_step,
            "flow_time_s": self.flow_time_s,
            "plate_residence_s": self.plate_residence_s,
            "settling_time_s": self.settling_time_s,
            "estimated_jet_parcels": self.jet_parcels,
            "estimated_plate_parcels": self.plate_parcels,
            "estimated_parcels": self.estimated_parcels,
        }

    def describe(self) -> list[str]:
        weight_key = ("case.yaml: dsmc.n_equivalent_particles"
                      if self.weight_source == "case.yaml"
                      else "derived from dsmc.target_particles_per_exit_cell")
        lines = [
            "Run settings",
            f"  deltaT                   {self.delta_t_s:.6e} s"
            f"   ({self.delta_t_source}) -> Courant {self.courant:.3f} at"
            f" U + 3 sigma = {self.characteristic_speed_m_per_s:.0f} m/s",
            f"                           = {self.time_step_over_collision_time:.3g}"
            f" exit mean collision times ({self.collision_time_s:.3e} s)",
            f"  nEquivalentParticles     {self.n_equivalent_particles:.6e}"
            f"   ({weight_key})",
            f"  exit cell occupancy      {self.exit_particles_per_cell:.3g}"
            f" parcels/cell",
            f"  injection                {self.parcels_per_step:.4g} parcels/step"
            f"   ({self.injection_rate_per_s:.4e} molecules/s)",
            f"  domain transit           {self.transit_s:.6e} s   (x_max / U)",
            f"  sampling starts at       {self.average_start_s:.6e} s"
            f"   = {self.average_start_s / self.transit_s:.3g} transits",
            f"  end time                 {self.end_time_s:.6e} s"
            f"   = {self.n_steps:,} steps",
            f"  sampling window          {self.sampling_time_s:.6e} s"
            f"   ({self.n_sampled_writes} of {self.n_writes} writes)",
            f"  writes                   every {self.write_interval_steps:,} steps"
            f" ({self.write_interval_s:.6e} s), {self.n_writes} in total,"
            f" the last AT endTime",
            f"  estimated parcels        {self.estimated_parcels:.3g}"
            f" at steady state   [ROUGH]",
            f"                           {self.jet_parcels:.3g} in the free jet"
            f" ({self.flow_time_s:.3e} s of injection)",
        ]
        if self.plate_residence_s > 0.0:
            lines += [
                f"                           {self.plate_parcels:.3g} held up at the"
                f" plate ({self.plate_residence_s:.3e} s: half its shorter side at",
                "                           the wall's mean thermal speed --"
                " an order of magnitude, not a prediction)",
                f"  settling time            ~{self.settling_time_s:.3e} s"
                f"   = {self.settling_time_s / self.transit_s:.2g} transits before"
                f" the plate flow is steady   [ROUGH]",
            ]
        return lines


def _floor_significant(value: float, digits: int) -> float:
    """``value`` rounded down to ``digits`` significant figures.

    The 1e-9 keeps a value that already has that few figures from losing its
    last one to floating-point multiplication.
    """
    scale = 10.0 ** (int(digits) - 1 - math.floor(math.log10(value)))
    return math.floor(value * scale + 1.0e-9) / scale


def derive(cfg, exit_state, geom, plan, nozzle) -> RunSettings:
    """Derive every solver time, and resolve the particle weight.

    Args:
        cfg: the case config.
        exit_state: the :class:`~plumetools.nozzleplume.exitstate.ExitState`.
        geom: the :class:`~plumetools.nozzleplume.geometry.PlumeGeometry`.
        plan: the :class:`~plumetools.nozzleplume.mesh.MeshPlan`, for the
            smallest cell (the time step) and the exit cell (the occupancy).
        nozzle: the :class:`~plumetools.nozzleplume.mesh.NozzlePrediction`, for
            the meshed exit area the solver actually injects through.

    Returns:
        A :class:`RunSettings`.

    Raises:
        ValueError: if the run would end before averaging starts.
    """
    dsmc = cfg.dsmc
    speed = exit_state.characteristic_speed_m_per_s
    smallest = plan.min_cell_size_m

    # --- time step -----------------------------------------------------------
    if dsmc.delta_t_s is not None:
        delta_t, delta_t_source = float(dsmc.delta_t_s), "case.yaml"
    else:
        # Rounded DOWN to three figures. controlDict then carries exactly the
        # step this module counted with -- so endTime really is a whole number
        # of them -- and the Courant number can only come out under the target.
        delta_t = _floor_significant(
            float(dsmc.courant_target) * smallest / speed, 3)
        delta_t_source = "derived"
    courant = speed * delta_t / smallest

    # --- particle weight -----------------------------------------------------
    exit_molecules = exit_state.number_density_per_m3 * plan.exit_cell_volume_m3
    if dsmc.n_equivalent_particles is not None:
        weight, weight_source = float(dsmc.n_equivalent_particles), "case.yaml"
    else:
        weight = exit_molecules / float(dsmc.target_particles_per_exit_cell)
        weight_source = "derived"

    # --- durations -----------------------------------------------------------
    transit = (geom.x_max_m - geom.x_min_m) / exit_state.speed_m_per_s
    if dsmc.average_start_s is not None:
        average_start = float(dsmc.average_start_s)
    else:
        average_start = float(dsmc.transient_domain_transits) * transit
    if dsmc.end_time_s is not None:
        end_time = float(dsmc.end_time_s)
    else:
        end_time = average_start + float(dsmc.sampling_domain_transits) * transit

    if not end_time > average_start:
        raise ValueError(
            f"the run ends at {end_time:g} s but averaging would start at "
            f"{average_start:g} s, so nothing would be averaged. Raise "
            f"dsmc.end_time_s or dsmc.sampling_domain_transits, or lower "
            f"dsmc.average_start_s / dsmc.transient_domain_transits.")

    # --- the write schedule, counted in STEPS -------------------------------
    n_steps_raw = max(1, int(math.ceil(end_time / delta_t - 1.0e-9)))
    if dsmc.write_interval_s is not None:
        write_steps = max(1, int(round(float(dsmc.write_interval_s) / delta_t)))
    else:
        # At least two writes inside the sampling window: one is enough to
        # post-process, two show whether the average has settled.
        target_writes = max(2, int(math.ceil(
            2.0 * end_time / (end_time - average_start))))
        write_steps = max(1, int(math.ceil(n_steps_raw / target_writes)))
    n_writes = max(1, int(math.ceil(n_steps_raw / write_steps)))
    n_steps = write_steps * n_writes

    # --- parcels -------------------------------------------------------------
    rate = exit_state.injection_rate_per_s(nozzle.area_m2)
    residence = 0.0
    if geom.plate is not None:
        escape = 0.5 * min(geom.plate.size_y_m, geom.plate.size_z_m)
        residence = escape / exit_state.species.mean_thermal_speed(
            geom.plate.temperature_K)

    return RunSettings(
        delta_t_s=delta_t,
        delta_t_source=delta_t_source,
        courant=courant,
        characteristic_speed_m_per_s=speed,
        collision_time_s=exit_state.collision_time_s,
        n_equivalent_particles=weight,
        weight_source=weight_source,
        exit_particles_per_cell=exit_molecules / weight,
        transit_s=transit,
        average_start_s=average_start,
        end_time_s=n_steps * delta_t,
        write_interval_steps=write_steps,
        n_steps_total=n_steps,
        injection_rate_per_s=rate,
        flow_time_s=geom.flow_length_m / exit_state.speed_m_per_s,
        plate_residence_s=residence,
    )
