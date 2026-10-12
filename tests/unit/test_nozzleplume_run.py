"""Time step, particle weight, durations, and the parcel estimate.

The weight is the centre of it. ``dsmc.n_equivalent_particles`` is the one
number a user of an expensive case edits by hand, so the properties protected
here are that it is written exactly as given, and that everything quoted about
its cost -- parcels in an exit cell, parcels per step, parcels in the domain --
follows from it in the obvious proportion.
"""

from __future__ import annotations

import math

import pytest

from _nozzleplume import R4D_CASE, derive, derived
from plumetools.nozzleplume.config import load_case_config


# --------------------------------------------------------------------------- #
# the particle weight
# --------------------------------------------------------------------------- #

def test_a_pinned_weight_is_used_exactly_as_given(tmp_path):
    *_, run = derived(tmp_path, dsmc={"n_equivalent_particles": 3.7e16})
    assert run.n_equivalent_particles == 3.7e16
    assert run.weight_source == "case.yaml"


def test_the_exit_occupancy_follows_from_the_pinned_weight(tmp_path):
    """n V / weight, with V the cube of the exit cell: 62.5 mm here."""
    *_, run = derived(tmp_path, dsmc={"n_equivalent_particles": 1.0e17})
    assert run.exit_particles_per_cell == pytest.approx(
        1.0e22 * 0.0625 ** 3 / 1.0e17)


def test_a_null_weight_is_derived_to_meet_the_occupancy_target(tmp_path):
    *_, run = derived(tmp_path, dsmc={"n_equivalent_particles": None,
                                      "target_particles_per_exit_cell": 25.0})
    assert run.weight_source == "derived"
    assert run.exit_particles_per_cell == pytest.approx(25.0)
    assert run.n_equivalent_particles == pytest.approx(1.0e22 * 0.0625 ** 3 / 25.0)


def test_the_target_is_ignored_when_the_weight_is_pinned(tmp_path):
    *_, a = derived(tmp_path / "a", dsmc={"target_particles_per_exit_cell": 5.0})
    *_, b = derived(tmp_path / "b", dsmc={"target_particles_per_exit_cell": 50.0})
    assert a.n_equivalent_particles == b.n_equivalent_particles == 1.0e17


def test_halving_the_weight_doubles_every_parcel_count(tmp_path):
    *_, heavy = derived(tmp_path / "a", dsmc={"n_equivalent_particles": 2.0e17})
    *_, light = derived(tmp_path / "b", dsmc={"n_equivalent_particles": 1.0e17})
    for name in ("exit_particles_per_cell", "parcels_per_step", "jet_parcels",
                 "plate_parcels", "estimated_parcels"):
        assert getattr(light, name) == pytest.approx(2.0 * getattr(heavy, name)), name


def test_the_weight_changes_no_time(tmp_path):
    """It is a statistical knob: the step, the end time and the write schedule
    are those of the mesh and the flow, not of the parcel count."""
    *_, heavy = derived(tmp_path / "a", dsmc={"n_equivalent_particles": 2.0e17})
    *_, light = derived(tmp_path / "b", dsmc={"n_equivalent_particles": 1.0e16})
    assert (heavy.delta_t_s, heavy.end_time_s, heavy.write_interval_steps) == (
        light.delta_t_s, light.end_time_s, light.write_interval_steps)


def test_injection_is_the_flux_through_the_meshed_nozzle(tmp_path):
    """Through the staircase patch, not the disk: that is what the solver
    injects across."""
    _, _, exit_state, _, nozzle, run = derived(tmp_path)
    assert run.injection_rate_per_s == pytest.approx(
        exit_state.number_flux_per_m2_s * nozzle.area_m2)
    assert run.parcels_per_step == pytest.approx(
        run.injection_rate_per_s * run.delta_t_s / 1.0e17)


def test_a_symmetric_domain_injects_its_share(tmp_path):
    *_, full = derived(tmp_path / "a")
    *_, quarter = derived(tmp_path / "b", domain={"symmetry_planes": ["y", "z"]})
    assert quarter.injection_rate_per_s == pytest.approx(
        full.injection_rate_per_s / 4.0)


# --------------------------------------------------------------------------- #
# the time step
# --------------------------------------------------------------------------- #

def test_the_derived_step_meets_the_courant_target_from_below(tmp_path):
    _, _, exit_state, plan, _, run = derived(tmp_path)
    exact = 0.2 * plan.min_cell_size_m / exit_state.characteristic_speed_m_per_s
    assert run.delta_t_source == "derived"
    assert run.delta_t_s <= exact
    assert run.delta_t_s == pytest.approx(exact, rel=0.01)
    assert run.courant <= 0.2


def test_the_derived_step_has_three_significant_figures(tmp_path):
    """So controlDict carries exactly the step the end time was counted in."""
    *_, run = derived(tmp_path)
    mantissa = run.delta_t_s / 10.0 ** math.floor(math.log10(run.delta_t_s))
    assert mantissa * 100.0 == pytest.approx(round(mantissa * 100.0), abs=1e-6)


def test_a_pinned_step_is_used_and_its_courant_number_reported(tmp_path):
    _, _, exit_state, plan, _, run = derived(tmp_path, dsmc={"delta_t_s": 1.0e-6})
    assert (run.delta_t_s, run.delta_t_source) == (1.0e-6, "case.yaml")
    assert run.courant == pytest.approx(
        exit_state.characteristic_speed_m_per_s * 1.0e-6 / plan.min_cell_size_m)


def test_the_step_is_compared_with_the_exit_collision_time(tmp_path):
    _, _, exit_state, _, _, run = derived(tmp_path)
    assert run.time_step_over_collision_time == pytest.approx(
        run.delta_t_s / exit_state.collision_time_s)


# --------------------------------------------------------------------------- #
# durations and the write schedule
# --------------------------------------------------------------------------- #

def test_durations_are_counted_in_transits_of_the_domain(tmp_path):
    """6 m at 3000 m/s is 2 ms; 4 transits of transient, 2 of sampling."""
    *_, run = derived(tmp_path)
    assert run.transit_s == pytest.approx(2.0e-3)
    assert run.average_start_s == pytest.approx(8.0e-3)
    assert run.end_time_s == pytest.approx(12.0e-3, rel=5e-3)


def test_the_end_time_is_whole_steps_and_whole_write_intervals(tmp_path):
    *_, run = derived(tmp_path)
    assert run.end_time_s == pytest.approx(run.n_steps * run.delta_t_s, rel=1e-12)
    assert run.n_steps % run.write_interval_steps == 0
    assert run.n_writes * run.write_interval_steps == run.n_steps


def test_the_end_time_is_never_rounded_down(tmp_path):
    """Rounding to whole write intervals may lengthen the run, not shorten it."""
    *_, run = derived(tmp_path)
    assert run.end_time_s >= 12.0e-3 * (1.0 - 1e-9)


def test_at_least_two_writes_fall_in_the_sampling_window(tmp_path):
    *_, run = derived(tmp_path)
    assert run.n_sampled_writes >= 2
    assert run.n_sampled_writes <= run.n_writes


def test_pinned_times_are_honoured(tmp_path):
    *_, run = derived(tmp_path, dsmc={
        "delta_t_s": 1.0e-6, "end_time_s": 4.0e-3, "average_start_s": 1.5e-3,
        "write_interval_s": 1.0e-3})
    assert run.average_start_s == 1.5e-3
    assert run.write_interval_steps == 1000
    assert (run.n_steps, run.n_writes) == (4000, 4)
    assert run.sampling_time_s == pytest.approx(2.5e-3)
    assert run.n_sampled_writes == 3            # the writes at 2, 3 and 4 ms


def test_averaging_from_the_start_counts_every_write(tmp_path):
    *_, run = derived(tmp_path, dsmc={
        "delta_t_s": 1.0e-6, "end_time_s": 4.0e-3, "average_start_s": 0.0,
        "write_interval_s": 1.0e-3})
    assert run.n_sampled_writes == run.n_writes == 4


def test_a_run_that_ends_before_averaging_starts_is_refused(tmp_path):
    with pytest.raises(ValueError, match="nothing would be averaged"):
        derived(tmp_path, dsmc={"end_time_s": 1.0e-3, "average_start_s": 2.0e-3})


# --------------------------------------------------------------------------- #
# the parcel estimate
# --------------------------------------------------------------------------- #

def test_the_free_jet_holds_its_flight_time_of_injection(tmp_path):
    """4 m to the plate at 3000 m/s."""
    *_, run = derived(tmp_path)
    assert run.flow_time_s == pytest.approx(4.0 / 3000.0)
    assert run.jet_parcels == pytest.approx(
        run.injection_rate_per_s * run.flow_time_s / 1.0e17)


def test_gas_is_held_up_at_the_plate_for_an_edge_crossing_at_wall_speed(tmp_path):
    """Half the plate's shorter side -- 1 m of a 4 x 2 m plate -- at the mean
    thermal speed of the 300 K wall."""
    _, _, exit_state, _, _, run = derived(tmp_path)
    assert run.plate_residence_s == pytest.approx(
        1.0 / exit_state.species.mean_thermal_speed(300.0))
    assert run.settling_time_s == pytest.approx(
        run.flow_time_s + run.plate_residence_s)


def test_without_a_plate_nothing_is_held_up(tmp_path):
    *_, run = derived(tmp_path, plate=None)
    assert run.plate_residence_s == 0.0
    assert run.plate_parcels == 0.0
    assert run.flow_time_s == pytest.approx(6.0 / 3000.0)   # the whole box


def test_the_report_says_where_the_weight_came_from(tmp_path):
    *_, pinned = derived(tmp_path / "a")
    *_, free = derived(tmp_path / "b", dsmc={"n_equivalent_particles": None})
    assert any("dsmc.n_equivalent_particles" in line for line in pinned.describe())
    assert any("target_particles_per_exit_cell" in line for line in free.describe())


# --------------------------------------------------------------------------- #
# the shipped case
# --------------------------------------------------------------------------- #

def test_r4d_runs_at_the_weight_in_its_case_yaml():
    cfg = load_case_config(R4D_CASE)
    *_, run = derive(cfg)
    assert run.n_equivalent_particles == cfg.dsmc.n_equivalent_particles
    assert run.weight_source == "case.yaml"


def test_r4d_averaging_starts_after_the_flow_has_settled():
    *_, run = derive(load_case_config(R4D_CASE))
    assert run.average_start_s > run.settling_time_s


def test_r4d_time_step_is_set_by_the_exit_cells():
    """31.25 mm cells and a 4536 m/s fast tail: 1.37 microseconds."""
    *_, run = derive(load_case_config(R4D_CASE))
    assert run.delta_t_s == pytest.approx(1.37e-6)
