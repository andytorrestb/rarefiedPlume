"""The derived exit state, the run settings, and the DSMC quality checks.

This is the traceability chain §12 of the case specification asks about: a
configuration labelled ``Kn = 0.1`` must never be able to run another case's
density, because no density is stored anywhere -- it only exists as a function of
the Knudsen number.
"""

from __future__ import annotations

import math

import pytest

from _cai2012 import derived, load_case
from plumetools.cai2012 import checks, gas, inflow, mesh
from plumetools.cai2012.geometry import from_config


# --------------------------------------------------------------------------- #
# the exit state
# --------------------------------------------------------------------------- #

def test_exit_state_derives_everything_from_kn_s0_and_T0(tmp_path):
    cfg = load_case(tmp_path)
    state = inflow.from_config(cfg)
    assert state.mean_free_path_m == pytest.approx(100.0 * 0.2)
    assert state.number_density_per_m3 == pytest.approx(
        gas.number_density_from_kn(100.0, 0.2, cfg.species(), 300.0))
    assert state.velocity_m_per_s == pytest.approx(706.95, rel=1e-4)


@pytest.mark.parametrize("kn", [100.0, 0.1, 0.01])
def test_the_knudsen_number_survives_into_the_exit_state(tmp_path, kn):
    state = inflow.from_config(load_case(tmp_path / str(kn), exit={"knudsen": kn}))
    assert state.knudsen == kn
    assert state.mean_free_path_m == pytest.approx(kn * 0.2)


def test_the_three_cases_have_three_different_exit_densities(tmp_path):
    densities = [
        inflow.from_config(load_case(tmp_path / str(kn), exit={"knudsen": kn}))
        .number_density_per_m3
        for kn in (100.0, 0.1, 0.01)]
    assert len(set(densities)) == 3


def test_the_speed_is_independent_of_the_knudsen_number(tmp_path):
    """S0 and T0 fix U0; the density is the only thing Kn changes."""
    speeds = {
        inflow.from_config(load_case(tmp_path / str(kn), exit={"knudsen": kn}))
        .velocity_m_per_s
        for kn in (100.0, 0.1, 0.01)}
    assert len(speeds) == 1


def test_beta0_matches_the_species(tmp_path):
    cfg = load_case(tmp_path)
    state = inflow.from_config(cfg)
    assert state.beta0 == pytest.approx(cfg.species().beta(300.0))
    assert state.most_probable_speed_m_per_s * math.sqrt(state.beta0) == (
        pytest.approx(1.0))


def test_characteristic_speed_is_the_drift_plus_three_sigma(tmp_path):
    """The mean would let the fast tail cross several cells per step unnoticed."""
    state = inflow.from_config(load_case(tmp_path))
    assert state.characteristic_speed_m_per_s == pytest.approx(
        state.velocity_m_per_s + 3.0 * state.thermal_sigma_m_per_s)
    assert state.characteristic_speed_m_per_s == pytest.approx(1457.0, rel=1e-3)


def test_number_flux_is_not_n_times_u(tmp_path):
    """n0 U0 omits the thermal spread. At S0 = 2 that is 0.04%; at S0 = 0 it is
    the whole quantity."""
    state = inflow.from_config(load_case(tmp_path))
    assert state.number_flux_per_m2_s != pytest.approx(
        state.number_density_per_m3 * state.velocity_m_per_s, rel=1e-6)
    assert state.number_flux_per_m2_s == pytest.approx(
        state.number_density_per_m3 * state.velocity_m_per_s, rel=1e-2)


def test_injection_rate_scales_with_the_meshed_area(tmp_path):
    state = inflow.from_config(load_case(tmp_path))
    assert state.injection_rate_per_s(2.0) == pytest.approx(
        2.0 * state.number_flux_per_m2_s)


def test_reference_state_is_always_kn_0p01(tmp_path):
    """Cai refers dx/lambda0 and dt/t0 to the Kn = 0.01 exit properties for
    EVERY case, not to the case's own."""
    reference = inflow.reference_state(load_case(tmp_path, exit={"knudsen": 100.0}))
    assert reference.knudsen == 0.01
    assert reference.mean_free_path_m == pytest.approx(0.002)


def test_exit_state_serialises_every_derived_quantity(tmp_path):
    document = inflow.from_config(load_case(tmp_path)).as_dict()
    for key in ("knudsen", "mean_free_path_m", "number_density_per_m3",
                "velocity_m_per_s", "number_flux_per_m2_s",
                "mean_free_path_convention"):
        assert key in document


def test_describe_marks_T0_as_an_assumption(tmp_path):
    lines = "\n".join(inflow.from_config(load_case(tmp_path)).describe())
    assert "[ASSUMPTION]" in lines
    assert "lambda0 = Kn L" in lines


# --------------------------------------------------------------------------- #
# run settings
# --------------------------------------------------------------------------- #

def test_the_time_step_hits_the_courant_target(tmp_path):
    _, _, _, _, run = derived(tmp_path, dsmc={"courant_target": 0.2})
    assert run.courant == pytest.approx(0.2)
    assert run.delta_t_source == "derived"


def test_an_explicit_time_step_is_used_and_labelled(tmp_path):
    _, _, _, _, run = derived(tmp_path, dsmc={"delta_t_s": 1.0e-7})
    assert run.delta_t_s == 1.0e-7
    assert run.delta_t_source == "case.yaml"


def test_the_particle_weight_meets_the_target_in_the_exit_cell(tmp_path):
    _, _, _, _, run = derived(tmp_path, resolution={
        "target_particles_per_cell": 20.0})
    assert run.exit_particles_per_cell == pytest.approx(20.0)
    assert run.weight_source == "derived"


def test_the_particle_weight_scales_with_the_exit_density(tmp_path):
    """Same mesh, ten times the density, ten times the weight -- so the parcel
    count stays put and only the physics changes."""
    _, _, _, _, thin = derived(tmp_path / "a", exit={"knudsen": 0.1})
    _, _, _, _, dense = derived(tmp_path / "b", exit={"knudsen": 0.01},
                                mesh={"core_cell_size_m": 0.05})
    assert dense.n_equivalent_particles / thin.n_equivalent_particles == (
        pytest.approx(10.0))


def test_an_explicit_weight_is_used_and_labelled(tmp_path):
    _, _, _, _, run = derived(tmp_path, dsmc={"n_equivalent_particles": 1.0e10})
    assert run.n_equivalent_particles == 1.0e10
    assert run.weight_source == "case.yaml"


def test_cais_reference_quantities_are_computed_for_every_case(tmp_path):
    _, _, _, _, run = derived(tmp_path, exit={"knudsen": 100.0})
    assert run.reference_mean_free_path_m == pytest.approx(0.002)
    assert run.reference_cell_size_m == pytest.approx(0.002)
    assert run.reference_collision_time_s == pytest.approx(5.01e-6, rel=1e-2)
    assert run.cai_delta_t_s == run.reference_collision_time_s


def test_cais_time_step_is_reported_not_adopted(tmp_path):
    """Cai's dt/t0 = 1 is computed for every case and reported next to the step
    actually used, rather than being adopted as if it were reconstructible."""
    _, _, state, plan, run = derived(tmp_path)
    assert run.delta_t_s != run.cai_delta_t_s
    assert run.cai_courant == pytest.approx(
        state.characteristic_speed_m_per_s * run.cai_delta_t_s
        / plan.min_cell_size_m)
    assert "Cai deltaT" in "\n".join(run.describe())


def test_cais_time_step_is_unusable_at_cais_own_cell_size(tmp_path):
    """This is why the step is derived instead. At dx = lambda_ref = 2 mm, one
    reference collision time carries a molecule at U0 + 3 sigma across 3.6
    cells -- a step at which the collision sampling is meaningless."""
    _, _, _, plan, run = derived(
        tmp_path, mesh={"core_cell_size_m": 0.002, "max_cells": 10 ** 12})
    assert plan.min_cell_size_m == pytest.approx(0.002)
    assert run.cai_courant == pytest.approx(3.65, rel=0.02)
    assert run.courant == pytest.approx(0.2)


def test_the_cai_transient_is_ten_thousand_reference_collision_times(tmp_path):
    _, _, _, _, run = derived(tmp_path, dsmc={"transient_basis": "cai"})
    assert run.transient_cai_s == pytest.approx(
        10000.0 * run.reference_collision_time_s)
    assert run.average_start_s == pytest.approx(run.transient_cai_s)


def test_the_transit_basis_is_shorter_and_says_so(tmp_path):
    _, _, _, _, run = derived(tmp_path, dsmc={"transient_basis": "transits"})
    assert run.average_start_s == pytest.approx(run.transient_transits_s)
    assert run.average_start_s < run.transient_cai_s
    assert "not the 10 000 Cai states" in "\n".join(run.describe())


def test_both_transients_are_always_reported(tmp_path):
    """Whichever basis was used, the other number is on the report."""
    for basis in ("cai", "transits"):
        _, _, _, _, run = derived(tmp_path / basis,
                                  dsmc={"transient_basis": basis})
        assert run.transient_cai_s > 0.0
        assert run.transient_transits_s > 0.0


def test_the_domain_transit_is_the_length_over_u0(tmp_path):
    _, geom, state, _, run = derived(tmp_path)
    assert run.domain_transit_s == pytest.approx(
        geom.x_max_m / state.velocity_m_per_s)


def test_the_run_ends_after_sampling_starts(tmp_path):
    _, _, _, _, run = derived(tmp_path)
    assert run.end_time_s > run.average_start_s
    assert run.n_steps == pytest.approx(run.end_time_s / run.delta_t_s)


def test_the_write_interval_divides_the_end_time(tmp_path):
    """OpenFOAM does NOT write at endTime as a special case -- it writes when the
    write index advances and nowhere else. An interval that does not divide the
    run leaves the final state unwritten, which makes the case look finished
    while having no record of its last stretch, and makes a resume re-run that
    stretch to no effect."""
    _, _, _, _, run = derived(tmp_path)
    assert run.n_steps % run.write_interval_steps == 0


def test_at_least_two_writes_fall_in_the_sampling_window(tmp_path):
    """One is enough to post-process; two show whether the average settled."""
    _, _, _, _, run = derived(tmp_path)
    writes = [k * run.write_interval_s for k in range(1, run.n_writes + 1)]
    sampled = [t for t in writes if t > run.average_start_s - 1e-15]
    assert len(sampled) >= 2
    assert sampled[-1] == pytest.approx(run.end_time_s)


def test_an_explicit_write_interval_is_rounded_to_whole_steps(tmp_path):
    """controlDict writes on a STEP count, so the configured seconds become the
    nearest whole number of steps -- and the run is still snapped up so the last
    write lands at endTime."""
    _, _, _, _, run = derived(tmp_path, dsmc={"write_interval_s": 1.0e-4})
    assert run.write_interval_steps == round(1.0e-4 / run.delta_t_s)
    assert run.write_interval_s == pytest.approx(
        run.write_interval_steps * run.delta_t_s)
    # Within half a time step of what was asked for -- that is the whole of the
    # error a step count can introduce.
    assert abs(run.write_interval_s - 1.0e-4) <= 0.5 * run.delta_t_s
    assert run.n_steps % run.write_interval_steps == 0


def test_the_end_time_is_a_whole_number_of_steps(tmp_path):
    """Otherwise the last step overshoots and the write index misses it."""
    _, _, _, _, run = derived(tmp_path)
    assert run.end_time_s == pytest.approx(run.n_steps * run.delta_t_s)
    assert isinstance(run.n_steps, int)


def test_the_write_schedule_survives_a_restart(tmp_path):
    """`writeControl timeStep` counts the GLOBAL step index, which each time
    directory stores in uniform/time. So the write times are the same whether
    the run went straight through or was resumed part-way -- which `runTime`
    control cannot promise, because it measures from the resume point.

    The property that matters: endTime is a whole number of write intervals, so
    a write lands on it from any resume point.
    """
    _, _, _, _, run = derived(tmp_path)
    assert run.n_steps % run.write_interval_steps == 0
    assert run.n_writes == run.n_steps // run.write_interval_steps
    assert run.n_writes >= 2


def test_a_run_that_would_average_nothing_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="nothing would be averaged"):
        derived(tmp_path, dsmc={"average_start_s": 1.0, "end_time_s": 0.5})


def test_a_non_positive_cell_size_is_rejected(tmp_path):
    cfg = load_case(tmp_path)
    geom = from_config(cfg)
    state = inflow.from_config(cfg)
    with pytest.raises(ValueError, match="must both be positive"):
        inflow.derive_run_settings(cfg, state, geom, min_cell_size_m=0.0,
                                   exit_cell_volume_m3=1.0)


# --------------------------------------------------------------------------- #
# checks
# --------------------------------------------------------------------------- #

def _report(tmp_path, **sections):
    cfg, geom, state, plan, run = derived(tmp_path, **sections)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    particles = mesh.estimate_particles(
        plan, geom, state, run.n_equivalent_particles, samples_per_axis=8)
    return checks.run_checks(cfg, geom, state, plan, run,
                             nozzle=nozzle, particles=particles)


def _status(report, name):
    return next(r.status for r in report.results if r.name == name)


def test_a_sane_case_passes_every_check(tmp_path):
    report = _report(tmp_path)
    assert not report.failures, [r.message for r in report.failures]


def test_an_over_large_time_step_fails_hard(tmp_path):
    """A molecule crossing more than one smallest cell per step skips cells
    without being offered a collision partner in them."""
    report = _report(tmp_path, dsmc={"delta_t_s": 1.0e-3})
    assert _status(report, "time_step_courant") == "fail"
    with pytest.raises(checks.DsmcCheckError, match="time_step_courant"):
        report.raise_if_failed()


def test_a_marginal_time_step_only_warns(tmp_path):
    report = _report(tmp_path, dsmc={"courant_target": 0.8})
    assert _status(report, "time_step_courant") == "warn"
    report.raise_if_failed()      # a warning must not stop the run


def test_a_cell_coarser_than_the_mean_free_path_warns(tmp_path):
    report = _report(tmp_path, exit={"knudsen": 0.001},
                     mesh={"core_cell_size_m": 0.05})
    assert _status(report, "cell_over_mean_free_path") == "warn"


def test_too_few_nozzle_faces_fails_hard(tmp_path):
    """A disk resolved by a handful of squares is not a circular exit."""
    report = _report(tmp_path, checks={"min_nozzle_faces": 10000})
    assert _status(report, "nozzle_faces") == "fail"


def test_an_over_large_nozzle_area_error_fails_hard(tmp_path):
    """The injected molecule flow is proportional to the meshed area."""
    report = _report(tmp_path, checks={"max_nozzle_area_error": 1.0e-6})
    assert _status(report, "nozzle_area") == "fail"


def test_a_thin_exit_cell_warns(tmp_path):
    report = _report(tmp_path, dsmc={"n_equivalent_particles": 1.0e30})
    assert _status(report, "particles_per_cell") == "warn"


def test_a_huge_parcel_count_warns(tmp_path):
    report = _report(tmp_path, dsmc={"n_equivalent_particles": 1.0})
    assert _status(report, "particle_count") == "warn"


def test_disabling_the_checks_is_itself_a_warning(tmp_path):
    report = _report(tmp_path, checks={"enabled": False})
    assert _status(report, "checks_enabled") == "warn"
    assert len(report.results) == 1


def test_patch_checks_reject_a_wall_boundary(tmp_path):
    """A wall would reflect: the plume would expand into a closed box."""
    from plumetools.mesh.boundary import PatchInfo

    cfg, geom, state, plan, run = derived(tmp_path)
    patches = {
        "nozzle": PatchInfo("nozzle", "patch", 300, 0),
        "upstreamVacuum": PatchInfo("upstreamVacuum", "patch", 100, 300),
        "vacuum": PatchInfo("vacuum", "wall", 500, 400),
    }
    report = checks.run_checks(cfg, geom, state, plan, run, patches=patches)
    assert _status(report, "no_wall_boundaries") == "fail"


def test_patch_checks_reject_an_empty_nozzle(tmp_path):
    from plumetools.mesh.boundary import PatchInfo

    cfg, geom, state, plan, run = derived(tmp_path)
    patches = {
        "nozzle": PatchInfo("nozzle", "patch", 0, 0),
        "upstreamVacuum": PatchInfo("upstreamVacuum", "patch", 100, 0),
        "vacuum": PatchInfo("vacuum", "patch", 500, 100),
    }
    report = checks.run_checks(cfg, geom, state, plan, run, patches=patches)
    assert _status(report, "nozzle_not_empty") == "fail"


def test_patch_checks_reject_a_missing_nozzle(tmp_path):
    from plumetools.mesh.boundary import PatchInfo

    cfg, geom, state, plan, run = derived(tmp_path)
    patches = {"upstreamVacuum": PatchInfo("upstreamVacuum", "patch", 100, 0)}
    report = checks.run_checks(cfg, geom, state, plan, run, patches=patches)
    assert _status(report, "required_patches") == "fail"


def test_the_report_is_machine_readable(tmp_path):
    document = _report(tmp_path).as_dict()
    assert document["passed"] is True
    assert document["n_failures"] == 0
    assert all({"name", "status", "message"} <= set(r) for r in document["results"])


def test_every_failure_is_listed_not_just_the_first(tmp_path):
    report = _report(tmp_path, dsmc={"delta_t_s": 1.0e-3},
                     checks={"min_nozzle_faces": 10000})
    with pytest.raises(checks.DsmcCheckError) as excinfo:
        report.raise_if_failed()
    assert "time_step_courant" in str(excinfo.value)
    assert "nozzle_faces" in str(excinfo.value)


def test_write_summary_is_json(tmp_path):
    import json

    path = checks.write_summary(tmp_path, {"case": {"name": "Kn100"}})
    assert path.name == "case-summary.json"
    assert json.loads(path.read_text(encoding="utf-8"))["case"]["name"] == "Kn100"
