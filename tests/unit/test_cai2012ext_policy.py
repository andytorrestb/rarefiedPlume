"""The Knudsen extension's fixed policy, its feasibility verdicts and provenance.

Nothing here meshes or runs anything. What is covered is what decides *which
sixteen cases exist* and *whether any of them may run*: the exclusion of Cai's
own Knudsen numbers (even disabled, even in a hand-edited manifest), the exact
matrix, the settings inherited from cases/cai2012-health and the five that are
deliberately not, the fixed Kn = 0.01 startup reference, the particle level
applied once, and the production verdicts for the dense cases.
"""

from __future__ import annotations

import copy
import dataclasses
import math
from pathlib import Path

import pytest
import yaml

from plumetools.cai2012 import checks, inflow, mesh
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012.study import load_study
from plumetools.cai2012ext import (
    EXCLUDED_KNUDSEN,
    EXTENSION_KNUDSEN,
    PARTICLE_MULTIPLIERS,
    feasibility,
    policy,
    provenance,
)

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"
HEALTH = REPO / "cases" / "cai2012-health"
PARENT = REPO / "cases" / "cai2012"


def shipped_study_document() -> dict:
    return yaml.safe_load((FAMILY / "study.yaml").read_text(encoding="utf-8"))


def write_study(tmp_path, document) -> Path:
    path = tmp_path / "study.yaml"
    path.write_text(yaml.safe_dump(document), encoding="utf-8")
    return path


def template_at(knudsen: float, multiplier: float = 1.0):
    cfg = load_case_config(FAMILY / "baseCase")
    return dataclasses.replace(
        cfg, exit=dataclasses.replace(cfg.exit, knudsen=knudsen),
        dsmc=dataclasses.replace(cfg.dsmc,
                                 numerical_particle_multiplier=multiplier))


def derived(cfg):
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom, min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    particles = mesh.estimate_particles(plan, geom, exit_state,
                                        run.n_equivalent_particles)
    report = checks.run_checks(cfg, geom, exit_state, plan, run, nozzle=nozzle,
                               particles=particles)
    return cfg, geom, exit_state, plan, run, nozzle, particles, report


# --------------------------------------------------------------------------- #
# the constants
# --------------------------------------------------------------------------- #

def test_the_excluded_knudsen_numbers_are_exactly_cais():
    assert EXCLUDED_KNUDSEN == (0.01, 0.1, 100.0)


def test_the_extension_and_the_exclusions_do_not_overlap():
    assert not any(policy.is_excluded(k) for k in EXTENSION_KNUDSEN)
    assert EXTENSION_KNUDSEN == (1.0e-4, 1.0e-3, 1.0, 10.0)
    assert PARTICLE_MULTIPLIERS == (1.0,)          # the baseline only
    assert policy.expected_case_count() == 4


def test_a_knudsen_number_typed_as_a_decimal_is_still_recognised():
    assert policy.is_excluded(0.1 + 1e-15)
    assert policy.is_excluded(1.0e2)
    assert not policy.is_excluded(0.011)


# --------------------------------------------------------------------------- #
# the study
# --------------------------------------------------------------------------- #

def test_the_shipped_study_is_exactly_the_extension(tmp_path):
    study = load_study(FAMILY / "study.yaml")
    policy.check_study(study)
    cases = study.expand()
    assert len(cases) == 4
    assert sorted(c.knudsen for c in cases) == sorted(EXTENSION_KNUDSEN)
    assert {c.multiplier for c in cases} == {1.0}
    assert not any(policy.is_excluded(c.knudsen) for c in study.all_cases())


@pytest.mark.parametrize("excluded", [0.01, 0.1, 100.0])
def test_an_excluded_kn_is_refused_even_when_disabled(tmp_path, excluded):
    document = shipped_study_document()
    document["kn_cases"].append(
        {"name": f"Kn{excluded:g}".replace(".", "p"), "Kn": excluded,
         "enabled": False})
    study = load_study(write_study(tmp_path, document))
    with pytest.raises(policy.ExtensionError, match="Cai's own Knudsen"):
        policy.check_study(study)


def test_a_missing_knudsen_number_is_refused(tmp_path):
    document = shipped_study_document()
    document["kn_cases"] = [c for c in document["kn_cases"] if c["Kn"] != 1.0]
    with pytest.raises(policy.ExtensionError, match="does not enable Kn"):
        policy.check_study(load_study(write_study(tmp_path, document)))


def test_a_second_particle_level_is_refused(tmp_path):
    document = shipped_study_document()
    document["particle_levels"].append({"name": "np2x", "multiplier": 2.0})
    with pytest.raises(policy.ExtensionError, match="no parcel-population sweep"):
        policy.check_study(load_study(write_study(tmp_path, document)))


@pytest.mark.parametrize("key", ["run_time_multiplier",
                                 "output_frequency_multiplier"])
def test_the_parent_study_multipliers_are_not_inherited(tmp_path, key):
    document = shipped_study_document()
    document[key] = 3.0
    with pytest.raises(policy.ExtensionError, match=key):
        policy.check_study(load_study(write_study(tmp_path, document)))


@pytest.mark.parametrize("section,key,value", [
    ("gas", "omega", 0.81),
    ("exit", "T0_K", 350.0),
    ("dsmc", "sampling_domain_transits", 1.5),
    ("dsmc", "transient_basis", "transits"),
    ("resolution", "target_particles_per_cell", 20.0),
])
def test_overrides_that_would_break_an_invariant_are_refused(tmp_path, section,
                                                             key, value):
    document = shipped_study_document()
    document["kn_cases"][0]["overrides"] = {section: {key: value}}
    with pytest.raises(policy.ExtensionError, match="overrides"):
        policy.check_study(load_study(write_study(tmp_path, document)))


def test_a_per_kn_mesh_budget_is_an_allowed_override(tmp_path):
    document = shipped_study_document()
    document["kn_cases"][2]["overrides"] = {"mesh": {"max_cells": 40000000},
                                            "dsmc": {"n_subdomains": 256}}
    policy.check_study(load_study(write_study(tmp_path, document)))


# --------------------------------------------------------------------------- #
# the template and the inheritance
# --------------------------------------------------------------------------- #

def test_the_shipped_template_passes_the_policy():
    cfg = load_case_config(FAMILY / "baseCase")
    policy.check_template(cfg)
    assert not policy.is_excluded(cfg.exit.knudsen)


def test_a_template_with_the_health_startup_is_refused():
    cfg = load_case_config(FAMILY / "baseCase")
    bad = dataclasses.replace(cfg, dsmc=dataclasses.replace(
        cfg.dsmc, transient_basis="transits"))
    with pytest.raises(policy.ExtensionError, match="transient_basis"):
        policy.check_template(bad)


def test_a_template_with_another_startup_reference_is_refused():
    cfg = load_case_config(FAMILY / "baseCase")
    bad = dataclasses.replace(cfg, resolution=dataclasses.replace(
        cfg.resolution, reference_knudsen=1.0))
    with pytest.raises(policy.ExtensionError, match="reference_knudsen"):
        policy.check_template(bad)


def test_a_template_at_an_excluded_kn_is_refused():
    with pytest.raises(policy.ExtensionError, match="Cai's own"):
        policy.check_template(template_at(100.0))


def test_everything_but_five_declared_keys_is_inherited_from_health():
    report = policy.inheritance_report(load_case_config(FAMILY / "baseCase"),
                                       load_case_config(HEALTH / "baseCase"))
    assert report["undeclared"] == {}
    assert set(report["overridden"]) == set(policy.DELIBERATE_OVERRIDES)
    assert len(policy.DELIBERATE_OVERRIDES) == 5
    overridden = report["overridden"]
    assert overridden["resolution.target_particles_per_cell"]["extension"] == 40.0
    assert overridden["dsmc.sampling_domain_transits"]["extension"] == 4.5
    assert overridden["dsmc.transient_basis"] == {
        "extension": "cai", "health": "transits",
        "reason": policy.DELIBERATE_OVERRIDES["dsmc.transient_basis"]}
    assert overridden["dsmc.n_subdomains"] == {
        "extension": 192, "health": 12,
        "reason": policy.DELIBERATE_OVERRIDES["dsmc.n_subdomains"]}
    for key in ("dsmc.write_interval_s", "dsmc.courant_target",
                "mesh.max_cells", "mesh.target_cell_over_mfp",
                "mesh.max_core_cell_over_D", "checks.max_cell_over_mfp",
                "post.centerline_radius_over_D"):
        assert key in report["inherited"], key


def test_the_two_swept_values_are_health_rows_ppc040_and_s4p5():
    assert policy.check_health_study_rows(HEALTH / "study.yaml") == []


# --------------------------------------------------------------------------- #
# the manifest gate
# --------------------------------------------------------------------------- #

def _manifest(*rows):
    return {"cases": [dict(row) for row in rows]}


def test_a_manifest_naming_an_excluded_kn_schedules_nothing():
    manifest = _manifest(
        {"name": "Kn1_np1x", "path": "Cases/Kn1_np1x", "Kn": 1.0,
         "feasibility": {"status": "production-candidate"}},
        {"name": "Kn0p01_np1x", "path": "Cases/Kn0p01_np1x", "Kn": 0.01,
         "feasibility": {"status": "production-candidate"}})
    with pytest.raises(policy.ExtensionError, match="0.01"):
        policy.schedulable(manifest, allow_unsuitable=True)


def test_unsuitable_cases_need_an_explicit_flag():
    manifest = _manifest(
        {"name": "Kn1_np1x", "path": "Cases/Kn1_np1x", "Kn": 1.0,
         "feasibility": {"status": "production-candidate"}},
        {"name": "Kn0p001_np1x", "path": "Cases/Kn0p001_np1x", "Kn": 0.001,
         "feasibility": {"status": "unsuitable-for-production",
                         "unresolved": ["cell_over_mfp: ..."]}})
    runnable, refused = policy.schedulable(manifest)
    assert [e["name"] for e in runnable] == ["Kn1_np1x"]
    assert refused[0][1] == ["cell_over_mfp: ..."]
    runnable, refused = policy.schedulable(manifest, allow_unsuitable=True)
    assert len(runnable) == 2 and not refused


def test_the_command_line_gate_exits_2_on_an_excluded_kn(tmp_path, capsys):
    path = tmp_path / "manifest.yaml"
    path.write_text(yaml.safe_dump(_manifest(
        {"name": "Kn100_np1x", "path": "Cases/Kn100_np1x", "Kn": 100.0})),
        encoding="utf-8")
    assert policy.main([str(path), "--all"]) == 2
    assert capsys.readouterr().out == ""


# --------------------------------------------------------------------------- #
# feasibility
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("knudsen", [10.0, 1.0])
def test_the_rarefied_cases_are_production_candidates(knudsen):
    verdict = feasibility.assess(*derived(template_at(knudsen)))
    assert verdict["status"] == feasibility.STATUS_CANDIDATE
    assert verdict["unresolved"] == []
    assert verdict["resolution"]["n_cells"] == 1_310_720
    assert verdict["resolution"]["cell_over_mean_free_path"] == pytest.approx(
        0.01 / (knudsen * 0.2))


@pytest.mark.parametrize("knudsen", [1.0e-3, 1.0e-4])
def test_the_dense_cases_are_unsuitable_for_stated_reasons(knudsen):
    verdict = feasibility.assess(*derived(template_at(knudsen)))
    assert verdict["status"] == feasibility.STATUS_UNSUITABLE
    failed = {c["name"] for c in verdict["criteria"] if not c["passed"]}
    assert failed == {"cell_over_mean_free_path", "mesh_budget",
                      "time_step_collision_time"}
    assert verdict["resolution"]["coarsened_for_budget"]
    # The physical input is exact: n0 is Kn's, however coarse the mesh.
    assert verdict["resolution"]["cells_to_meet_criterion"] > 1.0e10
    assert verdict["time_step"]["dt_over_exit_collision_time"] > 1.0


def test_the_criteria_are_the_nozzle_plume_checks_with_its_thresholds():
    from plumetools.nozzleplume.config import ChecksConfig

    verdict = feasibility.assess(*derived(template_at(1.0)))
    criteria = {c["name"]: c for c in verdict["criteria"]}
    assert set(criteria) == {
        "dsmc_checks", "mesh_budget", "cell_over_mean_free_path",
        "time_step_collision_time", "time_step_courant", "transient_length",
        "internal_energy_exchange", "particles_per_cell", "parcel_count"}
    nozzle = ChecksConfig()
    assert criteria["time_step_collision_time"]["limit"] == \
        nozzle.max_time_step_over_collision_time
    assert criteria["parcel_count"]["limit"] == nozzle.warn_parcels
    # Cai's startup is 17.7 domain transits; the jet needs one to settle.
    assert criteria["transient_length"]["passed"]
    assert verdict["time_step"]["average_start_over_settling_time"] == \
        pytest.approx(17.7245, rel=1e-4)


def test_more_parcels_do_not_change_any_numerical_criterion():
    one = feasibility.assess(*derived(template_at(1.0e-3, 1.0)))
    ten = feasibility.assess(*derived(template_at(1.0e-3, 10.0)))
    numerical = lambda v: {c["name"]: (c["value"], c["passed"])  # noqa: E731
                           for c in v["criteria"] if c["kind"] == "numerical"}
    assert numerical(one) == numerical(ten)
    assert ten["particles"]["estimated_parcels"] == pytest.approx(
        10.0 * one["particles"]["estimated_parcels"])


def test_the_cost_estimate_uses_the_health_measurements():
    model, origin = feasibility.load_cost_model(FAMILY, HEALTH)
    assert model.n_measurements >= 1
    verdict = feasibility.assess(*derived(template_at(10.0)), cost_model=model,
                                 cost_origin=origin)
    cost = verdict["cost"]
    assert cost["predicted_wall_hours"] == pytest.approx(model.seconds(
        parcels=cost["estimated_parcels"], cells=cost["n_cells"],
        steps=cost["n_steps"]) / 3600.0)
    assert "lower bound" in cost["cost_caveat"] or "longer" in cost["cost_caveat"]


# --------------------------------------------------------------------------- #
# provenance: the physical state, recomputed
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("knudsen", EXTENSION_KNUDSEN)
def test_n0_reproduces_the_target_kn_independently(knudsen):
    cfg = template_at(knudsen)
    exit_state = inflow.from_config(cfg)
    assert provenance.independent_knudsen(
        exit_state.number_density_per_m3, cfg) == pytest.approx(knudsen,
                                                                rel=1e-12)
    record = provenance.exit_conditions(cfg, from_config(cfg), exit_state)
    assert record["p0_Pa"] == pytest.approx(
        exit_state.number_density_per_m3 * 1.380649e-23 * 300.0)
    assert record["rho0_kg_per_m3"] == pytest.approx(
        exit_state.number_density_per_m3 * 6.63e-26)
    assert record["S0_recomputed"] == pytest.approx(2.0, rel=1e-12)


def test_the_density_scales_inversely_with_kn():
    densities = [inflow.from_config(template_at(k)).number_density_per_m3
                 for k in EXTENSION_KNUDSEN]
    for k, n in zip(EXTENSION_KNUDSEN, densities):
        assert n * k == pytest.approx(densities[0] * EXTENSION_KNUDSEN[0],
                                      rel=1e-12)


def test_the_startup_is_ten_thousand_kn_0p01_collision_times_for_every_kn():
    lambda_ref = 0.01 * 0.2
    mean_speed = math.sqrt(8.0 * 1.380649e-23 * 300.0 / (math.pi * 6.63e-26))
    expected = 1.0e4 * lambda_ref / mean_speed
    for knudsen in EXTENSION_KNUDSEN:
        _, _, exit_state, _, run, *_ = derived(template_at(knudsen))
        assert run.average_start_s == pytest.approx(expected, rel=1e-12)
        record = provenance.startup(template_at(knudsen), run, exit_state)
        assert record["achieved_reference_collision_times"] == pytest.approx(1e4)
        # ... which is a very different number of the case's OWN collision times.
        assert record["exit_collision_times"] == pytest.approx(
            1.0e4 * 0.01 / knudsen, rel=1e-9)


def test_averaging_is_four_and_a_half_transits_after_the_startup():
    _, _, _, _, run, *_ = derived(template_at(1.0))
    sampled = run.end_time_s - run.average_start_s
    assert 4.5 * run.domain_transit_s <= sampled + 1e-15
    assert sampled <= 4.5 * run.domain_transit_s + run.write_interval_s
    assert run.write_interval_s == pytest.approx(3.536e-4, rel=0.01)


def test_the_particle_level_is_applied_exactly_once():
    for multiplier in PARTICLE_MULTIPLIERS:
        cfg = template_at(1.0, multiplier)
        _, _, exit_state, plan, run, *_ = derived(cfg)
        assert run.n_equivalent_particles * multiplier * 40.0 == pytest.approx(
            exit_state.number_density_per_m3 * plan.exit_cell_volume_m3,
            rel=1e-12)
        assert run.exit_particles_per_cell == pytest.approx(40.0 * multiplier)


# --------------------------------------------------------------------------- #
# reuse, not forks
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", ["Allclean", "Allmesh", "Allpost", "Allrun",
                                  "postProcess.py", "runInflow.py"])
def test_the_case_scripts_are_cai2012s_byte_for_byte(name):
    assert (FAMILY / "baseCase" / name).read_bytes() == (
        PARENT / "baseCase" / name).read_bytes()


def test_the_audit_script_is_the_health_familys_byte_for_byte():
    assert (FAMILY / "baseCase" / "Allaudit.py").read_bytes() == (
        HEALTH / "baseCase" / "Allaudit.py").read_bytes()


def test_allcleancases_is_cai2012s_byte_for_byte():
    assert (FAMILY / "AllcleanCases").read_bytes() == (
        PARENT / "AllcleanCases").read_bytes()
