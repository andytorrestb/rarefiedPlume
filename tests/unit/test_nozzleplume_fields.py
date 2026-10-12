"""The ``0/`` fields, and the DSMC quality checks.

The field writer's load-bearing details are all about the plate: its temperature
has to be in ``boundaryT`` for the wall model to read, and ``fD`` and ``q`` have
to be ``calculated`` there or the loads the solver measures are thrown away at
write time. Both failures run to completion and look like results.
"""

from __future__ import annotations

import re

import pytest

from _nozzleplume import R4D_CASE, derive, derived
from plumetools.mesh.boundary import PatchInfo
from plumetools.nozzleplume import checks, foamfields
from plumetools.nozzleplume.config import load_case_config


def mesh_patches(nozzle_faces: int = 32, *, plate: str | None = "wall",
                 symmetry: bool = False, **types) -> dict:
    """The patches a built mesh has, in the order the pipeline creates them."""
    kinds = {"upstreamVacuum": "patch", "vacuum": "patch"}
    if symmetry:
        kinds["symmetry"] = "symmetry"
    if plate is not None:
        kinds["plate"] = plate
    kinds["nozzle"] = "patch"
    kinds.update(types)
    counts = {"nozzle": nozzle_faces}
    return {name: PatchInfo(name=name, type=kind, n_faces=counts.get(name, 100),
                            start_face=1000 * i)
            for i, (name, kind) in enumerate(kinds.items())}


def boundary_entry(text: str, patch: str) -> str:
    """One patch's entry in a field's boundaryField."""
    match = re.search(rf"^\s+{patch}\s*\{{([^}}]*)\}}", text.split("boundaryField")[1],
                      re.M)
    assert match, f"no boundaryField entry for {patch}"
    return " ".join(match.group(1).split())


def write_fields(tmp_path, patches=None, **sections) -> dict:
    cfg, geom, exit_state, *_ = derived(tmp_path, **sections)
    patches = patches if patches is not None else mesh_patches(
        plate=None if geom.plate is None else "wall",
        symmetry=bool(geom.symmetry_planes))
    written = foamfields.write_all(tmp_path, cfg, geom, exit_state, patches)
    return {path.name: path.read_text(encoding="utf-8") for path in written}


# --------------------------------------------------------------------------- #
# the inflow fields
# --------------------------------------------------------------------------- #

def test_all_twelve_fields_are_written(tmp_path):
    fields = write_fields(tmp_path)
    assert set(fields) == {
        "boundaryU", "boundaryT", "boundaryNumberDensity_exhaust", "dsmcRhoN",
        "fD", "iDof", "internalE", "linearKE", "momentum", "q", "rhoM", "rhoN"}


def test_the_nozzle_carries_the_exit_state(tmp_path):
    fields = write_fields(tmp_path)
    assert boundary_entry(fields["boundaryU"], "nozzle") == (
        "type fixedValue; value uniform (3000 0 0);")
    assert boundary_entry(fields["boundaryT"], "nozzle") == (
        "type fixedValue; value uniform 500;")
    assert boundary_entry(fields["boundaryNumberDensity_exhaust"], "nozzle") == (
        "type fixedValue; value uniform 1e+22;")


def test_the_r4d_nozzle_carries_the_given_values_to_full_precision(tmp_path):
    cfg = load_case_config(R4D_CASE)
    geom, exit_state, *_ = derive(cfg)
    written = foamfields.write_inflow_fields(tmp_path, cfg, geom, exit_state,
                                             mesh_patches(112))
    fields = {p.name: p.read_text(encoding="utf-8") for p in written}
    assert "uniform (3187.19 0 0);" in fields["boundaryU"]
    assert "uniform 498.04;" in fields["boundaryT"]
    assert "uniform 1.306e+22;" in fields["boundaryNumberDensity_exhaust"]


def test_boundary_t_is_a_scalar_field(tmp_path):
    """Standard dsmcFoam; the MNF fork reads a vector there."""
    assert "class       volScalarField;" in write_fields(tmp_path)["boundaryT"]


def test_the_plate_carries_its_temperature_and_no_velocity(tmp_path):
    """MaxwellianThermal reads both from HERE. zeroGradient would give 0 K, and
    a molecule re-emitted at zero thermal speed sticks to the wall."""
    fields = write_fields(tmp_path, plate={"temperature_K": 350.0})
    assert boundary_entry(fields["boundaryT"], "plate") == (
        "type fixedValue; value uniform 350;")
    assert boundary_entry(fields["boundaryU"], "plate") == (
        "type fixedValue; value uniform (0 0 0);")


def test_open_boundaries_carry_nothing(tmp_path):
    """They are not inflows: plumeFieldInflow has an explicit patch list."""
    fields = write_fields(tmp_path)
    for name in ("boundaryU", "boundaryT", "boundaryNumberDensity_exhaust"):
        for patch in ("vacuum", "upstreamVacuum"):
            assert boundary_entry(fields[name], patch) == "type zeroGradient;"


def test_only_the_nozzle_has_a_number_density(tmp_path):
    fields = write_fields(tmp_path)
    assert boundary_entry(fields["boundaryNumberDensity_exhaust"], "plate") == (
        "type zeroGradient;")


def test_a_symmetry_patch_gets_the_constraint_type(tmp_path):
    """OpenFOAM rejects anything else: 'inconsistent patch and patchField'."""
    fields = write_fields(tmp_path, domain={"symmetry_planes": ["y"]})
    for text in fields.values():
        assert boundary_entry(text, "symmetry") == "type symmetry;"


# --------------------------------------------------------------------------- #
# the measurement fields
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("field,zero", [("fD", "(0 0 0)"), ("q", "0")])
def test_the_plate_loads_are_calculated_on_the_plate(tmp_path, field, zero):
    """hitWallPatch accumulates into these. zeroGradient writes no values, so
    the load -- the result -- would be discarded at write time."""
    fields = write_fields(tmp_path)
    assert boundary_entry(fields[field], "plate") == (
        f"type calculated; value uniform {zero};")
    assert boundary_entry(fields[field], "vacuum") == "type zeroGradient;"


def test_the_other_measurement_fields_are_neutral_everywhere(tmp_path):
    fields = write_fields(tmp_path)
    for name in ("rhoN", "rhoM", "dsmcRhoN", "momentum", "linearKE", "internalE",
                 "iDof"):
        for patch in ("plate", "vacuum", "nozzle", "upstreamVacuum"):
            assert boundary_entry(fields[name], patch) == "type zeroGradient;"


def test_fields_have_the_classes_and_dimensions_dsmcfoam_declares(tmp_path):
    fields = write_fields(tmp_path)
    assert "class       volVectorField;" in fields["fD"]
    assert "dimensions      [1 -1 -2 0 0 0 0];" in fields["fD"]     # Pa
    assert "dimensions      [1 0 -3 0 0 0 0];" in fields["q"]       # W/m^2
    assert "dimensions      [0 -3 0 0 0 0 0];" in fields["rhoN"]


# --------------------------------------------------------------------------- #
# a mesh that is not the one case.yaml describes
# --------------------------------------------------------------------------- #

def test_a_mesh_without_the_nozzle_patch_is_refused(tmp_path):
    patches = mesh_patches()
    del patches["nozzle"]
    with pytest.raises(ValueError, match="nozzle patch"):
        write_fields(tmp_path, patches)


def test_an_empty_nozzle_patch_is_refused(tmp_path):
    with pytest.raises(ValueError, match="zero faces"):
        write_fields(tmp_path, mesh_patches(0))


def test_a_mesh_without_the_plate_case_yaml_has_is_refused(tmp_path):
    """It would run, with no plate in it."""
    with pytest.raises(ValueError, match="has a plate"):
        write_fields(tmp_path, mesh_patches(plate=None))


def test_a_wall_case_yaml_does_not_describe_is_refused(tmp_path):
    """There is no temperature to give it, and 0 K is not a temperature."""
    with pytest.raises(ValueError, match="wall patch"):
        write_fields(tmp_path, mesh_patches(), plate=None)


# --------------------------------------------------------------------------- #
# the quality checks
# --------------------------------------------------------------------------- #

def run_checks(tmp_path, patches=None, **sections) -> dict:
    cfg, geom, exit_state, plan, nozzle, run = derived(tmp_path, **sections)
    report = checks.run_checks(cfg, geom, exit_state, plan, run, nozzle=nozzle,
                               patches=patches)
    return {r.name: r for r in report.results}


def test_the_minimal_case_fails_nothing(tmp_path):
    results = run_checks(tmp_path)
    assert not [r for r in results.values() if r.status == "fail"]


def test_a_dense_exit_warns_about_the_cell_and_the_step_but_does_not_fail(tmp_path):
    """At Kn ~ 1e-4 neither criterion can be met; hiding that would be worse
    than reporting it, and refusing to run would be useless."""
    results = run_checks(tmp_path)
    assert results["cell_over_mean_free_path"].status == "warn"
    assert results["time_step_collision_time"].status == "warn"
    assert results["cell_over_mean_free_path"].value > 100.0


def test_the_cell_warning_says_how_many_levels_short_the_exit_is(tmp_path):
    cfg, geom, exit_state, plan, nozzle, run = derived(tmp_path)
    report = checks.run_checks(cfg, geom, exit_state, plan, run, nozzle=nozzle)
    message = {r.name: r for r in report.results}[
        "cell_over_mean_free_path"].message
    levels = int(re.search(r"take (\d+) more", message).group(1))
    assert plan.exit_cell_size_m / 2 ** levels <= exit_state.mean_free_path_m
    assert plan.exit_cell_size_m / 2 ** (levels - 1) > exit_state.mean_free_path_m


def test_a_rarefied_exit_meets_both_criteria(tmp_path):
    """The same checks pass where they can: Kn = 1 on the exit diameter."""
    results = run_checks(
        tmp_path, exit={"number_density_per_m3": None, "knudsen": 1.0})
    assert results["cell_over_mean_free_path"].status == "pass"
    assert results["time_step_collision_time"].status == "pass"


def test_a_step_that_skips_cells_is_a_hard_failure(tmp_path):
    results = run_checks(tmp_path, dsmc={"delta_t_s": 1.0e-4})
    assert results["time_step_courant"].status == "fail"
    assert "dsmc.delta_t_s below" in results["time_step_courant"].message


def test_a_starved_exit_cell_warns_with_the_weight_that_would_fix_it(tmp_path):
    cfg, geom, exit_state, plan, nozzle, run = derived(
        tmp_path, dsmc={"n_equivalent_particles": 1.0e18})
    report = checks.run_checks(cfg, geom, exit_state, plan, run, nozzle=nozzle)
    result = {r.name: r for r in report.results}["particles_per_cell"]
    assert result.status == "warn"
    suggested = float(re.search(r"nEquivalentParticles = ([\d.e+]+)",
                                result.message).group(1))
    assert exit_state.number_density_per_m3 * plan.exit_cell_volume_m3 / suggested == (
        pytest.approx(cfg.checks.min_particles_per_cell, rel=1e-3))


def test_too_many_parcels_warns(tmp_path):
    results = run_checks(tmp_path, dsmc={"n_equivalent_particles": 1.0e13})
    assert results["parcel_count"].status == "warn"


def test_averaging_before_the_flow_settles_warns(tmp_path):
    results = run_checks(tmp_path, dsmc={"average_start_s": 1.0e-4,
                                         "end_time_s": 1.0e-2})
    assert results["transient_length"].status == "warn"


def test_internal_modes_under_elastic_collisions_warn(tmp_path):
    """They would keep their injected energy for ever."""
    results = run_checks(tmp_path, gas={"internal_degrees_of_freedom": 2})
    assert results["internal_energy_exchange"].status == "warn"
    relaxing = run_checks(tmp_path / "b", gas={"internal_degrees_of_freedom": 2},
                          dsmc={"collision_model":
                                "LarsenBorgnakkeVariableHardSphere"})
    assert relaxing["internal_energy_exchange"].status == "pass"


def test_a_nozzle_of_a_handful_of_squares_is_a_hard_failure(tmp_path):
    results = run_checks(tmp_path, checks={"min_nozzle_faces": 1000})
    assert results["nozzle_faces"].status == "fail"


def test_a_nozzle_area_outside_tolerance_is_a_hard_failure(tmp_path):
    results = run_checks(tmp_path, checks={"max_nozzle_area_error": 0.001})
    assert results["nozzle_area"].status == "fail"


def test_disabling_the_checks_is_itself_reported(tmp_path):
    results = run_checks(tmp_path, checks={"enabled": False})
    assert list(results) == ["checks_enabled"]
    assert results["checks_enabled"].status == "warn"


# --------------------------------------------------------------------------- #
# the checks against a built mesh
# --------------------------------------------------------------------------- #

def test_the_planned_mesh_passes_the_patch_checks(tmp_path):
    results = run_checks(tmp_path, mesh_patches(32))
    for name in ("required_patches", "patch_types", "nozzle_matches_plan"):
        assert results[name].status == "pass", name


def test_a_plate_that_is_not_a_wall_is_a_hard_failure(tmp_path):
    """Declared 'patch' it would absorb the plume and feel no load."""
    results = run_checks(tmp_path, mesh_patches(32, plate="patch"))
    assert results["patch_types"].status == "fail"
    assert "feel no load" in results["patch_types"].message


def test_an_open_boundary_that_is_a_wall_is_a_hard_failure(tmp_path):
    """The plume would expand into a box."""
    results = run_checks(tmp_path, mesh_patches(32, vacuum="wall"))
    assert results["patch_types"].status == "fail"
    assert results["no_stray_walls"].status == "fail"


def test_a_nozzle_that_is_a_wall_is_a_hard_failure(tmp_path):
    results = run_checks(tmp_path, mesh_patches(32, nozzle="wall"))
    assert results["patch_types"].status == "fail"


def test_a_missing_plate_patch_is_a_hard_failure(tmp_path):
    results = run_checks(tmp_path, mesh_patches(32, plate=None))
    assert results["required_patches"].status == "fail"


def test_a_nozzle_with_the_wrong_face_count_is_not_the_planned_mesh(tmp_path):
    results = run_checks(tmp_path, mesh_patches(30))
    assert results["nozzle_matches_plan"].status == "fail"


def test_an_empty_patch_is_a_hard_failure(tmp_path):
    results = run_checks(tmp_path, mesh_patches(0))
    assert results["patches_not_empty"].status == "fail"


def test_raise_if_failed_lists_every_failure(tmp_path):
    cfg, geom, exit_state, plan, nozzle, run = derived(tmp_path)
    report = checks.run_checks(
        cfg, geom, exit_state, plan, run, nozzle=nozzle,
        patches=mesh_patches(30, vacuum="wall"))
    with pytest.raises(checks.DsmcCheckError) as excinfo:
        report.raise_if_failed()
    for name in ("patch_types", "no_stray_walls", "nozzle_matches_plan"):
        assert name in str(excinfo.value)


def test_the_summary_is_written_as_json(tmp_path):
    import json

    cfg, geom, exit_state, plan, nozzle, run = derived(tmp_path)
    report = checks.run_checks(cfg, geom, exit_state, plan, run, nozzle=nozzle)
    path = checks.write_summary(tmp_path, {
        "exit": exit_state.as_dict(), "mesh": plan.as_dict(), "run": run.as_dict(),
        "nozzle": nozzle.as_dict(), "checks": report.as_dict()})
    summary = json.loads(path.read_text(encoding="utf-8"))
    assert summary["run"]["n_equivalent_particles"] == 1.0e17
    assert summary["mesh"]["exit_level"] == 3
    assert summary["checks"]["passed"] is True
