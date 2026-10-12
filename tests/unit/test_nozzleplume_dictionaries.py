"""``dsmcProperties``, ``controlDict`` and the other generated dictionaries.

What these protect is the mapping from ``case.yaml`` to what the solver reads:
the particle weight written exactly as given, the species as given, the inflow
restricted to the nozzle patch, and the model names in the spelling ``dsmcFoam``
selects them by -- one of which is not the spelling anybody would guess.
"""

from __future__ import annotations

import re

import pytest

from _nozzleplume import R4D_CASE, derive, derived
from plumetools.nozzleplume import dictionaries
from plumetools.nozzleplume.config import load_case_config


def entry(text: str, key: str) -> str:
    """The value of the first top-level ``key value;`` entry in a dictionary."""
    match = re.search(rf"^{re.escape(key)}\s+([^;]+);", text, re.M)
    assert match, f"no '{key}' entry"
    return match.group(1).strip()


def properties(tmp_path, **sections) -> str:
    cfg, geom, exit_state, _, _, run = derived(tmp_path, **sections)
    return dictionaries.render_dsmc_properties(cfg, geom, exit_state, run)


# --------------------------------------------------------------------------- #
# dsmcProperties
# --------------------------------------------------------------------------- #

def test_the_particle_weight_is_written_exactly_as_given(tmp_path):
    text = properties(tmp_path, dsmc={"n_equivalent_particles": 3.7e16})
    assert float(entry(text, "nEquivalentParticles")) == 3.7e16
    assert "Set directly: dsmc.n_equivalent_particles" in text


def test_a_derived_weight_says_what_it_was_derived_from(tmp_path):
    text = properties(tmp_path, dsmc={"n_equivalent_particles": None})
    assert "Derived from dsmc.target_particles_per_exit_cell = 20" in text


def test_the_shipped_r4d_weight_reaches_dsmc_properties():
    cfg = load_case_config(R4D_CASE)
    geom, exit_state, _, _, run = derive(cfg)
    text = dictionaries.render_dsmc_properties(cfg, geom, exit_state, run)
    assert float(entry(text, "nEquivalentParticles")) == (
        cfg.dsmc.n_equivalent_particles)


def test_the_species_is_written_as_given(tmp_path):
    text = properties(tmp_path)
    assert entry(text, "typeIdList") == "(exhaust)"
    block = text.split("moleculeProperties")[1]
    assert re.search(r"mass\s+3\.4e-26;", block)
    assert re.search(r"diameter\s+3\.5e-10;", block)
    assert re.search(r"internalDegreesOfFreedom\s+0;", block)
    assert re.search(r"omega\s+0\.5;", block)


def test_the_r4d_molecular_mass_keeps_all_its_figures():
    """3.40055e-26 has six; a %g would write 3.40055e-26 only by luck."""
    cfg = load_case_config(R4D_CASE)
    geom, exit_state, _, _, run = derive(cfg)
    text = dictionaries.render_dsmc_properties(cfg, geom, exit_state, run)
    assert re.search(r"mass\s+3\.40055e-26;", text)


def test_inflow_is_restricted_to_the_nozzle_patch(tmp_path):
    """The whole reason for plumeFieldInflow: stock FreeStream would inject on
    every open boundary."""
    text = properties(tmp_path)
    assert entry(text, "InflowBoundaryModel") == "plumeFieldInflow"
    coeffs = text.split("plumeFieldInflowCoeffs")[1].split("typeIdList")[0]
    assert "patches ( nozzle );" in coeffs
    assert re.search(r"exhaust\s+boundaryNumberDensity_exhaust;", coeffs)


def test_a_renamed_nozzle_patch_is_the_one_injected_through(tmp_path):
    text = properties(tmp_path, mesh={"patch_names": {"nozzle": "exitPlane"}})
    assert "patches ( exitPlane );" in text


def test_elastic_vhs_is_the_default_collision_model(tmp_path):
    text = properties(tmp_path)
    assert entry(text, "BinaryCollisionModel") == "VariableHardSphere"
    assert re.search(r"VariableHardSphereCoeffs\s*\{\s*Tref\s+273;\s*\}", text)


def test_larsen_borgnakke_carries_its_relaxation_number(tmp_path):
    text = properties(tmp_path, gas={"internal_degrees_of_freedom": 2}, dsmc={
        "collision_model": "LarsenBorgnakkeVariableHardSphere",
        "relaxation_collision_number": 7.0})
    assert entry(text, "BinaryCollisionModel") == "LarsenBorgnakkeVariableHardSphere"
    coeffs = text.split("LarsenBorgnakkeVariableHardSphereCoeffs")[1].split("}")[0]
    assert re.search(r"relaxationCollisionNumber\s+7;", coeffs)
    assert re.search(r"internalDegreesOfFreedom\s+2;", text)


def test_the_collisionless_model_is_written_under_the_name_dsmcfoam_knows(tmp_path):
    """The class is NoBinaryCollision; its type name is 'none'. v2512 rejects
    the class name with "Unknown BinaryCollisionModel type"."""
    text = properties(tmp_path, dsmc={"collision_model": "NoBinaryCollision"})
    assert entry(text, "BinaryCollisionModel") == "none"
    assert not re.search(r"^BinaryCollisionModel\s+NoBinaryCollision", text, re.M)


def test_the_diffuse_wall_is_the_default_and_names_its_temperature(tmp_path):
    text = properties(tmp_path, plate={"temperature_K": 350.0})
    assert entry(text, "WallInteractionModel") == "MaxwellianThermal"
    assert "350 K" in text
    assert "MaxwellianThermalCoeffs" not in text      # it has none


def test_the_mixed_wall_carries_its_diffuse_fraction(tmp_path):
    text = properties(tmp_path, dsmc={
        "wall_interaction_model": "MixedDiffuseSpecular", "diffuse_fraction": 0.8})
    assert entry(text, "WallInteractionModel") == "MixedDiffuseSpecular"
    assert re.search(
        r"MixedDiffuseSpecularCoeffs\s*\{\s*diffuseFraction\s+0\.8;\s*\}", text)


def test_a_free_plume_still_names_a_wall_model(tmp_path):
    """DSMCCloud constructs one unconditionally."""
    text = properties(tmp_path, plate=None)
    assert entry(text, "WallInteractionModel") == "MaxwellianThermal"
    assert "NO wall patches" in text


# --------------------------------------------------------------------------- #
# controlDict
# --------------------------------------------------------------------------- #

def control(tmp_path, **sections) -> tuple:
    cfg, *_, run = derived(tmp_path, **sections)
    return dictionaries.render_control_dict(cfg, run), run


def test_control_dict_carries_the_step_and_the_end_the_run_counted(tmp_path):
    text, run = control(tmp_path)
    assert float(entry(text, "deltaT")) == run.delta_t_s
    assert float(entry(text, "endTime")) == pytest.approx(run.end_time_s, rel=1e-9)
    assert float(entry(text, "endTime")) / float(entry(text, "deltaT")) == (
        pytest.approx(run.n_steps, abs=1e-6))


def test_control_dict_writes_on_the_step_count(tmp_path):
    """Not runTime: that schedule restarts at a resume and can write nothing."""
    text, run = control(tmp_path)
    assert entry(text, "writeControl") == "timeStep"
    assert int(entry(text, "writeInterval")) == run.write_interval_steps
    assert entry(text, "purgeWrite") == "0"


def test_control_dict_resumes_and_loads_the_inflow_library(tmp_path):
    text, _ = control(tmp_path)
    assert entry(text, "startFrom") == "latestTime"
    assert entry(text, "libs") == '( "libplumeDsmcBoundaryModels.so" )'
    assert entry(text, "writeFormat") == "ascii"       # the post reads text


def test_control_dict_averages_the_plate_loads_from_the_sampling_start(tmp_path):
    text, run = control(tmp_path)
    functions = text.split("functions")[1]
    assert float(re.search(r"timeStart\s+([^;]+);", functions).group(1)) == (
        pytest.approx(run.average_start_s))
    for field in ("fD", "q", "rhoN", "rhoM", "momentum", "linearKE", "dsmcRhoN"):
        assert re.search(rf"^\s+{field}\s*$", functions, re.M), field


def test_only_the_number_densities_keep_their_variance(tmp_path):
    """It is the measured statistical error; a symmTensor per cell for every
    vector field would be written for nobody."""
    text, _ = control(tmp_path)
    blocks = dict(re.findall(r"^\s+(\w+)\s*\{[^}]*prime2Mean\s+(\w+);", text, re.M))
    assert {k for k, v in blocks.items() if v == "on"} == {"rhoN", "dsmcRhoN"}


# --------------------------------------------------------------------------- #
# the rest
# --------------------------------------------------------------------------- #

def test_the_initial_fill_is_the_configured_near_vacuum(tmp_path):
    cfg, *_ = derived(tmp_path, dsmc={"initial_number_density_per_m3": 5.0e9})
    text = dictionaries.render_dsmc_initialise_dict(cfg)
    assert re.search(r"exhaust\s+5e\+09;", text)
    assert entry(text, "velocity") == "(0 0 0)"


def test_scotch_is_the_default_decomposition(tmp_path):
    cfg, *_ = derived(tmp_path, dsmc={"n_subdomains": 12})
    text = dictionaries.render_decompose_par_dict(cfg)
    assert (entry(text, "numberOfSubdomains"), entry(text, "method")) == (
        "12", "scotch")
    assert "coeffs" not in text


def test_simple_decomposition_carries_a_split_that_multiplies_out(tmp_path):
    cfg, *_ = derived(tmp_path, dsmc={"n_subdomains": 12,
                                      "decomposition_method": "simple"})
    text = dictionaries.render_decompose_par_dict(cfg)
    nx, ny, nz = (int(v) for v in re.search(
        r"n\s+\((\d+) (\d+) (\d+)\);", text).groups())
    assert nx * ny * nz == 12
    assert nx <= ny and nx <= nz        # lateral cuts first


@pytest.mark.parametrize("n,split", [(1, (1, 1, 1)), (4, (1, 2, 2)),
                                     (6, (1, 3, 2)), (12, (2, 3, 2)),
                                     (7, (1, 7, 1))])
def test_simple_split_prefers_lateral_cuts(n, split):
    assert dictionaries.simple_split(n) == split


def test_simple_split_rejects_no_ranks():
    with pytest.raises(ValueError, match="n_subdomains"):
        dictionaries.simple_split(0)


def test_every_dictionary_is_written_with_unix_line_endings(tmp_path):
    cfg, geom, exit_state, _, _, run = derived(tmp_path)
    written = dictionaries.write_case_dictionaries(tmp_path, cfg, geom, exit_state,
                                                   run)
    assert [p.name for p in written] == [
        "dsmcProperties", "controlDict", "dsmcInitialiseDict", "decomposeParDict",
        "fvSchemes", "fvSolution"]
    for path in written:
        raw = path.read_bytes()
        assert b"\r" not in raw
        assert f"object      {path.name};".encode() in raw
