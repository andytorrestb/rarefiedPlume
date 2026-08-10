"""The generated ``0/`` fields and OpenFOAM dictionaries.

The assertions here are mostly about **what the solver reads**, not about
formatting: a ``zeroGradient`` on the nozzle aborts the run, a ``wall`` on a
vacuum boundary silently changes the physics, and a ``FreeStream`` block would
produce a case that runs and is wrong.
"""

from __future__ import annotations

import pytest

from _cai2012 import derived
from plumetools.cai2012 import dictionaries, foamfields
from plumetools.mesh.boundary import PatchInfo

PATCHES = {
    "nozzle": PatchInfo("nozzle", "patch", 316, 0),
    "upstreamVacuum": PatchInfo("upstreamVacuum", "patch", 15684, 316),
    "vacuum": PatchInfo("vacuum", "patch", 57344, 16000),
}

HALF_PATCHES = dict(PATCHES, symmetry=PatchInfo("symmetry", "symmetry", 400, 73344))


def _uncommented(text: str) -> str:
    return "\n".join(line.split("//", 1)[0] for line in text.splitlines())


# --------------------------------------------------------------------------- #
# 0/ inflow fields
# --------------------------------------------------------------------------- #

@pytest.fixture
def written(tmp_path):
    cfg, geom, state, plan, run = derived(tmp_path)
    paths = foamfields.write_all(tmp_path, cfg, state, PATCHES)
    return cfg, state, run, {p.name: p.read_text(encoding="utf-8") for p in paths}


def test_the_three_inflow_fields_and_nine_measurement_fields_are_written(written):
    _, _, _, files = written
    assert "boundaryU" in files
    assert "boundaryT" in files
    assert "boundaryNumberDensity_Ar" in files
    assert len(files) == 3 + len(foamfields.MEASUREMENT_FIELDS)


def test_boundary_T_is_a_volScalarField(written):
    """Standard dsmcFoam declares it as one; the MNF fork takes a vector.
    Writing the wrong class is a hard read failure, which is the good outcome."""
    _, _, _, files = written
    assert "class       volScalarField;" in files["boundaryT"]


def test_the_nozzle_carries_the_derived_exit_state(written):
    _, state, _, files = written
    assert f"uniform {state.T0_K:.10g};" in files["boundaryT"]
    assert f"uniform ({state.velocity_m_per_s:.10g} 0 0);" in files["boundaryU"]
    assert f"uniform {state.number_density_per_m3:.10g};" in (
        files["boundaryNumberDensity_Ar"])


def test_the_nozzle_entries_are_fixedValue(written):
    """zeroGradient evaluates to the zero internal field, and plumeFieldInflow
    aborts with 'Zero boundary temperature on inflow patch'."""
    _, _, _, files = written
    for name in ("boundaryT", "boundaryU", "boundaryNumberDensity_Ar"):
        block = files[name].split("    nozzle\n    {", 1)[1].split("}", 1)[0]
        assert "fixedValue" in block


def test_the_vacuum_boundaries_carry_nothing(written):
    """They are not inflows -- plumeFieldInflow has an explicit patch list -- so
    zeroGradient says plainly that what they hold is irrelevant."""
    _, _, _, files = written
    block = files["boundaryT"].split("    vacuum\n    {", 1)[1].split("}", 1)[0]
    assert "zeroGradient" in block
    assert "fixedValue" not in block


def test_the_exit_density_changes_with_the_knudsen_number(tmp_path):
    """The field the solver injects from cannot drift from the Kn that made it."""
    densities = set()
    for kn in (100.0, 0.1, 0.01):
        directory = tmp_path / str(kn)
        cfg, _, state, _, _ = derived(directory, exit={"knudsen": kn})
        foamfields.write_all(directory, cfg, state, PATCHES)
        text = (directory / "0" / "boundaryNumberDensity_Ar").read_text(
            encoding="utf-8")
        block = text.split("    nozzle\n    {", 1)[1].split("}", 1)[0]
        densities.add(block.split("uniform", 1)[1].split(";", 1)[0].strip())
    assert len(densities) == 3


def test_a_symmetry_patch_gets_a_symmetry_patch_field(tmp_path):
    """OpenFOAM rejects anything else with 'inconsistent patch and patchField
    types'."""
    cfg, _, state, _, _ = derived(tmp_path, geometry={"symmetry_mode": "half_y"})
    foamfields.write_all(tmp_path, cfg, state, HALF_PATCHES)
    text = (tmp_path / "0" / "boundaryT").read_text(encoding="utf-8")
    block = text.split("    symmetry\n    {", 1)[1].split("}", 1)[0]
    assert "type            symmetry;" in block


def test_a_missing_nozzle_patch_is_an_error(tmp_path):
    cfg, _, state, _, _ = derived(tmp_path)
    with pytest.raises(ValueError, match="not in the mesh"):
        foamfields.write_all(tmp_path, cfg, state,
                             {"vacuum": PatchInfo("vacuum", "patch", 10, 0)})


def test_an_empty_nozzle_patch_is_an_error(tmp_path):
    """Nothing would ever be injected."""
    cfg, _, state, _, _ = derived(tmp_path)
    patches = dict(PATCHES, nozzle=PatchInfo("nozzle", "patch", 0, 0))
    with pytest.raises(ValueError, match="zero faces"):
        foamfields.write_all(tmp_path, cfg, state, patches)


def test_the_measurement_fields_name_this_meshs_patches(written):
    """Generated rather than shipped in a 0.orig, because 'nozzle' does not
    exist until createPatch has run."""
    _, _, _, files = written
    for name, *_ in foamfields.MEASUREMENT_FIELDS:
        for patch in PATCHES:
            assert patch in files[name], (name, patch)


def test_fD_and_q_are_calculated_on_a_wall_if_one_ever_exists(tmp_path):
    """This case has no walls, but a variant that adds one must not lose its
    surface data: zeroGradient writes no values at all."""
    cfg, _, state, _, _ = derived(tmp_path)
    patches = dict(PATCHES, body=PatchInfo("body", "wall", 12, 80000))
    foamfields.write_all(tmp_path, cfg, state, patches)
    text = (tmp_path / "0" / "fD").read_text(encoding="utf-8")
    block = text.split("    body\n    {", 1)[1].split("}", 1)[0]
    assert "calculated" in block


def test_written_fields_use_lf_line_endings(written, tmp_path):
    for path in (tmp_path / "0").iterdir():
        assert b"\r\n" not in path.read_bytes()


# --------------------------------------------------------------------------- #
# constant/dsmcProperties
# --------------------------------------------------------------------------- #

@pytest.fixture
def properties(tmp_path):
    cfg, _, state, _, run = derived(tmp_path)
    return cfg, state, run, dictionaries.render_dsmc_properties(cfg, state, run)


def test_dsmc_properties_selects_vhs(properties):
    """[PAPER] Cai states Variable Hard Sphere."""
    _, _, _, text = properties
    assert "BinaryCollisionModel            VariableHardSphere;" in text
    assert "VariableHardSphereCoeffs" in text


def test_dsmc_properties_carries_the_reference_temperature(properties):
    cfg, _, _, text = properties
    assert f"Tref                        {cfg.gas.t_ref_K:g};" in text


def test_dsmc_properties_uses_the_repository_argon(properties):
    _, _, _, text = properties
    assert "mass                            6.63e-26;" in text
    assert "diameter                        4.17e-10;" in text
    assert "omega                           0.74;" in text
    assert "internalDegreesOfFreedom        0;" in text


def test_dsmc_properties_writes_the_derived_particle_weight(properties):
    _, _, run, text = properties
    assert f"nEquivalentParticles            " \
           f"{run.n_equivalent_particles:.6e};" in text


def test_dsmc_properties_injects_only_on_the_nozzle(properties):
    """The single most important line in the case."""
    _, _, _, text = properties
    assert "InflowBoundaryModel             plumeFieldInflow;" in text
    assert "patches ( nozzle );" in text


def test_dsmc_properties_names_the_per_face_density_field(properties):
    _, _, _, text = properties
    assert "Ar    boundaryNumberDensity_Ar;" in text


def test_dsmc_properties_declares_one_species(properties):
    _, _, _, text = properties
    assert "typeIdList                      (Ar);" in text


def test_collisions_can_be_switched_off_but_the_file_says_so(tmp_path):
    """A NoBinaryCollision run cannot be evidence that DSMC APPROACHES the
    collisionless solution -- it is the collisionless solution."""
    cfg, _, state, _, run = derived(tmp_path, dsmc={"collisions_enabled": False})
    text = dictionaries.render_dsmc_properties(cfg, state, run)
    assert "BinaryCollisionModel            NoBinaryCollision;" in text
    assert "COLLISIONLESS" in text


def test_a_non_positive_particle_weight_is_rejected(tmp_path):
    import dataclasses

    cfg, _, state, _, run = derived(tmp_path)
    with pytest.raises(ValueError, match="nEquivalentParticles"):
        dictionaries.render_dsmc_properties(
            cfg, state, dataclasses.replace(run, n_equivalent_particles=0.0))


# --------------------------------------------------------------------------- #
# system/controlDict
# --------------------------------------------------------------------------- #

def test_control_dict_loads_the_custom_inflow_library(tmp_path):
    """Without this line dsmcFoam reports plumeFieldInflow as an unknown model
    and stops -- which is the intended failure, not a fallback."""
    cfg, _, _, _, run = derived(tmp_path)
    text = dictionaries.render_control_dict(cfg, run)
    assert f'libs            ( "{dictionaries.PLUME_LIB}" );' in text


def test_control_dict_carries_the_derived_times(tmp_path):
    cfg, _, _, _, run = derived(tmp_path)
    text = dictionaries.render_control_dict(cfg, run)
    assert f"deltaT          {run.delta_t_s:g};" in text
    assert f"endTime         {run.end_time_s:g};" in text
    assert f"timeStart       {run.average_start_s:g};" in text


def test_control_dict_averages_rhoN_first(tmp_path):
    """A single-timestep rhoN is one step's worth of parcels. rhoNMean is what
    the validation reads."""
    cfg, _, _, _, run = derived(tmp_path)
    text = dictionaries.render_control_dict(cfg, run)
    fields = text.split("fields\n        (", 1)[1]
    assert fields.strip().startswith("rhoN")
    for name in ("rhoM", "momentum", "linearKE"):
        assert name in fields


def test_control_dict_records_cais_time_step_as_a_comment(tmp_path):
    cfg, _, _, _, run = derived(tmp_path)
    assert "Cai's dt/t0 = 1" in dictionaries.render_control_dict(cfg, run)


def test_control_dict_resumes_from_the_latest_time(tmp_path):
    """`startFrom latestTime` is what makes ./Allrun continue a case instead of
    restarting it. On a case with no results it resolves to 0, so a fresh run is
    unaffected -- and `purgeWrite 0` is what keeps something to resume from."""
    cfg, _, _, _, run = derived(tmp_path)
    text = dictionaries.render_control_dict(cfg, run)
    assert "startFrom       latestTime;" in text
    assert "startTime       0;" in text
    assert "purgeWrite      0;" in text


def test_control_dict_does_not_adjust_the_time_step(tmp_path):
    """DSMC sampling assumes a fixed step."""
    cfg, _, _, _, run = derived(tmp_path)
    assert "adjustTimeStep  no;" in dictionaries.render_control_dict(cfg, run)


# --------------------------------------------------------------------------- #
# the rest
# --------------------------------------------------------------------------- #

def test_initialise_dict_is_a_near_vacuum_at_rest(tmp_path):
    cfg, _, _, _, _ = derived(tmp_path)
    text = dictionaries.render_dsmc_initialise_dict(cfg)
    assert "velocity        (0 0 0);" in text
    assert "Ar          1e+10;" in text


@pytest.mark.parametrize("n,expected", [
    (1, (1, 1, 1)),
    (2, (1, 2, 1)),
    (4, (1, 2, 2)),
    (6, (1, 3, 2)),
    (8, (2, 2, 2)),
    (12, (2, 3, 2)),
    (7, (1, 7, 1)),
])
def test_the_subdomain_split_prefers_lateral_cuts(n, expected):
    """(n 1 1) hands the whole graded outer shell -- where most of the parcels
    are -- to the last rank. Measured 3x imbalance at Kn = 100."""
    assert dictionaries.simple_split(n) == expected


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 8, 9, 12, 16, 24, 32, 48, 64])
def test_the_subdomain_split_multiplies_back_to_n(n):
    nx, ny, nz = dictionaries.simple_split(n)
    assert nx * ny * nz == n


def test_a_non_positive_subdomain_count_is_rejected():
    with pytest.raises(ValueError, match="n_subdomains"):
        dictionaries.simple_split(0)


def test_decompose_par_dict_writes_the_three_dimensional_split(tmp_path):
    cfg, _, _, _, _ = derived(tmp_path, dsmc={"n_subdomains": 12})
    text = dictionaries.render_decompose_par_dict(cfg)
    assert "numberOfSubdomains  12;" in text
    assert "method              simple;" in text
    assert "n               (2 3 2);" in text


def test_fv_schemes_and_solution_are_empty_but_present(tmp_path):
    """fvMesh requires the files; dsmcFoam solves no PDE."""
    assert "default none;" in dictionaries.render_fv_schemes()
    assert "solvers         { }" in dictionaries.render_fv_solution()


def test_write_case_dictionaries_writes_all_six(tmp_path):
    cfg, _, state, _, run = derived(tmp_path)
    paths = dictionaries.write_case_dictionaries(tmp_path, cfg, state, run)
    assert sorted(p.name for p in paths) == [
        "controlDict", "decomposeParDict", "dsmcInitialiseDict",
        "dsmcProperties", "fvSchemes", "fvSolution"]
    for path in paths:
        assert path.is_file()
        assert b"\r\n" not in path.read_bytes()
        assert dictionaries.GENERATOR in path.read_text(encoding="utf-8")
