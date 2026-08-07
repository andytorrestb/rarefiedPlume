"""The generated OpenFOAM dictionaries and 0/ fields.

These are the files the solver actually reads, so the tests are about the
entries whose absence or wrong form would produce a case that runs and is wrong.
"""

from __future__ import annotations

import re
from dataclasses import replace

import numpy as np
import pytest

from plumetools.config import StagnationConfig
from plumetools.markelov1999 import foamdicts, foamfields
from plumetools.markelov1999.constants import PSI_TO_PA, molecular_mass_kg
from plumetools.mesh.boundary import PatchInfo

from test_markelov_mesh import make_cfg  # noqa: E402

WEIGHT = 3.622366595756453e10


@pytest.fixture
def cfg():
    base = make_cfg()
    return replace(
        base,
        stagnation=StagnationConfig(
            p0_pa=5 * PSI_TO_PA, T0_K=300.0, throat_radius_m=0.00041275),
        output=replace(base.output, inflow_model="plumeFieldInflow"),
    )


@pytest.fixture
def patches():
    return {
        "symmetry": PatchInfo("symmetry", "symmetry", 4, 0),
        "upstreamVacuum": PatchInfo("upstreamVacuum", "patch", 2, 4),
        "vacuum": PatchInfo("vacuum", "patch", 3, 6),
        "inflow": PatchInfo("inflow", "patch", 5, 9),
        "cylinder": PatchInfo("cylinder", "wall", 6, 14),
        "plate": PatchInfo("plate", "wall", 7, 20),
    }


class FakeInflow:
    def __init__(self, n=5):
        self.n_faces = n
        self.rhoN = np.linspace(2.2e17, 6.3e18, n)
        self.T = np.full(n, 300.0)
        self.U = np.tile([789.4849, 0.0, 0.0], (n, 1))
        self.v_limit = 789.4849
        self.areas = np.full(n, 4.6e-5)


# --------------------------------------------------------------------------- #
# constant/dsmcProperties
# --------------------------------------------------------------------------- #

def test_nitrogen_is_the_only_species(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    assert re.search(r"typeIdList\s+\(N2\);", text)
    assert re.search(r"moleculeProperties\s*\{\s*N2", text)
    assert "Ar" not in text, "the legacy case mixes an N2 model with an Ar solver"


def test_the_collision_model_family_the_paper_names(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    assert "BinaryCollisionModel            LarsenBorgnakkeVariableHardSphere;" in text
    assert re.search(r"relaxationCollisionNumber\s+5;", text), "ZR = 5 is a paper value"
    assert re.search(r"Tref\s+273;", text)


def test_two_rotational_degrees_of_freedom_and_no_vibration(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    assert re.search(r"internalDegreesOfFreedom\s+2;", text)
    assert "vibrational" not in text.lower()


def test_the_vhs_coefficients_are_exposed_and_sourced(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    assert re.search(r"diameter\s+4\.17e-10;", text)
    assert re.search(r"omega\s+0\.74;", text)
    assert "wedge15Ma5" in text, "the source of the non-paper coefficients"
    assert "ASSUMPTION" in text


def test_the_mass_is_derived_from_the_molar_mass(cfg):
    """So the analytical model and the solver cannot disagree about m. The
    tutorial's rounded 46.5e-27 is 0.04% away, and two masses that nearly agree
    are harder to debug than one."""
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    match = re.search(r"mass\s+([\d.e+-]+);", text)
    assert match
    assert float(match.group(1)) == pytest.approx(
        molecular_mass_kg(28.0134), rel=1e-9)
    assert float(match.group(1)) == pytest.approx(4.651735e-26, rel=1e-6)


def test_the_custom_inflow_model_is_selected_with_a_patch_list(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    assert "InflowBoundaryModel             plumeFieldInflow;" in text
    assert re.search(r"patches\s*\(\s*inflow\s*\);", text)
    assert re.search(r"N2\s+boundaryNumberDensity_N2;", text)


def test_free_stream_is_refused_rather_than_written(cfg):
    """It takes one number density per species and injects on every patch-type
    boundary. Writing FreeStreamCoeffs would produce a case that runs and is
    wrong, which is worse than one that does not run."""
    stock = replace(cfg, output=replace(cfg.output, inflow_model="FreeStream"))
    with pytest.raises(ValueError, match="cannot express this case"):
        foamdicts.render_dsmc_properties(stock, WEIGHT)


def test_the_particle_weight_is_written(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    match = re.search(r"nEquivalentParticles\s+([\d.e+-]+);", text)
    assert match
    assert float(match.group(1)) == pytest.approx(WEIGHT, rel=1e-6)


def test_a_nonpositive_weight_is_refused(cfg):
    with pytest.raises(ValueError, match="must be positive"):
        foamdicts.render_dsmc_properties(cfg, 0.0)


def test_the_wall_model_is_named_and_its_temperature_source_explained(cfg):
    text = foamdicts.render_dsmc_properties(cfg, WEIGHT)
    assert "WallInteractionModel            MaxwellianThermal;" in text
    assert "boundaryT" in text, "the wall temperature comes from there, not from here"


# --------------------------------------------------------------------------- #
# system/controlDict
# --------------------------------------------------------------------------- #

def test_the_library_is_loaded(cfg):
    """Without this, dsmcFoam reports plumeFieldInflow as an unknown model."""
    text = foamdicts.render_control_dict(cfg)
    assert 'libs            ( "libplumeDsmcBoundaryModels.so" );' in text


def test_the_solver_and_times_are_written(cfg):
    text = foamdicts.render_control_dict(cfg)
    assert "application     dsmcFoam;" in text
    assert re.search(r"deltaT\s+2e-07;", text)
    assert re.search(r"endTime\s+0\.004;", text)


def test_field_average_covers_the_surface_force_density(cfg):
    """fD is reset every timestep by DSMCCloud::resetFields, so a snapshot is a
    one-step momentum tally. fDMean is what the pressure is read from."""
    text = foamdicts.render_control_dict(cfg)
    assert "type            fieldAverage;" in text
    block = text[text.index("fieldAverage1"):]
    assert re.search(r"^\s+fD$", block, re.M)
    assert re.search(r"^\s+dsmcRhoN$", block, re.M)
    assert re.search(r"timeStart\s+0\.002;", text)


def test_averaging_starts_after_the_transient(cfg):
    assert cfg.dsmc.average_start_s < cfg.dsmc.end_time_s
    text = foamdicts.render_control_dict(cfg)
    start = float(re.search(r"timeStart\s+([\d.e+-]+);", text).group(1))
    end = float(re.search(r"endTime\s+([\d.e+-]+);", text).group(1))
    assert 0 < start < end


# --------------------------------------------------------------------------- #
# system/dsmcInitialiseDict
# --------------------------------------------------------------------------- #

def test_the_domain_starts_near_vacuum_and_at_rest(cfg):
    text = foamdicts.render_dsmc_initialise_dict(cfg)
    assert re.search(r"N2\s+1e\+14;", text)
    assert "velocity        (0 0 0);" in text
    assert "near-vacuum" in text.lower()


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #

def test_all_dictionaries_are_written(cfg, tmp_path):
    paths = foamdicts.write_case_dictionaries(tmp_path, cfg, WEIGHT)
    names = {p.name for p in paths}
    assert names == {"dsmcProperties", "controlDict", "dsmcInitialiseDict",
                     "decomposeParDict", "fvSchemes", "fvSolution"}
    assert (tmp_path / "constant" / "dsmcProperties").is_file()
    assert (tmp_path / "system" / "controlDict").is_file()


def test_written_dictionaries_are_lf_terminated_and_say_they_are_generated(
        cfg, tmp_path):
    for path in foamdicts.write_case_dictionaries(tmp_path, cfg, WEIGHT):
        assert b"\r\n" not in path.read_bytes(), path.name
        assert "GENERATED by" in path.read_text(encoding="utf-8")
        assert "do not edit" in path.read_text(encoding="utf-8")


# --------------------------------------------------------------------------- #
# 0/ fields
# --------------------------------------------------------------------------- #

def test_boundary_t_is_a_scalar_field_not_a_vector(cfg, tmp_path, patches):
    """Standard dsmcFoam declares boundaryT as a volScalarField; the MNF fork
    takes a vector holding per-component translational temperature. Writing the
    wrong class is a hard read failure."""
    foamfields.write_inflow_fields(tmp_path, FakeInflow(), cfg, patches)
    text = (tmp_path / "0" / "boundaryT").read_text(encoding="utf-8")
    assert "class       volScalarField;" in text
    assert "volVectorField" not in text


def test_the_inflow_patch_carries_a_nonuniform_list(cfg, tmp_path, patches):
    foamfields.write_inflow_fields(tmp_path, FakeInflow(n=5), cfg, patches)
    for name, kind in (("boundaryU", "vector"),
                       ("boundaryT", "scalar"),
                       ("boundaryNumberDensity_N2", "scalar")):
        text = (tmp_path / "0" / name).read_text(encoding="utf-8")
        assert f"nonuniform List<{kind}>" in text, name
        block = text[text.index("inflow"):]
        assert re.search(rf"List<{kind}>\s*\n\s*5\s*\n", block), name


def test_walls_carry_a_temperature_because_the_wall_model_reads_it(
        cfg, tmp_path, patches):
    """MaxwellianThermal takes the wall temperature from cloud.boundaryT() on the
    wall face. A zeroGradient there evaluates to the zero internal field, so
    reflected particles would leave with no thermal speed at all."""
    foamfields.write_inflow_fields(tmp_path, FakeInflow(), cfg, patches)
    text = (tmp_path / "0" / "boundaryT").read_text(encoding="utf-8")

    for wall in ("cylinder", "plate"):
        block = text[text.index(f"    {wall}\n"):]
        block = block[:block.index("}")]
        assert "fixedValue" in block, wall
        assert "uniform 300" in block, wall


def test_walls_carry_a_zero_velocity(cfg, tmp_path, patches):
    foamfields.write_inflow_fields(tmp_path, FakeInflow(), cfg, patches)
    text = (tmp_path / "0" / "boundaryU").read_text(encoding="utf-8")
    block = text[text.index("    cylinder\n"):]
    assert "uniform (0 0 0)" in block[:block.index("}")]


def test_constraint_patches_repeat_their_own_type(cfg, tmp_path, patches):
    """OpenFOAM rejects anything else with "inconsistent patch and patchField
    types"."""
    foamfields.write_inflow_fields(tmp_path, FakeInflow(), cfg, patches)
    text = (tmp_path / "0" / "boundaryT").read_text(encoding="utf-8")
    block = text[text.index("    symmetry\n"):]
    assert "type            symmetry;" in block[:block.index("}")]


def test_fd_and_q_are_calculated_on_walls(cfg, tmp_path, patches):
    """DSMCParcel::hitWallPatch writes into their boundary values. A zeroGradient
    patch field writes no values at all, so the surface pressure -- which IS that
    data -- would be discarded at write time."""
    foamfields.write_measurement_fields(tmp_path, patches)

    for field in ("fD", "q"):
        text = (tmp_path / "0" / field).read_text(encoding="utf-8")
        for wall in ("cylinder", "plate"):
            block = text[text.index(f"    {wall}\n"):]
            block = block[:block.index("}")]
            assert "calculated" in block, f"{field} on {wall}"
        # An open boundary has nothing to record, so zeroGradient is right there.
        block = text[text.index("    vacuum\n"):]
        assert "zeroGradient" in block[:block.index("}")]


def test_the_other_measurement_fields_stay_zero_gradient(cfg, tmp_path, patches):
    """Only fD and q are written by hitWallPatch; the rest are cell fields."""
    foamfields.write_measurement_fields(tmp_path, patches)
    text = (tmp_path / "0" / "rhoN").read_text(encoding="utf-8")
    block = text[text.index("    cylinder\n"):]
    assert "zeroGradient" in block[:block.index("}")]


def test_all_nine_measurement_fields_are_written(cfg, tmp_path, patches):
    """dsmcFoam aborts with `cannot find file "0/q"` if any is missing, and
    dsmcInitialise does not create them."""
    written = foamfields.write_measurement_fields(tmp_path, patches)
    assert {p.name for p in written} == {
        "dsmcRhoN", "fD", "iDof", "internalE", "linearKE", "momentum", "q",
        "rhoM", "rhoN"}


def test_the_boundary_entries_name_this_meshs_patches(cfg, tmp_path, patches):
    """Generated, not shipped as a 0.orig: the legacy code hard-coded a patch
    list and wrote cylinder and plate into meshes that had neither (AD-03)."""
    reduced = {k: v for k, v in patches.items() if k != "plate"}
    foamfields.write_measurement_fields(tmp_path, reduced)
    text = (tmp_path / "0" / "rhoN").read_text(encoding="utf-8")
    assert "plate" not in text
    assert "cylinder" in text


def test_a_face_count_mismatch_is_an_error(cfg, tmp_path, patches):
    with pytest.raises(ValueError, match="but patch 'inflow' declares"):
        foamfields.write_inflow_fields(tmp_path, FakeInflow(n=99), cfg, patches)


def test_an_absent_inflow_patch_is_an_error(cfg, tmp_path, patches):
    del patches["inflow"]
    with pytest.raises(ValueError, match="not in the mesh"):
        foamfields.write_inflow_fields(tmp_path, FakeInflow(), cfg, patches)


def test_written_fields_are_lf_terminated(cfg, tmp_path, patches):
    for path in foamfields.write_all(tmp_path, FakeInflow(), cfg, patches):
        assert b"\r\n" not in path.read_bytes(), path.name


def test_the_density_written_round_trips(cfg, tmp_path, patches):
    """Precision matters: the flux verification compares the analytical value
    against the values as written, so formatting loss shows up there."""
    inflow = FakeInflow(n=5)
    foamfields.write_inflow_fields(tmp_path, inflow, cfg, patches)
    values = foamfields.read_patch_field(
        tmp_path / "0" / "boundaryNumberDensity_N2", "inflow")
    np.testing.assert_allclose(values, inflow.rhoN, rtol=1e-9)


def test_the_velocity_written_round_trips(cfg, tmp_path, patches):
    inflow = FakeInflow(n=5)
    foamfields.write_inflow_fields(tmp_path, inflow, cfg, patches)
    values = foamfields.read_patch_field(tmp_path / "0" / "boundaryU", "inflow")
    np.testing.assert_allclose(values, inflow.U, rtol=1e-9)
