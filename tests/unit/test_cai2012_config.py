"""``case.yaml`` loading and validation for the Cai 2012 family.

Two things are being protected here:

* the shipped ``cases/cai2012/baseCase/case.yaml`` loads, and carries Cai's
  values -- so a schema change that breaks the real case fails a Tier-0 test
  rather than a mesh run;
* every cross-check that would otherwise surface as a meshing failure, a solver
  abort, or (worse) a case that runs cleanly at the wrong conditions.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from _cai2012 import MINIMAL_CASE, load_case, write_case
from plumetools.cai2012 import MODEL
from plumetools.cai2012.config import CaiConfigError, load_case_config

REPO = Path(__file__).resolve().parents[2]
BASE_CASE = REPO / "cases" / "cai2012" / "baseCase"


# --------------------------------------------------------------------------- #
# the shipped case
# --------------------------------------------------------------------------- #

def test_the_shipped_base_case_loads():
    cfg = load_case_config(BASE_CASE)
    assert cfg.model == MODEL


def test_the_shipped_base_case_carries_cais_paper_values():
    """D, S0 and the gas are printed in the paper; these must not drift."""
    cfg = load_case_config(BASE_CASE)
    assert cfg.nozzle.diameter_m == 0.2
    assert cfg.exit.speed_ratio == 2.0
    assert cfg.gas.species_name == "Ar"
    assert cfg.dsmc.binary_collision_model == "VariableHardSphere"
    assert cfg.exit.characteristic_length == "diameter"


def test_the_shipped_base_case_uses_the_repository_argon():
    cfg = load_case_config(BASE_CASE)
    assert (cfg.gas.mass_kg, cfg.gas.diameter_m, cfg.gas.omega) == (
        6.63e-26, 4.17e-10, 0.74)


def test_the_shipped_base_case_documents_T0_as_an_assumption():
    """300 K is an implementation assumption, and the file must say so."""
    text = (BASE_CASE / "case.yaml").read_text(encoding="utf-8")
    assert "T0_K: 300.0" in text
    assert "[ASSUMPTION]" in text


def test_the_shipped_base_case_keeps_collisions_on():
    """Even at Kn = 100. See docs/cai2012-case.md."""
    assert load_case_config(BASE_CASE).dsmc.collisions_enabled is True


def test_the_shipped_base_case_derives_rather_than_pins():
    """Cell size, particle weight and time step are all null -> derived."""
    cfg = load_case_config(BASE_CASE)
    assert cfg.mesh.core_cell_size_m is None
    assert cfg.dsmc.n_equivalent_particles is None
    assert cfg.dsmc.delta_t_s is None


def test_the_shipped_base_case_defaults_to_cais_transient():
    """The template says what the paper says; the study may override it."""
    cfg = load_case_config(BASE_CASE)
    assert cfg.dsmc.transient_basis == "cai"
    assert cfg.dsmc.transient_collision_times == 10000.0


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #

def test_missing_file_is_an_error(tmp_path):
    with pytest.raises(CaiConfigError, match="no case.yaml"):
        load_case_config(tmp_path)


def test_a_case_yaml_path_is_accepted_as_well_as_a_directory(tmp_path):
    path = write_case(tmp_path)
    assert load_case_config(path).model == MODEL


def test_the_wrong_model_is_rejected(tmp_path):
    """A source-flow case.yaml must not load here, and vice versa."""
    write_case(tmp_path)
    data = yaml.safe_load((tmp_path / "case.yaml").read_text(encoding="utf-8"))
    data["model"] = "markelov1999_axisymmetric"
    (tmp_path / "case.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(CaiConfigError, match="model is"):
        load_case_config(tmp_path)


def test_unknown_top_level_key_is_rejected(tmp_path):
    write_case(tmp_path)
    data = yaml.safe_load((tmp_path / "case.yaml").read_text(encoding="utf-8"))
    data["stagnation"] = {"p0_pa": 1.0}
    (tmp_path / "case.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(CaiConfigError, match="unknown top-level key"):
        load_case_config(tmp_path)


def test_unknown_section_key_is_rejected(tmp_path):
    """A typo in a physical parameter must fail, not fall back to a default."""
    with pytest.raises(CaiConfigError, match="unknown key"):
        load_case(tmp_path, exit={"speed_rato": 2.0})


def test_yaml_unsigned_exponent_strings_are_coerced(tmp_path):
    """YAML 1.1 parses 1.0e14 as a STRING; 1.0e+14 is a float.

    Without coercion the value looks right in the file, loads without complaint,
    and fails deep inside a dictionary renderer.
    """
    cfg = load_case(tmp_path, dsmc={"initial_number_density_per_m3": "1.0e10"})
    assert cfg.dsmc.initial_number_density_per_m3 == pytest.approx(1.0e10)


def test_a_non_numeric_string_on_a_numeric_key_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="not a number"):
        load_case(tmp_path, nozzle={"diameter_m": "wide"})


# --------------------------------------------------------------------------- #
# derived quantities on the config
# --------------------------------------------------------------------------- #

def test_radius_is_derived_from_the_diameter(tmp_path):
    cfg = load_case(tmp_path)
    assert cfg.nozzle.radius_m == pytest.approx(0.1)
    assert cfg.nozzle.area_m2 == pytest.approx(3.14159265e-2, rel=1e-6)


def test_characteristic_length_follows_the_configured_choice(tmp_path):
    assert load_case(tmp_path).characteristic_length_m == pytest.approx(0.2)
    cfg = load_case(tmp_path, exit={"characteristic_length": "radius"})
    assert cfg.characteristic_length_m == pytest.approx(0.1)


def test_species_round_trips_through_the_gas_section(tmp_path):
    species = load_case(tmp_path).species()
    assert (species.name, species.mass_kg, species.omega) == ("Ar", 6.63e-26, 0.74)


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("section,key", [
    ("nozzle", "diameter_m"), ("exit", "T0_K"), ("exit", "knudsen"),
    ("gas", "mass_kg"), ("geometry", "x_max_over_D"),
])
def test_non_positive_physical_inputs_are_rejected(tmp_path, section, key):
    with pytest.raises(CaiConfigError, match="finite and positive"):
        load_case(tmp_path, **{section: {key: 0.0}})


def test_a_negative_derived_key_is_rejected_but_null_is_not(tmp_path):
    load_case(tmp_path, dsmc={"delta_t_s": None})       # null means "derive"
    with pytest.raises(CaiConfigError, match="delta_t_s"):
        load_case(tmp_path, dsmc={"delta_t_s": -1.0e-6})


def test_unknown_characteristic_length_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="characteristic_length"):
        load_case(tmp_path, exit={"characteristic_length": "throat"})


def test_unknown_mean_free_path_convention_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="mean_free_path_convention"):
        load_case(tmp_path, gas={"mean_free_path_convention": "grasp"})


def test_unknown_symmetry_mode_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="symmetry_mode"):
        load_case(tmp_path, geometry={"symmetry_mode": "quarter"})


def test_half_y_without_a_symmetry_patch_is_rejected(tmp_path):
    """The y = 0 plane would be written open and delete every molecule."""
    with pytest.raises(CaiConfigError, match="symmetry"):
        load_case(tmp_path,
                  geometry={"symmetry_mode": "half_y"},
                  mesh={"patch_names": {"nozzle": "nozzle",
                                        "upstream_vacuum": "upstreamVacuum",
                                        "outer": "vacuum"}})


def test_half_y_with_a_symmetry_patch_is_accepted(tmp_path):
    cfg = load_case(tmp_path, geometry={"symmetry_mode": "half_y"})
    assert cfg.geometry.symmetry_mode == "half_y"


def test_the_mnf_dialect_is_rejected(tmp_path):
    """This family writes boundaryT as a volScalarField; the fork reads a vector."""
    with pytest.raises(CaiConfigError, match="dialect"):
        load_case(tmp_path, output={"dialect": "mnf"})


def test_free_stream_is_rejected_with_the_reason(tmp_path):
    """Stock FreeStream would turn every vacuum boundary into an inlet."""
    with pytest.raises(CaiConfigError, match="every patch-type boundary"):
        load_case(tmp_path, output={"inflow_model": "FreeStream"})


def test_a_patch_name_serving_two_roles_is_rejected(tmp_path):
    """The vacuum boundary would inject."""
    with pytest.raises(CaiConfigError, match="more than one role"):
        load_case(tmp_path, mesh={"patch_names": {
            "nozzle": "nozzle", "upstream_vacuum": "nozzle", "outer": "vacuum"}})


def test_a_missing_patch_role_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="missing role"):
        load_case(tmp_path, mesh={"patch_names": {"nozzle": "nozzle"}})


def test_internal_degrees_of_freedom_must_be_zero_for_argon(tmp_path):
    with pytest.raises(CaiConfigError, match="monatomic"):
        load_case(tmp_path, gas={"internal_degrees_of_freedom": 2})


def test_larsen_borgnakke_is_rejected_for_a_monatomic_gas(tmp_path):
    with pytest.raises(CaiConfigError, match="binary_collision_model"):
        load_case(tmp_path, dsmc={
            "binary_collision_model": "LarsenBorgnakkeVariableHardSphere"})


def test_a_core_narrower_than_the_nozzle_is_rejected(tmp_path):
    """The exit disk would straddle the grading boundary and inject through
    faces of two different sizes."""
    with pytest.raises(CaiConfigError, match="core_half_over_D"):
        load_case(tmp_path, mesh={"core_half_over_D": 0.4})


def test_a_core_reaching_past_the_domain_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="core_x_over_D"):
        load_case(tmp_path, mesh={"core_x_over_D": 20.0})


def test_a_centreline_reaching_past_the_domain_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="centerline_x_over_D_max"):
        load_case(tmp_path, post={"centerline_x_over_D_max": 20.0})


def test_a_contracting_outer_expansion_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="outer_expansion"):
        load_case(tmp_path, mesh={"outer_expansion": 0.5})


def test_a_core_cell_cap_coarser_than_a_quarter_diameter_is_rejected(tmp_path):
    """Fewer than four cells across the disk is a square, not a circle."""
    with pytest.raises(CaiConfigError, match="max_core_cell_over_D"):
        load_case(tmp_path, mesh={"max_core_cell_over_D": 0.5})


def test_warn_courant_above_max_courant_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="warn_courant"):
        load_case(tmp_path, checks={"warn_courant": 2.0, "max_courant": 1.0})


def test_an_unknown_post_plane_is_rejected(tmp_path):
    with pytest.raises(CaiConfigError, match="post.plane|plane"):
        load_case(tmp_path, post={"plane": "yz"})


def test_defaults_cover_every_section(tmp_path):
    """A case.yaml with only 'model' loads: every section has defaults."""
    (tmp_path / "case.yaml").write_text(f"model: {MODEL}\n", encoding="utf-8")
    cfg = load_case_config(tmp_path)
    assert cfg.nozzle.diameter_m == 0.2
    assert cfg.exit.speed_ratio == 2.0


def test_the_minimal_fixture_matches_the_schema_sections():
    """The test fixture must not drift into naming a section the schema lost."""
    from plumetools.cai2012.config import _SECTIONS
    assert set(MINIMAL_CASE) - {"model", "meta"} == set(_SECTIONS)
