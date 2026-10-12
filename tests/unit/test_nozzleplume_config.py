"""``case.yaml`` loading and validation for the nozzle-plume family.

Two things are being protected here:

* the shipped ``cases/r-4d/baseCase/case.yaml`` loads, and carries the R-4D
  exit-plane data it was set up from -- so a schema change that breaks the real
  case, or an edit that moves one of its given values, fails a Tier-0 test
  rather than a mesh run;
* every cross-check that would otherwise surface as a meshing failure, a solver
  abort, or (worse) a case that runs cleanly at the wrong conditions.
"""

from __future__ import annotations

import pytest
import yaml

from _nozzleplume import MINIMAL_CASE, R4D_CASE, load_case, write_case
from plumetools.nozzleplume import MODEL
from plumetools.nozzleplume.config import ConfigError, load_case_config


# --------------------------------------------------------------------------- #
# the shipped case
# --------------------------------------------------------------------------- #

def test_the_shipped_r4d_case_loads():
    assert load_case_config(R4D_CASE).model == MODEL


def test_the_shipped_r4d_case_carries_the_given_exit_plane_data():
    """The values the case was set up from. These must not drift."""
    cfg = load_case_config(R4D_CASE)
    assert cfg.nozzle.exit_diameter_m == 0.376936
    assert cfg.exit.temperature_K == 498.04
    assert cfg.exit.speed_m_per_s == 3187.19
    assert cfg.exit.number_density_per_m3 == 1.306e22
    assert cfg.gas.molecular_mass_kg == 3.40055e-26
    assert cfg.gas.molecular_diameter_m == 3.5e-10


def test_the_shipped_r4d_case_carries_the_given_plate():
    cfg = load_case_config(R4D_CASE)
    assert (cfg.plate.size_y_m, cfg.plate.size_z_m, cfg.plate.thickness_m) == (
        34.14, 11.88, 0.25)


def test_the_shipped_r4d_case_pins_the_particle_weight():
    """The point of the case: nEquivalentParticles is edited directly."""
    cfg = load_case_config(R4D_CASE)
    assert cfg.dsmc.n_equivalent_particles is not None
    assert cfg.dsmc.n_equivalent_particles > 0.0


def test_the_shipped_r4d_case_labels_what_was_not_given():
    """The stand-off, the wall temperature and the gas model are choices, and
    the file has to say so where they are made."""
    text = (R4D_CASE / "case.yaml").read_text(encoding="utf-8")
    for line in ("standoff_m: 20.0", "temperature_K: 300.0", "omega: 0.5",
                 "internal_degrees_of_freedom: 0"):
        tagged = [row for row in text.splitlines() if line in row]
        assert tagged and "[ASSUMPTION]" in tagged[0], line


# --------------------------------------------------------------------------- #
# structure
# --------------------------------------------------------------------------- #

def test_the_minimal_case_loads(tmp_path):
    cfg = load_case(tmp_path)
    assert cfg.model == MODEL
    assert cfg.plate is not None
    assert cfg.species().name == "exhaust"


def test_a_case_directory_or_the_file_itself_both_load(tmp_path):
    path = write_case(tmp_path)
    assert load_case_config(path).nozzle == load_case_config(tmp_path).nozzle


def test_a_missing_file_is_reported(tmp_path):
    with pytest.raises(ConfigError, match="no case.yaml"):
        load_case_config(tmp_path)


def test_the_wrong_model_is_rejected(tmp_path):
    """A Cai case.yaml belongs to the other loader; the two share no keys."""
    with pytest.raises(ConfigError, match="cai2012_circular_plume"):
        load_case(tmp_path, model="cai2012_circular_plume")


def test_an_unknown_top_level_key_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="unknown top-level"):
        load_case(tmp_path, geometry={"x_max_m": 1.0})


def test_an_unknown_key_in_a_section_is_rejected(tmp_path):
    """A typo in a physical parameter must stop the run, not fall back to a
    default."""
    with pytest.raises(ConfigError, match="temperatur_K"):
        load_case(tmp_path, exit={"temperatur_K": 500.0})


@pytest.mark.parametrize("section", ["nozzle", "exit", "gas", "domain", "mesh"])
def test_a_missing_required_section_is_reported(tmp_path, section):
    with pytest.raises(ConfigError, match=f"required section.*{section}"):
        load_case(tmp_path, **{section: None})


def test_a_missing_required_key_is_named(tmp_path):
    data = {k: dict(v) if isinstance(v, dict) else v for k, v in MINIMAL_CASE.items()}
    del data["gas"]["molecular_mass_kg"]
    (tmp_path / "case.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(ConfigError, match="molecular_mass_kg"):
        load_case_config(tmp_path)


def test_a_numeric_string_is_coerced(tmp_path):
    """YAML 1.1 reads 1.0e22 -- no sign on the exponent -- as a string."""
    (tmp_path / "case.yaml").write_text(
        yaml.safe_dump(MINIMAL_CASE).replace("1.0e+22", "1.0e22"), encoding="utf-8")
    cfg = load_case_config(tmp_path)
    assert cfg.exit.number_density_per_m3 == 1.0e22


def test_a_non_numeric_string_names_the_exponent_trap(tmp_path):
    with pytest.raises(ConfigError, match="signed exponent"):
        load_case(tmp_path, exit={"temperature_K": "warm"})


def test_no_plate_section_means_a_free_plume(tmp_path):
    assert load_case(tmp_path, plate=None).plate is None


def test_a_patch_name_override_keeps_the_other_defaults(tmp_path):
    cfg = load_case(tmp_path, mesh={"patch_names": {"plate": "panel"}})
    assert cfg.mesh.patch_names["plate"] == "panel"
    assert cfg.mesh.patch_names["nozzle"] == "nozzle"


def test_an_unknown_patch_role_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="unknown patch role"):
        load_case(tmp_path, mesh={"patch_names": {"cylinder": "cylinder"}})


def test_two_roles_sharing_a_patch_name_are_rejected(tmp_path):
    """Merging two boundaries merges two physics."""
    with pytest.raises(ConfigError, match="more than one role"):
        load_case(tmp_path, mesh={"patch_names": {"plate": "vacuum"}})


# --------------------------------------------------------------------------- #
# exactly one of
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("extra", [{"mass_density_kg_per_m3": 3.4e-4},
                                   {"knudsen": 1e-3}])
def test_two_exit_densities_are_rejected(tmp_path, extra):
    with pytest.raises(ConfigError, match="exactly one"):
        load_case(tmp_path, exit=extra)


def test_two_exit_speeds_are_rejected(tmp_path):
    with pytest.raises(ConfigError, match="exactly one"):
        load_case(tmp_path, exit={"speed_ratio": 5.0})


def test_no_exit_density_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="none is set"):
        load_case(tmp_path, exit={"number_density_per_m3": None})


@pytest.mark.parametrize("exit_section", [
    {"number_density_per_m3": None, "mass_density_kg_per_m3": 3.4e-4},
    {"number_density_per_m3": None, "knudsen": 1e-3},
    {"speed_m_per_s": None, "speed_ratio": 5.0},
])
def test_each_alternative_form_is_accepted_alone(tmp_path, exit_section):
    load_case(tmp_path, exit=exit_section)


# --------------------------------------------------------------------------- #
# values
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("section,key", [
    ("nozzle", "exit_diameter_m"), ("exit", "temperature_K"),
    ("gas", "molecular_mass_kg"), ("gas", "molecular_diameter_m"),
    ("plate", "standoff_m"), ("plate", "thickness_m"),
    ("domain", "x_max_m"), ("mesh", "background_cell_size_m"),
])
def test_a_non_positive_physical_input_is_rejected(tmp_path, section, key):
    with pytest.raises(ConfigError, match=f"{section}.{key}"):
        load_case(tmp_path, **{section: {key: -1.0}})


def test_a_non_positive_particle_weight_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="n_equivalent_particles"):
        load_case(tmp_path, dsmc={"n_equivalent_particles": 0.0})


def test_a_null_particle_weight_is_accepted_and_means_derive(tmp_path):
    assert load_case(
        tmp_path, dsmc={"n_equivalent_particles": None}
    ).dsmc.n_equivalent_particles is None


@pytest.mark.parametrize("key,value", [
    ("collision_model", "HardSphere"),
    ("wall_interaction_model", "Sticky"),
    ("decomposition_method", "metis"),
])
def test_an_unknown_model_name_is_rejected(tmp_path, key, value):
    with pytest.raises(ConfigError, match=key):
        load_case(tmp_path, dsmc={key: value})


def test_a_diffuse_fraction_outside_zero_to_one_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="diffuse_fraction"):
        load_case(tmp_path, dsmc={"diffuse_fraction": 1.5})


def test_a_species_name_that_is_not_an_identifier_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="species_name"):
        load_case(tmp_path, gas={"species_name": "R-4D"})


def test_omega_below_a_hard_sphere_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="omega"):
        load_case(tmp_path, gas={"omega": 0.3})


# --------------------------------------------------------------------------- #
# symmetry
# --------------------------------------------------------------------------- #

def test_the_plume_axis_is_never_a_symmetry_plane(tmp_path):
    with pytest.raises(ConfigError, match="symmetry_planes"):
        load_case(tmp_path, domain={"symmetry_planes": ["x"]})


def test_a_repeated_symmetry_plane_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="repeats"):
        load_case(tmp_path, domain={"symmetry_planes": ["y", "y"]})


def test_an_off_axis_plate_cannot_be_modelled_with_its_symmetry_plane(tmp_path):
    """Halving the domain about y = 0 would mirror the offset plate."""
    with pytest.raises(ConfigError, match="centre_y_m"):
        load_case(tmp_path, domain={"symmetry_planes": ["y"]},
                  plate={"centre_y_m": 0.5})


def test_an_off_axis_plate_is_fine_about_the_other_plane(tmp_path):
    load_case(tmp_path, domain={"symmetry_planes": ["z"]},
              plate={"centre_y_m": 0.5})


# --------------------------------------------------------------------------- #
# refinement regions
# --------------------------------------------------------------------------- #

def _region(**changes) -> dict:
    return {"mesh": {"refinement_regions": [
        {"name": "r", "level": 2, "x_max_m": 1.0, "radius_m": 0.5, **changes}]}}


def test_a_region_defaults_to_a_cylinder_from_the_exit_plane(tmp_path):
    region = load_case(tmp_path, **_region()).mesh.refinement_regions[0]
    assert (region.shape, region.x_min_m) == ("cylinder", 0.0)


def test_a_cylinder_without_a_radius_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="radius_m"):
        load_case(tmp_path, **_region(radius_m=None))


def test_a_box_needs_all_four_lateral_bounds(tmp_path):
    with pytest.raises(ConfigError, match="z_max_m"):
        load_case(tmp_path, **_region(shape="box", radius_m=None, y_min_m=-1.0,
                                      y_max_m=1.0, z_min_m=-1.0))


def test_a_box_region_loads(tmp_path):
    region = load_case(tmp_path, **_region(
        shape="box", radius_m=None, y_min_m=-1.0, y_max_m=1.0, z_min_m=-1.0,
        z_max_m=1.0)).mesh.refinement_regions[0]
    assert region.shape == "box"


def test_a_box_given_a_radius_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="cylinder"):
        load_case(tmp_path, **_region(shape="box", y_min_m=-1.0, y_max_m=1.0,
                                      z_min_m=-1.0, z_max_m=1.0))


def test_level_zero_is_the_background_and_needs_no_region(tmp_path):
    with pytest.raises(ConfigError, match="level"):
        load_case(tmp_path, **_region(level=0))


def test_an_inverted_region_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match="x_max_m"):
        load_case(tmp_path, **_region(x_min_m=2.0))


def test_two_regions_with_one_name_are_rejected(tmp_path):
    region = {"name": "r", "level": 2, "x_max_m": 1.0, "radius_m": 0.5}
    with pytest.raises(ConfigError, match="repeats the name"):
        load_case(tmp_path, mesh={"refinement_regions": [region, dict(region)]})


def test_a_region_named_after_a_patch_is_rejected(tmp_path):
    """snappyHexMesh keeps both in one geometry dictionary."""
    with pytest.raises(ConfigError, match="also a patch name"):
        load_case(tmp_path, **_region(name="plate"))
