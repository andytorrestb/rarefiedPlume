"""case.yaml loading and validation."""

from __future__ import annotations

import textwrap

import pytest
import yaml

from _helpers import SQUARE_POINTS, TRI_FACES, make_case, patch
from plumetools.config import ConfigError, ConfigWarning, load_case_config

MINIMAL = {
    "model": "source_flow",
    "mesh": {"sphere_radius_m": 0.5},
    "geometry": {"patch": "inflow", "sphere_radius_m": 0.5},
    "gas": {"gamma": 1.4, "molar_mass_g_per_mol": 28.0134, "species_name": "Ar"},
    "stagnation": {"p0_pa": 3275009.71275, "T0_K": 300.0},
}


def write_case(tmp_path, config, *, with_mesh=True, dsmc=True):
    root = tmp_path / "case"
    root.mkdir(exist_ok=True)
    (root / "case.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    if with_mesh:
        make_case(root, points=SQUARE_POINTS, faces=TRI_FACES,
                  patches=[patch("inflow", 2, 0), patch("vacuum", 0, 2)])
    if dsmc:
        (root / "constant").mkdir(parents=True, exist_ok=True)
        (root / "constant" / "dsmcProperties").write_text(
            textwrap.dedent("""\
                typeIdList                      (Ar);
                moleculeProperties { Ar { mass 6.63e-26; } }
                """), encoding="utf-8")
    return root


def test_loads_a_valid_config(tmp_path):
    cfg = load_case_config(write_case(tmp_path, MINIMAL))
    assert cfg.gas.gamma == 1.4
    assert cfg.stagnation.p0_pa == 3275009.71275


def test_defaults_fill_in_omitted_sections(tmp_path):
    cfg = load_case_config(write_case(tmp_path, MINIMAL))
    assert cfg.angular.exponent_offset == 0.41
    assert cfg.angular.quadrature_points == 500


def test_missing_file_raises(tmp_path):
    with pytest.raises(ConfigError, match="no case.yaml"):
        load_case_config(tmp_path)


def test_unknown_top_level_key_is_rejected(tmp_path):
    cfg = {**MINIMAL, "sourceflow": {}}
    with pytest.raises(ConfigError, match="unknown top-level key"):
        load_case_config(write_case(tmp_path, cfg))


def test_unknown_section_key_is_rejected(tmp_path):
    """A typo in a physical parameter must fail, not silently take the default."""
    cfg = {**MINIMAL, "stagnation": {**MINIMAL["stagnation"], "p0_psi": 475}}
    with pytest.raises(ConfigError, match="p0_psi"):
        load_case_config(write_case(tmp_path, cfg))


def test_mismatched_sphere_radii_are_rejected(tmp_path):
    """The generator and the model must mesh and evaluate the same sphere (SM-09)."""
    cfg = {**MINIMAL, "geometry": {"patch": "inflow", "sphere_radius_m": 0.4}}
    with pytest.raises(ConfigError, match="sphere_radius_m"):
        load_case_config(write_case(tmp_path, cfg))


def test_inflow_patch_must_exist_in_the_mesh(tmp_path):
    cfg = {**MINIMAL, "geometry": {"patch": "nope", "sphere_radius_m": 0.5}}
    with pytest.raises(ConfigError, match="not in the mesh"):
        load_case_config(write_case(tmp_path, cfg))


def test_species_must_be_in_the_dsmc_type_list(tmp_path):
    """Catches writing boundaryNumberDensity_<X> for a species the solver lacks."""
    cfg = {**MINIMAL, "gas": {**MINIMAL["gas"], "species_name": "He"}}
    with pytest.raises(ConfigError, match="typeIdList"):
        load_case_config(write_case(tmp_path, cfg))


def test_absent_output_patch_warns_when_legacy_flag_is_set(tmp_path):
    """AD-03 frozen: the original wrote cylinder/plate into every field file."""
    cfg = {**MINIMAL,
           "output": {"patches": {"cylinder": {"type": "calculated", "value": "uniform"}}},
           "legacy": {"emit_absent_patches": True}}
    with pytest.warns(ConfigWarning, match="cylinder"):
        load_case_config(write_case(tmp_path, cfg))


def test_absent_output_patch_is_an_error_without_the_legacy_flag(tmp_path):
    """New cases should not be able to name a patch the mesh does not have."""
    cfg = {**MINIMAL,
           "output": {"patches": {"cylinder": {"type": "calculated"}}},
           "legacy": {"emit_absent_patches": False}}
    with pytest.raises(ConfigError, match="not in the mesh"):
        load_case_config(write_case(tmp_path, cfg))


def test_validation_is_skipped_before_the_mesh_exists(tmp_path):
    """A config must load so the mesh can be generated from it."""
    cfg = load_case_config(write_case(tmp_path, MINIMAL, with_mesh=False, dsmc=False))
    assert cfg.mesh.sphere_radius_m == 0.5


def test_rhoN_T0_can_be_unset_to_use_one_temperature(tmp_path):
    """`null` means 'use stagnation.T0_K', i.e. the corrected SM-02 behaviour."""
    cfg = {**MINIMAL, "legacy": {"rhoN_T0_K": None}}
    assert load_case_config(write_case(tmp_path, cfg)).legacy.rhoN_T0_K is None
