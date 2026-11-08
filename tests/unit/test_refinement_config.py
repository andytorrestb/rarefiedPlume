from __future__ import annotations

import importlib
import pytest
import yaml
from _refinement import write_cartesian_case


def loader(path):
    return importlib.import_module("plumetools.refinement.config").load_settings(path)


def test_custom_gas_and_fields_do_not_require_nozzle_inputs(tmp_path):
    case = write_cartesian_case(tmp_path / "outside study")
    config = loader(case / "refinement.yaml")
    assert config.species.name == "Xe"
    assert config.fields["temperature"] == "TObserved"
    assert config.targets.edge_over_mfp == 0.8
    assert config.targets.min_mean_parcels == 20


@pytest.mark.parametrize("section,key,value", [(None, "nozzle", {}), ("gas", "diametre", 1),
                                              ("fields", "temperatur", "x"), ("targets", "mfp_ratio", 1)])
def test_unknown_keys_are_errors(tmp_path, section, key, value):
    case = write_cartesian_case(tmp_path / "case")
    path = case / "refinement.yaml"
    raw = yaml.safe_load(path.read_text())
    (raw if section is None else raw.setdefault(section, {}))[key] = value
    path.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match=key):
        loader(path)


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_target_is_rejected(tmp_path, value):
    case = write_cartesian_case(tmp_path / "case")
    path = case / "refinement.yaml"
    raw = yaml.safe_load(path.read_text())
    raw["targets"] = {"edge_over_mfp": value}
    path.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="edge_over_mfp"):
        loader(path)
