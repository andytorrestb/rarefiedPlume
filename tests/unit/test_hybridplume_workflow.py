"""hybridplume: configuration, generated dictionaries and refusals."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from plumetools.hybridplume import cfd, config
from plumetools.hybridplume import mesh as hm

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases" / "cai2012-hybrid" / "hybrid.yaml"


def test_every_profile_validates_and_is_full_3d():
    raw = yaml.safe_load(STUDY.read_text())
    for name in ["base"] + list(raw["profiles"]):
        cfg = config.load(STUDY, name)
        src = config.build_source(cfg)
        plan = hm.plan(cfg, src)
        text = hm.render_block_mesh_dict(plan)
        for forbidden in ("wedge", "empty", "symmetry"):
            assert f"type {forbidden}" not in text
        assert plan.divisions[1] == plan.divisions[2]   # no half domain


@pytest.mark.parametrize("key", ["mesh.axisymmetric", "mesh.wedge", "dsmc.radial_weighting",
                                 "domain.symmetry"])
def test_symmetry_reductions_are_refused(key):
    section, item = key.split(".")
    with pytest.raises(config.ConfigError):
        config.load(STUDY, "smoke-kn0p01", {section: {item: True}})


def test_background_floor_must_be_a_floor():
    with pytest.raises(config.ConfigError):
        config.load(STUDY, "smoke-kn0p01", {"cfd": {"background_pressure_ratio": 0.5}})


def test_departures_from_paper_are_labelled():
    cfg = config.load(STUDY, "smoke-kn1em3")
    dep = config.departures_from_paper(cfg)
    assert dep["source.knudsen"]["label"].startswith("extension")


def test_cfd_case_uses_vhs_transport_and_caps_the_viscous_step(tmp_path):
    cfg = config.load(STUDY, "smoke-kn0p01")
    src = config.build_source(cfg)
    plan = hm.plan(cfg, src)
    out = cfd.write_case(tmp_path, cfg, src, plan)
    assert out["transport_table_max_rel_error"] < 1e-4
    thermo = (tmp_path / "constant" / "thermophysicalProperties").read_text()
    assert "transport       tabulated;" in thermo and "sutherland" not in thermo
    ctrl = (tmp_path / "system" / "controlDict").read_text()
    assert 'libs ("libhybridPlumeThermo.so");' in ctrl
    tm = out["timing"]
    assert tm["max_delta_t_s"] <= tm["viscous_dt_limit_s"] * (1 + 1e-12)
    p = (tmp_path / "0" / "p").read_text()
    assert "waveTransmissive" in p and f"{src.pressure:.10g}" in p
