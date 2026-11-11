"""Kinematic-viscosity cap for the coarse CFD profile; no solver."""
from pathlib import Path

import pytest

from plumetools.hybridplume import cfd, config, mesh

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases" / "cai2012-hybrid" / "hybrid.yaml"


def load(name, overlay=None):
    return config.load(STUDY, name, overlay)


def test_nu_max_equals_background_kinematic_viscosity():
    cfg = load("hpc-kn1em3-capsule")
    src = config.build_source(cfg)
    tm = cfd.timing(cfg, src, mesh.plan(cfg, src))
    bg = cfd.background_state(cfg, src)
    assert tm["nu_max_m2_s"] == pytest.approx(
        src.gas.viscosity(bg["temperature_K"]) / bg["density_kg_m3"], rel=1e-12)


def test_thermo_dict_selects_capped_transport_only_when_enabled():
    cfg = load("hpc-kn1em3-capsule")
    src = config.build_source(cfg)
    nu = cfd.timing(cfg, src, mesh.plan(cfg, src))["nu_max_m2_s"]
    capped = cfd.render_thermophysical_properties(src, nu_max=nu)
    assert "transport       cappedTabulated;" in capped
    assert "nuMax" in capped
    ell = load("hpc-kn1em3-ellipsoid")
    esrc = config.build_source(ell)
    assert cfd.timing(ell, esrc, mesh.plan(ell, esrc))["nu_max_m2_s"] is None
    plain = cfd.render_thermophysical_properties(esrc)
    assert "nuMax" not in plain and "cappedTabulated" not in plain
    assert cfd.render_thermophysical_properties(esrc, nu_max=None) == plain


@pytest.mark.parametrize("factor", [0, -1, float("nan"), float("inf")])
def test_viscosity_cap_factor_must_be_positive(factor):
    with pytest.raises(config.ConfigError):
        load("hpc-kn1em3-capsule", {"cfd": {"viscosity_cap": {"factor": factor}}})
