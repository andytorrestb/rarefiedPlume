"""Fitted-capsule mode routing and the coarse 30D x 15D profile; no solver."""
from pathlib import Path

import pytest

from plumetools.hybridplume import config
from plumetools.hybridplume import mesh

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases" / "cai2012-hybrid" / "hybrid.yaml"


def load(name, overlay=None):
    return config.load(STUDY, name, overlay)


def test_capsule_profile_mesh_matches_spec():
    cfg = load("hpc-kn1em3-capsule")
    src = config.build_source(cfg)
    p = mesh.plan(cfg, src)
    assert 1.90e6 < p.n_cells < 2.02e6
    assert p.max_expansion_ratio < 1.11
    assert abs(mesh.predicted_source(p, src).area_error) < 0.01
    assert cfg.get("interface.mode") == "fitted-capsule"


def test_is_spatial_fit():
    assert config.is_spatial_fit(load("hpc-kn1em3-capsule"))
    assert config.is_spatial_fit(load("hpc-kn1em3-ellipsoid"))
    assert not config.is_spatial_fit(load("planar-laptop-kn0p01"))
    assert config.is_fitted(load("planar-laptop-kn0p01"))
    assert not config.is_fitted(load("smoke-kn0p01"))


def test_fitted_capsule_requires_cartesian():
    with pytest.raises(config.ConfigError):
        load("hpc-kn1em3-capsule", {"mesh": {"geometry": "planar"}})


def test_allrun_has_no_literal_ellipsoid_mode_tests():
    text = (ROOT / "cases/cai2012-hybrid/Allrun").read_text()
    assert "= fitted-ellipsoid ]" not in text and "spatial_fit()" in text
