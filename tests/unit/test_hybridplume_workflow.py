"""hybridplume: configuration and its refusals."""

from __future__ import annotations

from pathlib import Path

import pytest

from plumetools.hybridplume import config

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "cases" / "cai2012-hybrid" / "hybrid.yaml"


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
