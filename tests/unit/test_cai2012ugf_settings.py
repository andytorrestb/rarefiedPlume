"""unigas.yaml: the uniGasFoam settings beside each cai2012-unigas case.yaml.

The defaults are the values the approved design fixes
(docs/superpowers/specs/2026-11-19-cai-kn-sweep-unigasfoam-design.md, §4), and
the owner's decisions 4 and 5 are refusals, not options.
"""

from __future__ import annotations

import pytest

from plumetools.cai2012ugf.settings import (
    UniGasConfigError,
    load_settings,
    settings_from_dict,
)


def test_defaults_are_the_spec_values():
    s = settings_from_dict({}, "test")
    assert s.cell_weighted is False
    assert (s.collision.model, s.collision.bgk_model) == ("hybrid", "unifiedStochasticParticleSBGK")
    assert (s.collision.dsmc_partner_model, s.collision.dsmc_model) == ("noTimeCounter", "variableHardSphere")
    assert s.collision.macro_interpolation is True and s.collision.theta == 1e-3
    a = s.adaptive
    assert (a.adaptation_interval, a.max_sub_cell_size_mfp_ratio, a.max_time_step_mct_ratio,
            a.max_courant_number, a.dsmc_cells_time_step_limited, a.max_sub_cell_levels) == (20, 0.5, 0.2, 0.5, True, 10)
    d = s.decomposition
    assert (d.breakdown_max, d.theta, d.smoothing_passes, d.decomposition_interval) == (0.05, 0.2, 5, 100)
    assert (s.writes.startup_interval_transits, s.writes.keep_lagrangian) == (1.0, False)
    assert (s.fields.mean_field, s.fields.window_field) == ("Ar", "ArWindow")
    v = s.verdict
    assert v.max_effective_spacing_over_mfp is None and v.min_particles_per_dsmc_cell is None
    assert (v.max_time_step_mct_ratio, v.min_particles_per_sub_cell, v.occupancy_floor_n_over_n0) == (0.2, 2.0, 1e-3)


def test_cell_weighting_is_refused():           # decision 5
    with pytest.raises(UniGasConfigError, match="cell_weighted"):
        settings_from_dict({"cell_weighted": True}, "x")


def test_unknown_key_is_named():
    with pytest.raises(UniGasConfigError, match="adaptive.maxCourant"):
        settings_from_dict({"adaptive": {"maxCourant": 1}}, "x")


def test_non_usp_or_non_hybrid_is_refused():
    for bad in ({"collision": {"model": "dsmc"}}, {"collision": {"bgk_model": "stochasticParticleSBGK"}}):
        with pytest.raises(UniGasConfigError):
            settings_from_dict(bad, "x")


def test_time_step_limit_off_is_refused():      # decision 4
    with pytest.raises(UniGasConfigError, match="dsmc_cells_time_step_limited"):
        settings_from_dict({"adaptive": {"dsmc_cells_time_step_limited": False}}, "x")


def test_same_mean_and_window_field_is_refused():
    with pytest.raises(UniGasConfigError):
        settings_from_dict({"fields": {"window_field": "Ar"}}, "x")


def test_missing_file_names_its_path(tmp_path):
    with pytest.raises(UniGasConfigError, match="unigas.yaml"):
        load_settings(tmp_path)
