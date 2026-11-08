"""Acceptance uses measured local MFP; already adequate cells need no sweep."""
from __future__ import annotations

from dataclasses import replace
import importlib

import numpy as np
import pytest

from _refinement import unit_inputs


def assess(*args, **kwargs):
    return importlib.import_module("plumetools.refinement.assess").assess(*args, **kwargs)


def test_measured_edge_target_and_twenty_parcels_are_inclusive():
    g, m, gas = unit_inputs(edges=[0.8, 0.81, 0.4], lambdas=[1, 1, 1], parcels=[20, 40, 20])
    result = assess(g, m, gas)
    assert result.mesh_status.tolist() == [1, 0, 1]
    assert result.occupancy_status.tolist() == [1, 1, 1]
    assert result.lambda_m == pytest.approx([1, 1, 1])
    assert result.target_edge_m == pytest.approx([0.8, 0.8, 0.8])
    assert result.summary["action"] == "needs_mesh"
    assert result.summary["mesh_failed_cells"] == 1


def test_passing_cells_request_no_additional_refinement():
    g, m, gas = unit_inputs(edges=[0.79, 0.3], lambdas=[1, 1], parcels=[20, 80])
    result = assess(g, m, gas)
    assert result.summary["action"] == "criteria_satisfied"
    assert result.mesh_status.tolist() == [1, 1]
    assert result.target_edge_m == pytest.approx([0.8, 0.8])


@pytest.mark.parametrize("role,value", [
    ("number_density", 0), ("number_density", np.nan),
    ("temperature_K", -1), ("temperature_K", np.inf),
    ("mean_parcels", 0), ("mean_parcels", np.nan),
])
def test_unreliable_measurements_remain_unresolved(role, value):
    g, m, gas = unit_inputs(edges=[0.2, 0.2], lambdas=[1, 1], parcels=[40, 40])
    values = getattr(m, role).copy()
    values[1] = value
    result = assess(g, replace(m, **{role: values}), gas)
    assert result.mesh_status[1] == -1
    assert not result.valid[1]
    assert result.summary["action"] == "unresolved"
    assert result.summary["unresolved_volume_fraction"] == pytest.approx(0.5)


def test_missing_window_cannot_become_a_passing_result():
    g, m, gas = unit_inputs(edges=[0.2], lambdas=[1], parcels=[40])
    result = assess(g, replace(m, window={}), gas)
    assert result.summary["action"] == "unresolved"
    assert result.mesh_status.tolist() == [-1]


def test_low_occupancy_requests_population_before_mesh_correction():
    g, m, gas = unit_inputs(edges=[2.0, 0.2], lambdas=[1, 1], parcels=[10, 40])
    result = assess(g, m, gas)
    assert result.summary["action"] == "needs_parcels"
    assert result.mesh_status.tolist() == [-1, 1]
    assert result.occupancy_status.tolist() == [0, 1]


def test_placeholder_empty_cell_does_not_supply_trusted_mfp():
    g, m, gas = unit_inputs(edges=[0.2], lambdas=[1], parcels=[0.001])
    result = assess(g, replace(m, number_density=np.array([0.0])), gas)
    assert result.mesh_status.tolist() == [-1]
    assert result.summary["action"] == "unresolved"


def test_custom_targets_are_used_without_changing_the_data():
    from plumetools.refinement.model import Targets
    g, m, gas = unit_inputs(edges=[0.5], lambdas=[1], parcels=[40])
    result = assess(g, m, gas, targets=Targets(edge_over_mfp=0.4, min_mean_parcels=40))
    assert result.mesh_status.tolist() == [0]
    assert result.occupancy_status.tolist() == [1]
    assert m.mean_parcels.tolist() == [40]


def test_mismatched_geometry_and_measurement_counts_are_refused():
    g, m, gas = unit_inputs(edges=[0.5], lambdas=[1], parcels=[40])
    with pytest.raises(ValueError, match="cell"):
        assess(g, replace(m, number_density=np.array([1.0, 2.0])), gas)


def test_zero_cell_volume_is_refused():
    from plumetools.refinement.model import Geometry
    with pytest.raises(ValueError, match="volume"):
        Geometry(np.zeros((1, 3)), np.zeros(1), np.ones(1), "mesh")
