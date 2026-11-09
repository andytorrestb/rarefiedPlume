"""Small numerical and setup checks for the opt-in planar hybrid prototype."""
from pathlib import Path

import numpy as np
import pytest


def ellipse_grid():
    x, y = np.meshgrid(np.arange(20) / 10, np.arange(-10, 10) / 10, indexing="ij")
    lower = np.column_stack((x.ravel(), y.ravel()))
    upper = lower + 0.1
    corners = np.maximum(np.abs(lower), np.abs(upper))
    retained = (corners[:, 0] / 1.4)**2 + (corners[:, 1] / 0.45)**2 < 1
    return lower, upper, retained


def test_fit_respects_resolved_ellipse():
    from plumetools.hybridplume.ellipse import fit_ellipse, ellipse_points
    lower, upper, retained = ellipse_grid()
    fit = fit_ellipse(lower, upper, retained, source_radius=0.5, clearance=0.01)
    assert fit["status"] == "fitted"
    assert fit["area_m2"] > 0.45
    arc = ellipse_points(fit, 301)
    assert np.all((arc[:, 0] / 1.4)**2 + (arc[:, 1] / 0.45)**2 <= 1 + 1e-7)
    assert np.allclose(arc[[0, -1], 0], 0, atol=1e-9)
    assert np.max(np.abs(arc[[0, -1], 1])) <= 0.5


def test_excluded_interior_pocket_limits_the_whole_ellipse():
    from plumetools.hybridplume.ellipse import fit_ellipse
    lower, upper, retained = ellipse_grid()
    hole = (np.isclose(lower[:, 0], 0.6)
            & (lower[:, 1] >= -0.1 - 1e-8) & (upper[:, 1] <= 0.1 + 1e-8))
    retained[hole] = False
    fit = fit_ellipse(lower, upper, retained, source_radius=0.5, clearance=0.01)
    assert fit["status"] == "fitted"
    assert fit["xc"] + fit["a"] <= 0.59 + 1e-6


def test_no_source_connected_ellipse_is_reported():
    from plumetools.hybridplume.ellipse import fit_ellipse
    lower, upper, retained = ellipse_grid()
    retained[lower[:, 0] < 0.3] = False
    fit = fit_ellipse(lower, upper, retained, source_radius=0.5, clearance=0.01)
    assert fit["status"] == "no_feasible_ellipse"
    assert "a" not in fit
