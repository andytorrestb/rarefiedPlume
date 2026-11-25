"""Collisionless plate loads from Cai's exit Maxwellian (spec §9)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from plumetools.cai2012 import analytical
from plumetools.plateugf import caiplate

K = 1.380649e-23
M = 6.63e-26
T0 = 300.0
S0 = 2.0
R = 0.1
N0 = 6.62e16
VP = math.sqrt(2 * K * T0 / M)           # most probable speed
U0 = S0 * VP


def _moments(points):
    return caiplate.moments(np.asarray(points, float), radius_m=R, number_density=N0,
                            speed=U0, temperature_K=T0, mass_kg=M)


def test_density_matches_cai_free_plume():
    pts = np.array([[1.0, 0.3, 0.0], [5.0, 2.0, 1.0], [2.0, 0.0, 0.0],
                    [10.0, 4.0, 3.0], [0.5, 0.2, 0.1], [20.0, 5.0, 0.0]])
    n = _moments(pts)["n"] / N0
    expected = analytical.density_ratio(pts[:, 0], pts[:, 1], pts[:, 2], R, S0)
    assert n == pytest.approx(expected, rel=2e-3)


def test_axial_speed_matches_cai_centreline():
    x = np.array([2.0, 5.0, 10.0])
    m = _moments(np.column_stack([x, np.zeros(3), np.zeros(3)]))
    speed_ratio = m["flux"] / m["n"] / VP
    assert speed_ratio == pytest.approx(analytical.centerline_axial_speed_ratio(x, R, S0), rel=2e-3)


def test_far_field_pressure_falls_as_inverse_square():
    loads = caiplate.plate_loads(np.array([[20.0, 0, 0], [40.0, 0, 0]]), radius_m=R,
                                 number_density=N0, speed=U0, temperature_K=T0,
                                 mass_kg=M, wall_temperature_K=300.0)
    assert loads.pressure[1] / loads.pressure[0] == pytest.approx(0.25, rel=0.03)
    assert np.all(loads.pressure > 0) and np.all(loads.heat_flux > 0)


def test_shear_is_zero_on_axis_and_points_outward():
    loads = caiplate.plate_loads(np.array([[5.0, 0, 0], [5.0, 1.0, 0], [5.0, 0, 2.0]]),
                                 radius_m=R, number_density=N0, speed=U0,
                                 temperature_K=T0, mass_kg=M, wall_temperature_K=300.0)
    assert abs(loads.shear[0]).max() < 1e-6 * loads.pressure[0]
    assert loads.shear[1, 0] > 0 and abs(loads.shear[1, 1]) < 1e-9 * loads.shear[1, 0]
    assert loads.shear[2, 1] > 0 and abs(loads.shear[2, 0]) < 1e-9 * loads.shear[2, 1]


def test_diffuse_reemission_adds_to_the_incident_pressure():
    loads = caiplate.plate_loads(np.array([[5.0, 0, 0]]), radius_m=R, number_density=N0,
                                 speed=U0, temperature_K=T0, mass_kg=M, wall_temperature_K=300.0)
    m = _moments([[5.0, 0, 0]])
    reemitted = M * m["flux"][0] * math.sqrt(math.pi * K * 300.0 / (2 * M))
    assert loads.pressure[0] == pytest.approx(m["pxx"][0] + reemitted, rel=1e-12)
    assert loads.heat_flux[0] == pytest.approx(m["energy"][0] - 2 * K * 300.0 * m["flux"][0], rel=1e-12)
