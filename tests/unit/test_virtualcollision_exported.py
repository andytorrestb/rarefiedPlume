"""Export adapter limits and preservation of legacy/raw temperature semantics."""
import numpy as np
import pytest
from plumetools.cai2012 import analytical
from plumetools.virtualcollision import DiskSource, NumericalOptions, VirtualCollisionModel
from plumetools.virtualcollision.kernel import ZeroKernel, VHSKernel
from plumetools.virtualcollision.exported import BallisticCache, PointPrediction, correction_grid, _profile
from plumetools.virtualcollision.moments import derive


def options():
    return NumericalOptions(samples=8, batches=2, hazard_steps=2,
        partner_radial=2, partner_azimuthal=4, partner_speed=8)


def test_zero_signed_correction_retains_exact_cai_points_after_interpolation():
    source = DiskSource(.1, 1e17, 300., 2.)
    model = VirtualCollisionModel(source, ZeroKernel(), options())
    points = np.array([[.3, .08, 0], [.8, .15, 0]])
    prediction = PointPrediction(model, correction_grid(.2, refinement=0), BallisticCache()).evaluate(points)
    expected = analytical.moments(*points.T, .1, 2.)
    np.testing.assert_allclose(prediction['fields']['n']/source.number_density_m3,
                               expected['n_over_n0'], rtol=1e-14)
    np.testing.assert_allclose(prediction['fields']['U']/source.thermal_speed,
                               expected['U_sqrt_beta0'], rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(prediction['fields']['T_trans']/300., expected['T_over_T0'], rtol=1e-13)


def test_local_temperature_average_and_raw_aggregate_are_kept_distinct():
    source = DiskSource(.1, 1., 300., 2.)
    velocity = np.array([[100., 0, 0], [-100., 0, 0]])
    fields = derive(np.ones(2), velocity, velocity[:, :, None]*velocity[:, None, :], source.species.mass_kg)
    prediction = dict(fields=fields, ballistic=fields, batch_fields=[fields, fields])
    geometry = dict(volumes_m3=np.ones(2), bin_index=np.zeros(2, dtype=int),
                    n_bins=1, export_bin_index=np.array([0]))
    columns, _ = _profile(prediction, geometry, source)
    assert columns['T_over_T0_virtual'][0] == 0
    assert columns['T_over_T0_raw_aggregated_virtual'][0] == pytest.approx(
        10000/(3*source.species.gas_constant_j_per_kg_k*300))


def test_raw_correction_api_agrees_with_full_distribution_increment():
    source = DiskSource(.1, 1e17, 300., 2.)
    model = VirtualCollisionModel(source, VHSKernel(), options())
    xb, rb = np.array([[.1, .3], [.5, .8]]), np.array([[0, .1], [.1, .3]])
    delta = model.raw_correction_bins(xb, rb)
    full = model.evaluate_bins(xb, rb)
    for key in ('M0', 'M1', 'M2'):
        np.testing.assert_allclose(delta[key], full[key]-full['ballistic_raw_moments'][key],
                                   rtol=1e-10, atol=source.number_density_m3*1e-12)
