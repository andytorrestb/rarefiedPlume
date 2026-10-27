"""Search acceleration preserves the physical residence estimator."""
from dataclasses import replace
import numpy as np

from plumetools.virtualcollision import DiskSource, NumericalOptions, VirtualCollisionModel
from plumetools.virtualcollision.kernel import VHSKernel, ZeroKernel


def model(**changes):
    options = NumericalOptions(samples=8, batches=2, hazard_steps=4,
        partner_radial=2, partner_azimuthal=4, partner_speed=8,
        domain_length_D=4, domain_radius_D=4, cai_radial=8,
        cai_azimuthal=16, **changes)
    source = DiskSource(.1, 1e17, 300., 2.)
    return VirtualCollisionModel(source, VHSKernel(source.species), options)


def test_ray_search_envelope_preserves_all_residence_moments():
    transport = model()
    xedges, redges = np.linspace(0, .8, 9), np.linspace(0, .8, 9)
    xb = np.repeat(np.column_stack((xedges[:-1], xedges[1:])), 8, axis=0)
    rb = np.tile(np.column_stack((redges[:-1], redges[1:])), (8, 1))
    velocities = np.array([[1., 0., 0.], [-1., 0., 0.], [0., 1., 0.],
                           [.3, -.8, .2], [.3, .1, -.9], [-.7, -.3, .2]])
    position = np.array([.35, .27, .05])
    for velocity in velocities:
        unfiltered = transport._ray_moments(position, velocity, 3., xb, rb)
        selected = transport._candidate_bins(position, velocity, xb, rb)
        accelerated = [np.zeros_like(a) for a in unfiltered]
        subset = transport._ray_moments(position, velocity, 3., xb[selected], rb[selected])
        for destination, values, expected in zip(accelerated, subset, unfiltered):
            destination[selected] = values
            np.testing.assert_array_equal(destination, expected)
        assert np.all(np.isin(np.flatnonzero(unfiltered[0]), selected))


def test_chunking_does_not_change_transport_raw_moments():
    first = model(chunk_size=3)
    second = VirtualCollisionModel(first.source, first.kernel, replace(first.options, chunk_size=17))
    xedges, redges = np.linspace(.2, .8, 5), np.linspace(0, .5, 5)
    xb = np.repeat(np.column_stack((xedges[:-1], xedges[1:])), 4, axis=0)
    rb = np.tile(np.column_stack((redges[:-1], redges[1:])), (4, 1))
    a, b = first.evaluate_bins(xb, rb), second.evaluate_bins(xb, rb)
    for key in ('M0', 'M1', 'M2'):
        np.testing.assert_allclose(a[key], b[key], rtol=2e-14, atol=1e-14*np.max(np.abs(a[key])))
    np.testing.assert_array_equal(a['valid'], b['valid'])
