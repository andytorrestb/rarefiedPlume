"""Plot exported comparisons without mutating or repairing scientific arrays."""
from pathlib import Path
import json

import numpy as np
import pytest

from plumetools.virtualcollision.export_figures import plot_case


def test_three_method_plots_keep_invalid_raw_data_and_missing_quantities(tmp_path):
    pytest.importorskip('matplotlib')
    x, z = np.meshgrid(np.linspace(.1, 2., 4), np.linspace(-1., 1., 4))
    n = .2*np.exp(-z.ravel()**2)/(1+x.ravel())
    virtual = n*1.01
    virtual[3] = -0.01
    plane = dict(x_over_D=x.ravel(), z_over_D=z.ravel(), normal_m=np.full(x.size, .005),
                 n_over_n0=n, n_over_n0_analytical=n*.98, n_over_n0_virtual=virtual,
                 virtual_valid=np.arange(x.size) != 3)
    profile = dict(x_over_D=np.linspace(.1, 2., 4),
        n_over_n0=np.array([.8, .4, .2, .1]),
        n_over_n0_analytical=np.array([.82, .42, .21, .11]),
        n_over_n0_virtual=np.array([.81, .41, -.001, .105]),
        n_over_n0_virtual_se=np.full(4, .001), virtual_valid=np.array([1, 1, 0, 1]))
    angular = dict(theta_deg=np.array([1.25, 3.75, 6.25, 8.75]),
        n_over_n0=profile['n_over_n0'], n_over_n0_analytical=profile['n_over_n0_analytical'],
        n_over_n0_virtual=profile['n_over_n0_virtual'], virtual_valid=profile['virtual_valid'])
    raw = virtual.copy()
    profile_raw = profile['n_over_n0_virtual'].copy()
    paths = plot_case(tmp_path, plane, profile, angular,
                      dict(case_name='synthetic test fixture', kn=10, quality='exploratory',
                           estimator_note='Synthetic fixture only; not DSMC validation evidence.'))
    assert {Path(path).name for path in paths} == {
        'density_contour_comparison.png', 'density_contour_comparison.svg',
        'density_three_way.png', 'density_errors.png',
        'centerline_density_three_way.png', 'angular_density_three_way.png'}
    assert all(Path(path).is_file() and Path(path).stat().st_size > 1000 for path in paths)
    np.testing.assert_array_equal(plane['n_over_n0_virtual'], raw)
    np.testing.assert_array_equal(profile['n_over_n0_virtual'], profile_raw)
    contours = json.loads((tmp_path/'density_contour_comparison.json').read_text())
    assert contours['levels'] == [.001, .01, .1]
    assert contours['common_masked_nodes'] == 1
    assert contours['methods']['Cai collisionless, same cells']['level_status']['0.001'] == 'outside-displayed-density-range'
