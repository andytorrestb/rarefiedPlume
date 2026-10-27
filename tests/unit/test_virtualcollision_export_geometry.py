"""Independent geometry/count checks for read-only compact-export recovery."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from plumetools.cai2012 import config, geometry, inflow, mesh
from plumetools.virtualcollision.export_geometry import (
    UnsupportedExportGeometry, load_export_context, point_volumes,
    centerline_geometry, angular_geometry,
)


def _context(tmp_path):
    return dict(case=tmp_path, cfg=SimpleNamespace(post=SimpleNamespace(
        centerline_radius_over_D=.1, centerline_x_over_D_max=.4,
        centerline_points=2)), geom=SimpleNamespace(diameter_m=1.),
        centres_axes=(np.array([.1, .3]), np.array([-.05, .05]), np.array([-.05, .05])),
        widths_axes=(np.array([.2, .2]), np.array([.1, .1]), np.array([.1, .1])),
        verification={})


def test_cell_volumes_recover_nonuniform_cartesian_products_and_reject_noncentres(tmp_path):
    context = _context(tmp_path)
    context['widths_axes'] = (np.array([.2, .4]), np.array([.1, .3]), np.array([.1, .5]))
    points = np.array([[.1, -.05, -.05], [.3, .05, .05]])
    np.testing.assert_allclose(point_volumes(context, points), [.002, .06], rtol=1e-14)
    with pytest.raises(UnsupportedExportGeometry, match='mesh centre'):
        point_volumes(context, [[.11, -.05, -.05]])


def test_tube_geometry_recovers_known_cell_counts_and_physical_volume(tmp_path):
    result = centerline_geometry(_context(tmp_path))
    np.testing.assert_array_equal(result['n_cells'], [4, 4])
    np.testing.assert_allclose(result['volume_sum_m3'], [.008, .008])
    np.testing.assert_allclose(result['bin_centres_over_D'], [.1, .3])
    assert len(result['points_m']) == 8
    assert np.all(np.hypot(result['points_m'][:, 1], result['points_m'][:, 2]) < .1)


def test_shell_geometry_uses_volume_reach_and_known_angular_footprints(tmp_path):
    context = _context(tmp_path)
    context.update(centres_axes=(np.array([.2, .8]), np.array([-.2, .2]), np.array([-.2, .2])),
                   widths_axes=(np.array([.5, .5]), np.array([.4, .4]), np.array([.4, .4])))
    result = angular_geometry(context, radius_over_D=1., half_width_over_D=.3, bin_deg=45.)
    assert len(result['points_m']) == 4
    np.testing.assert_allclose(result['points_m'][:, 0], .8)
    np.testing.assert_allclose(result['volumes_m3'], .08)
    np.testing.assert_allclose(result['n_cells'], [4., 0.])
    assert result['solid_angle_sr'].sum() == pytest.approx(2*np.pi)


def _archived_case(tmp_path):
    root = Path(__file__).resolve().parents[2]
    document = yaml.safe_load((root/'cases/cai2012-kn-extension/baseCase/case.yaml').read_text())
    document['mesh'].update(max_core_cell_over_D=.2, core_x_over_D=.5,
                            core_half_over_D=.6, outer_x_cells=4, outer_lateral_cells=4,
                            max_cells=200000)
    case = tmp_path/'case'
    (case/'system').mkdir(parents=True)
    (case/'results').mkdir()
    (case/'case.yaml').write_text(yaml.safe_dump(document))
    cfg = config.load_case_config(case)
    geom, state = geometry.from_config(cfg), inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, state)
    text = mesh.render_block_mesh_dict(cfg, geom, plan)
    (case/'system/blockMeshDict').write_text(text)
    (case/'results/comparison.yaml').write_text(yaml.safe_dump({'case': {'n_cells': plan.n_cells}}))
    return case, text


def test_archived_dictionary_is_verified_and_changed_grading_is_refused(tmp_path):
    case, text = _archived_case(tmp_path)
    recovered = load_export_context(case)
    assert recovered['verification']['status'] == 'verified'
    assert np.prod([len(a) for a in recovered['centres_axes']]) == recovered['plan'].n_cells
    # Changing a generated physical mesh input must invalidate recovery.
    (case/'system/blockMeshDict').write_text(text.replace('scale   1;', 'scale   2;'))
    with pytest.raises(UnsupportedExportGeometry, match='mesh differ: vertices'):
        load_export_context(case)


def test_exported_total_count_mismatch_is_refused(tmp_path):
    case, _ = _archived_case(tmp_path)
    (case/'results/comparison.yaml').write_text('case:\n  n_cells: 99\n')
    with pytest.raises(UnsupportedExportGeometry, match='total cell count'):
        load_export_context(case)
