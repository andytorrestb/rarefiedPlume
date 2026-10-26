"""Read-only recovery of comparison geometry from verified Cartesian mesh inputs.

The archived comparison exports omit C/V, but generated tensor-product meshes
have recoverable centres and volumes. Recovery is accepted only after checking
the actual blockMeshDict and the surviving exported point/count metadata.
No solver, mesher, case generator or physical-data reconstruction is invoked.
"""
from __future__ import annotations

from pathlib import Path
import math
import re
import numpy as np
import yaml

from plumetools.cai2012 import config, geometry, inflow, mesh
from plumetools.cai2012ext import comparison
from .source import DiskSource


class UnsupportedExportGeometry(ValueError):
    """The archived source cannot support a checked geometric reconstruction."""


def _uncomment(text):
    return re.sub(r'/\*.*?\*/|//[^\n]*', '', text, flags=re.S)


def _mesh_signature(text):
    text = _uncomment(text)
    vertices = re.search(r'\bvertices\s*\((.*?)\)\s*;', text, re.S)
    if vertices is None:
        raise UnsupportedExportGeometry('missing blockMeshDict vertices')
    points = [tuple(map(float, item.split())) for item in
              re.findall(r'\(([^()]*)\)', vertices.group(1))]
    blocks = re.findall(r'\bhex\s*\(([^()]*)\)\s*\(([^()]*)\)\s*simpleGrading\s*\(([^()]*)\)', text)
    if not blocks:
        raise UnsupportedExportGeometry('only simpleGrading Cartesian hex blocks are supported')
    connectivity = np.array([list(map(int, block[0].split())) for block in blocks])
    counts = np.array([list(map(int, block[1].split())) for block in blocks])
    grading = np.array([list(map(float, block[2].split())) for block in blocks])
    scale = re.search(r'\b(?:scale|convertToMeters)\s+([^;]+);', text)
    scale = float(scale.group(1)) if scale else 1.
    if counts.ndim != 2 or counts.shape[1] != 3 or grading.shape != counts.shape:
        raise UnsupportedExportGeometry('unsupported mesh grading format')
    return np.array(points)*scale, connectivity, counts, grading


def _verify_dictionary(case, cfg, geom, plan):
    path = case/'system'/'blockMeshDict'
    if not path.is_file():
        raise UnsupportedExportGeometry('system/blockMeshDict is missing')
    observed = _mesh_signature(path.read_text(encoding='utf-8', errors='replace'))
    expected = _mesh_signature(mesh.render_block_mesh_dict(cfg, geom, plan))
    for label, actual, target in zip(('vertices', 'connectivity', 'cell_counts', 'grading'), observed, expected):
        if actual.shape != target.shape:
            raise UnsupportedExportGeometry('configured and archived mesh differ: '+label)
        equal = (np.array_equal(actual, target) if label in ('connectivity', 'cell_counts')
                 else np.allclose(actual, target, rtol=1e-8, atol=1e-12))
        if not equal:
            raise UnsupportedExportGeometry('configured and archived mesh differ: '+label)
    return {'status': 'verified', 'dictionary': str(path), 'method':
            'parsed vertices/connectivity/cell counts/simpleGrading agree with pure mesh planner',
            'blocks': len(observed[2]), 'n_cells': int(np.prod(observed[2], axis=1).sum())}


def load_export_context(case):
    """Return source, axes and verified geometry, without generating a mesh."""
    case = Path(case)
    if (case/'interface.yaml').is_file():
        raise UnsupportedExportGeometry('interface sources are not direct-nozzle geometry')
    cfg = config.load_case_config(case)
    if cfg.mesh.type != 'graded_cartesian':
        raise UnsupportedExportGeometry('only verified graded_cartesian exports are supported')
    geom, state = geometry.from_config(cfg), inflow.from_config(cfg)
    if geom.symmetry_mode != 'none':
        raise UnsupportedExportGeometry('only the full Cartesian plume domain is supported')
    plan = mesh.plan(cfg, geom, state)
    check = _verify_dictionary(case, cfg, geom, plan)
    if check['n_cells'] != plan.n_cells:
        raise UnsupportedExportGeometry('archived total cell count disagrees with mesh planner')
    comparison_path = case/'results'/'comparison.yaml'
    if comparison_path.is_file():
        document = yaml.safe_load(comparison_path.read_text(encoding='utf-8')) or {}
        actual_count = document.get('case', {}).get('n_cells')
        if actual_count is not None and int(actual_count) != plan.n_cells:
            raise UnsupportedExportGeometry('exported result total cell count disagrees with mesh planner')
        check['exported_n_cells'] = actual_count
    segments = (plan.x_segments, plan.y_segments, plan.z_segments)
    edges = tuple(np.r_[np.concatenate([s.cell_edges()[:-1] for s in axis]), axis[-1].hi]
                  for axis in segments)
    centres = tuple(.5*(axis[1:]+axis[:-1]) for axis in edges)
    widths = tuple(np.diff(axis) for axis in edges)
    source = DiskSource(geom.radius_m, state.number_density_per_m3, state.T0_K,
                        state.speed_ratio, state.species)
    context = dict(case=case, cfg=cfg, geom=geom, exit_state=state, source=source,
                   plan=plan, edges_axes=edges, centres_axes=centres,
                   widths_axes=widths, verification=check)
    path = case/'results'/'meridional_plane.csv'
    if path.is_file():
        data = meridional_geometry(case, context)
        check['meridional_points_verified'] = len(data['points_m'])
    return context


def _nearest_indices(axis, values, tolerance):
    indices = np.searchsorted(axis, values)
    upper = np.minimum(indices, len(axis)-1)
    lower = np.maximum(indices-1, 0)
    ids = np.where(np.abs(axis[upper]-values) < np.abs(axis[lower]-values), upper, lower)
    error = np.abs(axis[ids]-values)
    if np.any(error > tolerance):
        raise UnsupportedExportGeometry('exported point is not an archived mesh centre; max error '+str(float(error.max())))
    return ids


def point_volumes(context, points):
    """Recover physical cell volumes only for verified mesh-centre positions.

    CSV coordinates have ten significant digits. Tolerance is 2e-8 D, explicitly
    a serialization tolerance; it is not a physical spatial smoothing width.
    """
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    tolerance = 2e-8*context['geom'].diameter_m
    ids = [_nearest_indices(axis, points[:, j], tolerance)
           for j, axis in enumerate(context['centres_axes'])]
    return np.prod([axis[index] for axis, index in zip(context['widths_axes'], ids)], axis=0)


def meridional_geometry(case, context=None):
    context = load_export_context(case) if context is None else context
    path = Path(case)/'results'/'meridional_plane.csv'
    if not path.is_file():
        raise UnsupportedExportGeometry('meridional_plane.csv is missing')
    table = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
    d = context['geom'].diameter_m
    if 'normal_m' not in table.dtype.names:
        raise UnsupportedExportGeometry('meridional export has no plane-normal coordinate')
    cfg = context['cfg']
    lateral_name = 'z_over_D' if cfg.post.plane == 'xz' else 'y_over_D'
    if lateral_name not in table.dtype.names:
        raise UnsupportedExportGeometry('meridional export has no '+lateral_name)
    points = np.column_stack((table['x_over_D']*d, table['normal_m'], table[lateral_name]*d))
    if cfg.post.plane != 'xz':
        points = points[:, [0, 2, 1]]
    return dict(table=table, points_m=points, volumes_m3=point_volumes(context, points))


def centerline_geometry(context):
    """Recover the exact legacy centre-assigned tube and axial bin weights."""
    cfg, geom = context['cfg'], context['geom']
    x, y, z = context['centres_axes']
    wx, wy, wz = context['widths_axes']
    yy, zz = np.meshgrid(y, z, indexing='ij')
    vv = wy[:, None]*wz[None, :]
    tube = np.hypot(yy, zz) <= float(cfg.post.centerline_radius_over_D)*geom.diameter_m
    x_max = float(cfg.post.centerline_x_over_D_max)*geom.diameter_m
    keep_x = (x >= 0.) & (x <= x_max)
    count = int(cfg.post.centerline_points)
    edges = np.linspace(0., x_max, count+1)
    bin_x = np.clip(np.digitize(x[keep_x], edges)-1, 0, count-1)
    lateral_y, lateral_z, lateral_volume = yy[tube], zz[tube], vv[tube]
    points = np.column_stack((np.repeat(x[keep_x], len(lateral_y)),
                              np.tile(lateral_y, keep_x.sum()), np.tile(lateral_z, keep_x.sum())))
    volume = (wx[keep_x, None]*lateral_volume).ravel()
    index = np.repeat(bin_x, len(lateral_y))
    n_cells = np.bincount(index, minlength=count)
    volume_sum = np.bincount(index, weights=volume, minlength=count)
    result = dict(points_m=points, volumes_m3=volume, bin_index=index, n_bins=count,
                  bin_edges_m=edges, bin_centres_over_D=.5*(edges[:-1]+edges[1:])/geom.diameter_m,
                  n_cells=n_cells, volume_sum_m3=volume_sum, valid=n_cells > 0)
    path = context['case']/'results'/'centerline_extended.csv'
    if path.is_file():
        table = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        ids = _nearest_indices(result['bin_centres_over_D'], table['x_over_D'], 2e-8)
        if not np.allclose(table['n_cells'], n_cells[ids], rtol=1e-8, atol=1e-8):
            raise UnsupportedExportGeometry('exported centreline counts disagree with reconstructed tube')
        result.update(table=table, export_bin_index=ids)
        context['verification']['centerline_export_bins_verified'] = len(ids)
    return result


def angular_geometry(context, *, radius_over_D=comparison.ANGULAR_RADIUS_OVER_D,
                     half_width_over_D=comparison.SHELL_HALF_WIDTH_OVER_D,
                     bin_deg=comparison.ANGULAR_BIN_DEG):
    """Recover shell cells and angular footprints with the exact legacy rule.

    Build one x-slab at a time, avoiding a full 3-D coordinate/volume allocation.
    Returned fractions have cells x angular_bins shape and match footprint_fractions.
    """
    geom = context['geom']
    x, y, z = context['centres_axes']
    wx, wy, wz = context['widths_axes']
    yy, zz = np.meshgrid(y, z, indexing='ij')
    lateral_volume = wy[:, None]*wz[None, :]
    target, half = radius_over_D*geom.diameter_m, half_width_over_D*geom.diameter_m
    points, volumes = [], []
    for xx, width in zip(x, wx):
        volume = width*lateral_volume
        big_r = np.sqrt(xx*xx+yy*yy+zz*zz)
        keep = (np.abs(big_r-target) <= half+.5*np.cbrt(volume)) & (xx > 0.)
        if np.any(keep):
            points.append(np.column_stack((np.full(int(keep.sum()), xx), yy[keep], zz[keep])))
            volumes.append(volume[keep])
    points = np.concatenate(points) if points else np.empty((0, 3))
    volume = np.concatenate(volumes) if volumes else np.empty(0)
    radius = np.linalg.norm(points, axis=1)
    theta = np.arctan2(np.hypot(points[:, 1], points[:, 2]), points[:, 0])
    edges = np.arange(0., 90.+.5*bin_deg, bin_deg)
    edges[-1] = min(edges[-1], 90.)
    fractions = comparison.footprint_fractions(theta, radius, volume, np.radians(edges))
    result = dict(points_m=points, volumes_m3=volume, fractions=fractions,
                  theta_rad=theta, radius_m=radius, theta_edges_deg=edges,
                  theta_deg=.5*(edges[:-1]+edges[1:]), n_cells=fractions.sum(axis=0),
                  solid_angle_sr=2*math.pi*(np.cos(np.radians(edges[:-1]))-np.cos(np.radians(edges[1:]))))
    path = context['case']/'results'/'angular_r10D.csv'
    if path.is_file() and radius_over_D == 10. and half_width_over_D == .5 and bin_deg == 2.5:
        table = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))
        if not (np.allclose(table['theta_deg'], result['theta_deg'], atol=1e-8, rtol=0)
                and np.allclose(table['n_cells'], result['n_cells'], atol=2e-6, rtol=2e-8)):
            raise UnsupportedExportGeometry('exported angular footprints disagree with reconstructed shell')
        result['table'] = table
        context['verification']['angular_export_bins_verified'] = len(table)
    return result
