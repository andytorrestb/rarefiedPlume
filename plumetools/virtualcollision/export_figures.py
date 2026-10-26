"""Three-method figures from preserved Cai/DSMC CSVs plus model predictions.

These routines only read arrays. Display masks needed by logarithmic figures
are separate from raw predictions, and their excluded counts are printed on
each figure. Reference quality and finite-grid limitations remain visible.
"""
from __future__ import annotations

from pathlib import Path
import textwrap

import numpy as np


VIRTUAL_LABEL = 'Virtual-collision model'
DENSITY_CONTOUR_LEVELS = (0.001, 0.01, 0.1)


def _density_contours_on_axes(ax, data, metadata):
    """Overlay the legacy density levels using one shared display mask."""
    import matplotlib.tri as mtri
    from matplotlib.lines import Line2D
    ids = _plane_layer(data)
    required = ('n_over_n0', 'n_over_n0_analytical', 'n_over_n0_virtual')
    if len(ids) < 3 or any(key not in data for key in required):
        ax.text(.5, .5, 'Density contour data unavailable', transform=ax.transAxes, ha='center')
        return dict(status='missing-density-data', levels=list(DENSITY_CONTOUR_LEVELS)), []
    x, z = np.asarray(data['x_over_D'])[ids], np.asarray(data['z_over_D'])[ids]
    arrays = [np.asarray(data[key], dtype=float)[ids] for key in required]
    tri = mtri.Triangulation(x, z)
    bad_by_method = [~np.isfinite(values) | (values <= 0) for values in arrays]
    common_bad = np.logical_or.reduce(bad_by_method)
    tri.set_mask(np.any(common_bad[tri.triangles], axis=1))
    displayed_nodes = np.zeros(len(ids), dtype=bool)
    displayed_nodes[np.unique(tri.triangles[~tri.mask])] = True
    styles = [('tab:red', '-', 1.4), ('black', '--', 1.05), ('tab:orange', '-.', 1.6)]
    labels = _labels(metadata, 'plane')
    records, handles = {}, []
    for index, (values, label, (color, style, width)) in enumerate(zip(arrays, labels, styles)):
        handles.append(Line2D([], [], color=color, linestyle=style, lw=width, label=label))
        supported = values[displayed_nodes]
        levels = [level for level in DENSITY_CONTOUR_LEVELS
                  if supported.size and supported.min() < level < supported.max()]
        statuses = {str(level): 'outside-displayed-density-range' for level in DENSITY_CONTOUR_LEVELS}
        if levels:
            # Replacements occur only on masked display nodes, never on a
            # scientific array or a displayed triangle.
            display = np.where(common_bad, 1., values)
            contours = ax.tricontour(tri, display, levels=levels, colors=color,
                                     linestyles=style, linewidths=width)
            for level, segments in zip(contours.levels, contours.allsegs):
                statuses[str(float(level))] = ('present' if any(len(segment) >= 2 for segment in segments)
                                              else 'no-contour-crossing')
            # Label the smooth Cai curves once; labels on every noisy model
            # fragment would obscure the shared-level comparison.
            if index == 1:
                ax.clabel(contours, fmt=lambda value: f'{value:g}', fontsize=8, inline=True)
        records[label] = dict(level_status=statuses,
                              nonpositive_or_nonfinite=int(bad_by_method[index].sum()))
    ax.set(xlabel='$X/D$', ylabel='$z/D$')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(alpha=.15)
    return (dict(status='available' if displayed_nodes.any() else 'no-common-displayable-triangles',
                levels=list(DENSITY_CONTOUR_LEVELS), methods=records,
                common_masked_nodes=int(common_bad.sum()), selected_nodes=len(ids),
                virtual_moment_invalid=int((~_valid(data, len(data['n_over_n0']))[ids]).sum())),
            handles)


def plot_density_contours(output, data, metadata):
    """Add an overlaid three-method density contour comparison, without transport."""
    from .assessment import write_json
    plt = _pyplot()
    if plt is None:
        return []
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.7, 7.2))
    status, handles = _density_contours_on_axes(ax, data, metadata)
    if handles:
        ax.legend(handles=handles, fontsize=8, loc='upper right')
    note = ('Same plane layer and common nonpositive/nonfinite display mask; '
            f'{status.get("common_masked_nodes", 0)} masked nodes. '
            f'Virtual moment-invalid: {status.get("virtual_moment_invalid", 0)}/'
            f'{status.get("selected_nodes", 0)}. '
            'Labels mark Cai density levels; raw arrays are unchanged. Missing crossings are recorded in JSON.')
    _finish(fig, metadata, r'Density contours: $n/n_0=0.1,\;0.01,\;0.001$', note)
    paths = [output/'density_contour_comparison.png', output/'density_contour_comparison.svg']
    for path in paths:
        fig.savefig(path, dpi=180)
    plt.close(fig)
    write_json(output/'density_contour_comparison.json', status)
    return paths


def _pyplot():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt


def _metadata_title(metadata):
    name = metadata.get('case_name', 'Cai plume')
    kn = metadata.get('kn', metadata.get('nominal_kn', 'unknown'))
    quality = str(metadata.get('quality', 'exploratory')).upper()
    return f'{name} | Kn={kn} | {quality}'


def _finish(fig, metadata, subtitle, note='', *, colorbar=False):
    fig.suptitle(_metadata_title(metadata)+'\n'+subtitle, fontsize=10)
    estimator = metadata.get('estimator_note',
        'Virtual prediction uses a finite-grid residence estimator; spatial convergence is not established.')
    footer = str(estimator)
    if note:
        footer += '\n'+note
    fig.text(.5, .012, '\n'.join(textwrap.wrap(footer, width=150)),
             ha='center', va='bottom', fontsize=7)
    if colorbar:
        fig.subplots_adjust(left=.065, right=.88, bottom=.23, top=.77, wspace=.18)
    else:
        fig.tight_layout(rect=(0, .09, 1, .86))


def _valid(data, size):
    value = data.get('virtual_valid')
    if value is None:
        return np.ones(size, dtype=bool)
    value = np.asarray(value)
    if value.shape != (size,):
        return np.ones(size, dtype=bool)
    return np.isfinite(value) & value.astype(bool)


def _labels(metadata, kind):
    name = metadata.get('case_name', 'DSMC')
    baseline = ('Cai collisionless, same shell cells' if kind == 'angular'
                else 'Cai collisionless, same cells')
    return (f'dsmcFoam, {name}', baseline, VIRTUAL_LABEL)


def _plane_layer(data):
    """Select the same positive plane layer used by legacy plane figures."""
    x, z = np.asarray(data['x_over_D']), np.asarray(data['z_over_D'])
    selected = np.isfinite(x) & np.isfinite(z)
    if 'normal_m' in data:
        normal = np.asarray(data['normal_m'])
        available = np.unique(normal[selected & np.isfinite(normal)])
        if available.size:
            positive = available[available >= 0]
            layer = positive.min() if positive.size else available[np.argmin(np.abs(available))]
            selected &= np.isclose(normal, layer, rtol=0, atol=max(abs(layer)*1e-8, 1e-13))
    ids = np.flatnonzero(selected)
    # Different normal layers can project to duplicate positions; never average
    # their derived temperatures. The explicitly selected layer supplies one row.
    if ids.size:
        _, first = np.unique(np.column_stack((x[ids], z[ids])), axis=0, return_index=True)
        ids = ids[np.sort(first)]
    return ids


def _plane_figure(output, data, metadata, key, name, ylabel, *, log=False,
                  vmin=None, vmax=None):
    import matplotlib.tri as mtri
    from matplotlib.colors import LogNorm, Normalize
    from matplotlib.cm import ScalarMappable

    plt = _pyplot()
    required = (key, key+'_analytical', key+'_virtual')
    if any(k not in data for k in required):
        return []
    ids = _plane_layer(data)
    if len(ids) < 3:
        return []
    x, z = np.asarray(data['x_over_D'])[ids], np.asarray(data['z_over_D'])[ids]
    try:
        triangulation = mtri.Triangulation(x, z)
    except (RuntimeError, ValueError):
        return []
    if not len(triangulation.triangles):
        return []
    arrays = [np.asarray(data[k], dtype=float)[ids] for k in required]
    finite = np.concatenate([values[np.isfinite(values)] for values in arrays])
    if vmin is None:
        vmin = float(finite.min()) if finite.size else 0.
    if vmax is None:
        vmax = float(finite.max()) if finite.size else 1.
    if vmax <= vmin:
        vmax = vmin+1.
    levels = np.geomspace(vmin, vmax, 17) if log else np.linspace(vmin, vmax, 17)
    norm = LogNorm(vmin=vmin, vmax=vmax) if log else Normalize(vmin=vmin, vmax=vmax)
    labels = _labels(metadata, 'plane')
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.9), sharex=True, sharey=True)
    excluded = []
    for ax, values, label in zip(axes, arrays, labels):
        bad = ~np.isfinite(values)
        if log:
            bad |= values <= 0
        excluded.append(int(bad.sum()))
        tri = mtri.Triangulation(x, z, triangles=triangulation.triangles)
        tri.set_mask(np.any(bad[tri.triangles], axis=1))
        if not np.all(tri.mask):
            # Nonfinite values are only replaced in the display copy; every
            # adjacent triangle is masked, so that replacement is never drawn.
            display = np.where(bad, vmin, values)
            ax.tricontourf(tri, display, levels=levels, norm=norm, cmap='viridis', extend='both')
        else:
            ax.text(.5, .5, 'No displayable values', transform=ax.transAxes,
                    ha='center', va='center')
        ax.set_title(label, fontsize=9)
        ax.set_xlabel('$X/D$')
        ax.set_aspect('equal', adjustable='box')
    axes[0].set_ylabel('$z/D$')
    moment_bad = int((~_valid(data, len(np.asarray(data[key])))[ids]).sum())
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap='viridis'), ax=list(axes),
                           fraction=.025, pad=.025)
    colorbar.set_label(ylabel)
    # Exclude colorbar from tight_layout's automatic panel layout.
    colorbar.ax.set_in_layout(False)
    note = ('Display masks (DSMC, Cai, virtual): '+str(excluded)+
            (' nonpositive/nonfinite values.' if log else ' nonfinite values.')+
            f' Virtual moment-invalid rows: {moment_bad}/{len(ids)}. Raw values remain in the CSV/NPZ.')
    _finish(fig, metadata, ylabel+' | common contour levels', note, colorbar=True)
    colorbar.ax.set_position([.91, .24, .017, .49])
    path = Path(output)/(name+'_three_way.png')
    fig.savefig(path, dpi=160)
    plt.close(fig)

    differences = [arrays[1]-arrays[0], arrays[2]-arrays[0]]
    finite_error = np.concatenate([a[np.isfinite(a)] for a in differences])
    extent = max(float(np.max(np.abs(finite_error))) if finite_error.size else 0., 1e-15)
    fig, axes = plt.subplots(1, 2, figsize=(10., 4.9), sharex=True, sharey=True)
    error_norm = Normalize(vmin=-extent, vmax=extent)
    for ax, values, label in zip(axes, differences, labels[1:]):
        bad = ~np.isfinite(values)
        tri = mtri.Triangulation(x, z, triangles=triangulation.triangles)
        tri.set_mask(np.any(bad[tri.triangles], axis=1))
        if not np.all(tri.mask):
            ax.tripcolor(tri, np.where(bad, 0., values), shading='gouraud',
                         norm=error_norm, cmap='coolwarm')
        else:
            ax.text(.5, .5, 'No finite errors', transform=ax.transAxes, ha='center')
        ax.set_title(label+' minus DSMC', fontsize=9)
        ax.set_xlabel('$X/D$')
        ax.set_aspect('equal', adjustable='box')
    axes[0].set_ylabel('$z/D$')
    colorbar = fig.colorbar(ScalarMappable(norm=error_norm, cmap='coolwarm'), ax=list(axes),
                           fraction=.03, pad=.025)
    colorbar.set_label('Signed '+ylabel+' error')
    colorbar.ax.set_in_layout(False)
    _finish(fig, metadata, 'Signed error maps | '+ylabel,
            'Common signed color range; nonfinite display holes counted above. '+note, colorbar=True)
    colorbar.ax.set_position([.91, .24, .017, .49])
    error_path = Path(output)/(name+'_errors.png')
    fig.savefig(error_path, dpi=160)
    plt.close(fig)
    return [path, error_path]


def _standard_error(data, key):
    aliases = {'n_over_n0': 'virtual_n_se', 'U_sqrt_beta0': 'virtual_U_se',
               'T_over_T0': 'virtual_T_se', 'uR': 'virtual_uR_se', 'F': 'virtual_F_se'}
    candidates = (key+'_virtual_se', 'virtual_'+key+'_se',
                  'se_'+key+'_virtual', aliases.get(key, ''))
    for candidate in candidates:
        if candidate and candidate in data:
            values = np.asarray(data[candidate], dtype=float)
            if values.ndim == 1:
                return values
    return None


def _profile(output, data, metadata, key, name, ylabel, *, angular=False, log=False):
    required = (key, key+'_analytical', key+'_virtual')
    xkey = 'theta_deg' if angular else 'x_over_D'
    if xkey not in data or any(k not in data for k in required):
        return []
    plt = _pyplot()
    x = np.asarray(data[xkey], dtype=float)
    if not len(x):
        return []
    fig, ax = plt.subplots(figsize=(7.7, 5.3))
    labels = _labels(metadata, 'angular' if angular else 'centerline')
    arrays = [np.asarray(data[k], dtype=float) for k in required]
    excluded = []
    for values, label, color, style in zip(arrays, labels,
        ('tab:red', 'black', 'tab:orange'), ('o', '-', '-')):
        bad = ~np.isfinite(x) | ~np.isfinite(values)
        if log:
            bad |= values <= 0
        excluded.append(int(bad.sum()))
        display = np.ma.masked_array(values, bad)
        ax.plot(x, display, style, label=label, color=color,
                markersize=3., linewidth=1.5)
    se = _standard_error(data, key)
    if se is not None and se.shape == arrays[2].shape:
        good = np.isfinite(se) & (se >= 0) & np.isfinite(arrays[2])
        if log:
            good &= arrays[2]-se > 0
        # Sampling-only bands are not DSMC confidence intervals, and lower
        # nonpositive bands are omitted rather than clipped onto a log axis.
        ax.fill_between(x, arrays[2]-se, arrays[2]+se, where=good,
                        color='tab:orange', alpha=.2, label='Virtual: one batch standard error')
    invalid = ~_valid(data, len(x))
    mark = invalid & np.isfinite(x) & np.isfinite(arrays[2])
    if log:
        mark &= arrays[2] > 0
    if mark.any():
        ax.plot(x[mark], arrays[2][mark], 'x', color='crimson', ms=6,
                label='Virtual moment-invalid')
    if log:
        ax.set_yscale('log')
    ax.set_xlabel(r'Polar angle $\theta$ from axis [deg]' if angular else '$X/D$')
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.25, which='both')
    ax.legend(fontsize=8, loc='best')
    note = ('Display-excluded rows (DSMC, Cai, virtual): '+str(excluded)+
            f'; virtual moment-invalid: {int(invalid.sum())}/{len(x)}. '
            'Raw values are retained; model bands show sampling uncertainty only; DSMC intervals unavailable.')
    if key == 'T_over_T0' and metadata.get('temperature_note'):
        note += ' '+str(metadata['temperature_note'])
    subtitle = ('Angular profiles at R/D=10' if angular else 'Centerline profiles')+' | '+ylabel
    _finish(fig, metadata, subtitle, note)
    path = Path(output)/(name+'_three_way.png')
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return [path]


def plot_case(output, plane, centerline, angular, metadata):
    """Write three-way comparison PNGs, returning their paths.

    Inputs are mappings of CSV column names to arrays. Missing quantities are
    skipped. No input arrays are mutated. A logarithmic density/pressure image
    needs a display mask for nonpositive values; its excluded count and all
    model realizability failures are printed explicitly, never repaired.
    """
    if _pyplot() is None:
        return []
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    written = []
    if plane is not None:
        written += plot_density_contours(output, plane, metadata)
        for key, name, label, log, low, high in (
            ('n_over_n0', 'density', r'$n/n_0$', True, 1e-4, 1.),
            ('ux', 'plane_velocity_axial', r'$u_x\sqrt{\beta_0}$', False, 0., 2.6),
            ('ur', 'plane_velocity_radial', r'$u_r\sqrt{\beta_0}$ (cylindrical)', False, 0., 2.6),
            ('T_over_T0', 'plane_temperature', r'$T/T_0$', False, 0., 1.1),
            ('p_over_p0', 'plane_pressure', r'$p/p_0$', True, 1e-5, 1.)):
            written += _plane_figure(output, plane, metadata, key, name, label,
                                     log=log, vmin=low, vmax=high)
    if centerline is not None:
        for key, name, label, log in (
            ('n_over_n0', 'centerline_density', r'$n/n_0$', True),
            ('U_sqrt_beta0', 'centerline_velocity', r'$u_x\sqrt{\beta_0}$', False),
            ('T_over_T0', 'centerline_temperature', r'Local-cell $T/T_0$ (legacy average)', False),
            ('p_over_p0', 'centerline_pressure', r'$p/p_0$', True)):
            written += _profile(output, centerline, metadata, key, name, label, log=log)
    if angular is not None:
        for key, name, label, log in (
            ('n_over_n0', 'angular_density', r'$n/n_0$', True),
            ('uR', 'angular_velocity', r'$u_R\sqrt{\beta_0}$ (spherical outward)', False),
            ('F', 'angular_flux', r'$F=n u_R R^2/(\Gamma_0\pi R_0^2)$ [1/sr]', False),
            ('T_over_T0', 'angular_temperature', r'Local-cell $T/T_0$ (legacy average)', False)):
            written += _profile(output, angular, metadata, key, name, label,
                                angular=True, log=log)
    return written
