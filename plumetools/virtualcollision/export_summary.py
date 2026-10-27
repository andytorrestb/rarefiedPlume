"""Readable index and four-Kn views of completed, saved three-method runs."""
from pathlib import Path
import json
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .assessment import write_json


def read_table(path):
    table = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True, dtype=None, encoding='utf-8'))
    return {key: np.asarray(table[key]) for key in table.dtype.names}


def summarize(study):
    study = Path(study).resolve()
    output = study/'results'/'virtual_collision'
    output.mkdir(parents=True, exist_ok=True)
    reports = []
    for path in sorted(study.rglob('results/virtual_collision/report.json')):
        report = json.loads(path.read_text(encoding='utf-8'))
        case = path.parents[2]
        report['case_path'] = str(case)
        report['output_path'] = str(path.parent)
        tables = {name: read_table(path.parent/(name+'_three_way.csv'))
                  for name in ('plane', 'centerline', 'angular')}
        from .export_figures import plot_density_contours
        new_contours = plot_density_contours(path.parent, tables['plane'], report)
        report['plots'] = list(dict.fromkeys(report.get('plots', [])+[str(p) for p in new_contours]))
        # Check that the known near-lip quadrature discrepancy does not leak
        # into the density metric's common mask. The exact saved signed delta
        # recovers the actual ballistic field used, without rerunning transport.
        grid = np.load(path.parent/'correction_grid_raw.npz')
        raw = np.load(path.parent/'plane_model_raw.npz')
        points = raw['points_m']
        r = np.hypot(points[:, 1], points[:, 2])
        delta = grid['batch_delta_M0'].mean(axis=0).reshape(len(grid['x_centres_m']), len(grid['r_centres_m']))
        interpolator = RegularGridInterpolator((grid['x_centres_m'], grid['r_centres_m']), delta)
        correction = interpolator(np.column_stack((points[:, 0], np.maximum(r, grid['r_centres_m'][0]))))
        # Exit density follows the configured Kn convention, not DSMC values.
        from .export_geometry import load_export_context
        context = load_export_context(case)
        n0 = context['source'].number_density_m3
        baseline = (raw['M0']-correction)/n0
        expected = tables['plane']['n_over_n0_analytical']
        common = ~tables['plane']['nozzle_lip_zone'].astype(bool)
        common &= expected >= .001
        report['plane_cai_density_outside_lip'] = dict(
            common_points=int(common.sum()),
            max_absolute=float(np.max(np.abs(baseline[common]-expected[common]))),
            normalized_L2=float(np.linalg.norm(baseline[common]-expected[common])/np.linalg.norm(expected[common])),
            note='Original Cai columns retained; meridian disk quadrature has an additional localized near-lip azimuth error.')
        write_json(path, report)
        flags = report['quality']
        pointer = case/'results'/'THREE_METHOD_COMPARISON.md'
        pointer.write_text('# DSMC, Cai and virtual-collision comparison\n\n'+flags+'\n\n'
            '[All tables, model raw arrays and diagnostics](virtual_collision/README.md).\n\n'
            '![Density contour comparison](virtual_collision/density_contour_comparison.png)\n\n'
            '![Density contours](virtual_collision/density_three_way.png)\n\n'
            '![Centerline density](virtual_collision/centerline_density_three_way.png)\n\n'
            '![Angular flux](virtual_collision/angular_flux_three_way.png)\n', encoding='utf-8')
        reports.append((report, tables))
    reports.sort(key=lambda pair: -pair[0]['kn'])
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None
    if plt is not None and reports:
        from .export_figures import _density_contours_on_axes
        contour_rows = int(np.ceil(len(reports)/2))
        fig, contour_axes = plt.subplots(contour_rows, 2, figsize=(13, 6*contour_rows), squeeze=False)
        legend_handles = []
        for ax, (report, tables) in zip(contour_axes.ravel(), reports):
            status, legend_handles = _density_contours_on_axes(ax, tables['plane'], report)
            short = 'non-production DSMC' if report['kn'] < .01 else 'stationarity inconclusive'
            if report['diagnostics']['multiple_event_warning']:
                short += '; one-event exposure exceeded'
            ax.set_title(f'Kn={report["kn"]:g}\n{short}', fontsize=9)
        for ax in contour_axes.ravel()[len(reports):]:
            ax.set_visible(False)
        if legend_handles:
            legend_handles[0].set_label('DSMC (dsmcFoam)')
        fig.legend(handles=legend_handles, loc='upper center', ncol=3, bbox_to_anchor=(.5, .97), fontsize=9)
        fig.suptitle('Density contour comparison | n/n0 = 0.1, 0.01, 0.001', y=.997, fontsize=12)
        fig.text(.5, .008, 'Labels mark Cai levels. Shared display masks; raw nonphysical values and reference/model validity flags are retained in each case report.',
                 ha='center', fontsize=8)
        fig.tight_layout(rect=(0, .03, 1, .93))
        fig.savefig(output/'density_contours_all_kn_three_way.png', dpi=180)
        fig.savefig(output/'density_contours_all_kn_three_way.svg')
        plt.close(fig)
        for kind, coordinate, fields, labels in (
            ('centerline', 'x_over_D', ('n_over_n0', 'U_sqrt_beta0', 'T_over_T0'),
             ('Density n/n0', 'Axial velocity U sqrt(beta0)', 'Local-cell T/T0 (legacy average)')),
            ('angular', 'theta_deg', ('n_over_n0', 'uR', 'F'),
             ('Density n/n0', 'Spherical outward velocity sqrt(beta0)', 'Particle flux / sr'))):
            fig, axes = plt.subplots(len(reports), 3, figsize=(14, 3.1*len(reports)), squeeze=False)
            for row, (report, tables) in enumerate(reports):
                data = tables[kind]
                x = data[coordinate]
                for ax, key, label in zip(axes[row], fields, labels):
                    for suffix, style, color, legend in (
                        ('', '.', 'tab:red', 'DSMC'), ('_analytical', '-', 'black', 'Cai collisionless'),
                        ('_virtual', '-', 'tab:orange', 'Virtual-collision model')):
                        y = np.asarray(data[key+suffix], dtype=float)
                        show = np.isfinite(y) & ((y > 0) if key == 'n_over_n0' else True)
                        ax.plot(x[show], y[show], style, color=color, ms=3, lw=1.2, label=legend)
                    y = data[key+'_virtual']
                    se = data.get(key+'_virtual_se')
                    if se is not None:
                        low, high = y-se, y+se
                        show = np.isfinite(low) & np.isfinite(high)
                        if key == 'n_over_n0':
                            show &= low > 0
                        ax.fill_between(x, low, high, where=show, color='tab:orange', alpha=.15)
                    bad = ~data['virtual_valid'].astype(bool) & np.isfinite(y)
                    if key == 'n_over_n0':
                        bad &= y > 0
                        ax.set_yscale('log')
                    ax.plot(x[bad], y[bad], 'rx', ms=4)
                    ax.set(xlabel='X/D' if kind == 'centerline' else 'theta [deg]', ylabel=label)
                    ax.grid(alpha=.2)
                short = 'non-production reference' if report['kn'] < .01 else 'stationarity inconclusive'
                if report['diagnostics']['multiple_event_warning']:
                    short += '; one-event exposure exceeded'
                axes[row, 0].set_title(f'Kn={report["kn"]:g}: {short}', fontsize=9, loc='left')
            handles, legends = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles[:3], legends[:3], loc='upper center', ncol=3, bbox_to_anchor=(.5, .985))
            fig.suptitle('Exploratory three-method comparison | '+kind, y=.999, fontsize=12)
            fig.text(.5, .005, 'Shading: one model batch standard error; crosses: moment-invalid. '
                'Raw values retained. Numerical smoothing/convergence and DSMC uncertainty remain unresolved.',
                ha='center', fontsize=8)
            fig.tight_layout(rect=(0, .025, 1, .955))
            fig.savefig(output/(kind+'_all_kn_three_way.png'), dpi=160)
            plt.close(fig)
    lines = ['# All extension DSMC/Cai/virtual-collision comparisons', '',
        'All four available direct-nozzle result sets have been processed. Original DSMC/Cai files are unchanged.', '',
        '[Usefulness review and refinement plan](review/README.md).', '',
        '![Overlaid density contours](density_contours_all_kn_three_way.png)', '',
        '![Centerline overview](centerline_all_kn_three_way.png)', '',
        '![Angular overview](angular_all_kn_three_way.png)', '',
        '| Kn | Case | Three-method plots |', '|---:|---|---|']
    for report, _ in reports:
        case = Path(report['case_path'])
        link = Path(os.path.relpath(case/'results'/'THREE_METHOD_COMPARISON.md', output)).as_posix()
        lines.append(f'| {report["kn"]:g} | {case.parent.name}/{case.name} | [Open]({link}) |')
    lines += ['', 'Each case contains 19 PNG plots plus the density-contour SVG, augmented profile/plane CSVs, raw model moments/tensors, '
        'batch correction grids, metrics and validity/conservation diagnostics.', '',
        'Missing: new uniform `Cases/Kn10_np1x` and `Cases/Kn1_np1x` have no DSMC comparison CSVs. '
        'The available archived graded runs above were processed instead. DSMC tensor components and '
        'independent instantaneous time blocks are also absent.', '',
        'All comparisons are exploratory. Kn 0.001 and 0.0001 DSMC references remain non-production; '
        'their one-event predictions are outside the modeled exposure range. Kn 1/10 stationarity is inconclusive. '
        'Invalid predictions are preserved and flagged, not replaced with Cai.', '',
        'The correction is a finite-ring raw-moment interpolant added to Cai at verified original cell locations. '
        'Centerline/angular plots preserve the legacy local-cell temperature average; '
        'raw-aggregate model temperatures are separately stored. See each report for geometric and Cai-reproduction checks.', '']
    (output/'README.md').write_text('\n'.join(lines), encoding='utf-8')
    write_json(output/'completed_cases.json', [report for report, _ in reports])
    return len(reports)
