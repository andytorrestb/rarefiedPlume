#!/usr/bin/env python3
"""Review saved comparisons; never rerun transport or alter predictions."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import sys
import json
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from plumetools.virtualcollision.export_summary import read_table
from plumetools.virtualcollision.export_geometry import load_export_context, centerline_geometry
from plumetools.virtualcollision.assessment import dsmc_reference_metrics, write_json, write_rows

STUDY = ROOT/'cases'/'cai2012-kn-extension'
OUTPUT = STUDY/'results'/'virtual_collision'/'review'


def score(table, field, weights, common, *, scale, bounds=None, coordinate=None):
    mask = np.asarray(common).copy()
    if bounds is not None:
        mask &= (coordinate >= bounds[0]) & (coordinate < bounds[1])
    ref, base, vc = [np.asarray(table[field+suffix], float) for suffix in ('', '_analytical', '_virtual')]
    result = dsmc_reference_metrics(ref, base, vc, weights=weights, common_mask=mask, scale=scale)
    valid = np.asarray(table['virtual_valid'], bool)
    result['invalid_model_on_common_mask'] = int((mask & ~valid).sum())
    result['headline_eligible'] = False
    if mask.any():
        w = weights[mask]
        def rms(values):
            values = np.asarray(values)[mask]
            return float(np.sqrt(np.sum(w*values**2)/w.sum())) if np.isfinite(values).all() else None
        result['correction_RMS'] = rms(vc-base)
        se = table.get(field+'_virtual_se')
        if se is not None:
            se = np.asarray(se, float)
            result['reported_model_SE_RMS'] = rms(se)
            result['correction_to_reported_SE_RMS'] = (result['correction_RMS']/result['reported_model_SE_RMS']
                if result['reported_model_SE_RMS'] else None)
            result['zero_reported_SE_weight_fraction'] = float(weights[mask & (se == 0)].sum()/w.sum())
            improvement = np.abs(base-ref)-np.abs(vc-ref)
            result['weight_fraction_gap_gt_2_model_SE_and_realizable'] = float(
                weights[mask & valid & np.isfinite(se) & (improvement > 2*se)].sum()/w.sum())
        else:
            result['reported_model_SE_RMS'] = None
        result['improved_weight_fraction'] = float(weights[mask & (np.abs(vc-ref) < np.abs(base-ref))].sum()/w.sum())
        result['first_point_fraction_of_baseline_squared_error'] = None
        indices = np.flatnonzero(mask)
        squared = weights[indices]*(base[indices]-ref[indices])**2
        if squared.sum() > 0:
            result['first_point_fraction_of_baseline_squared_error'] = float(squared[0]/squared.sum())
    return result


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    global_rows, region_rows, diagnostic_rows = [], [], []
    for path in sorted(STUDY.rglob('results/virtual_collision/report.json')):
        report = json.loads(path.read_text(encoding='utf-8'))
        case = path.parents[2]
        kn = float(report['kn'])
        context = load_export_context(case)
        geometry = centerline_geometry(context)
        diagnostics = report['diagnostics']
        remaining = np.asarray(diagnostics['remaining_exposure'])
        pair_weights = np.repeat(np.asarray(diagnostics['symmetric_pair_weights']), 2)
        grid = np.load(path.parent/'correction_grid_raw.npz')
        nonzero_batches = np.count_nonzero(grid['batch_delta_M0'], axis=0)
        diagnostic_rows.append(dict(kn=kn, quality=report['quality'],
            exposure_mean=diagnostics['exposure_mean'],
            omitted_opportunity_proxy=float(np.average(-np.expm1(-remaining), weights=pair_weights)),
            weighted_remaining_above_one=float(np.average(remaining > 1, weights=pair_weights)),
            uncollided_fraction=diagnostics['ensemble_uncollided_source_number_fraction'],
            grid_bins=len(nonzero_batches), all_zero_batch_corrections=int((nonzero_batches == 0).sum()),
            nonzero_correction_in_at_most_two_batches=int((nonzero_batches <= 2).sum()),
            caveat='nonzero batch corrections measure estimator support, not literal event counts; zero SE can be unresolved zero-hit sampling'))
        for kind in ('plane', 'centerline', 'angular'):
            table = read_table(path.parent/(kind+'_three_way.csv'))
            if kind == 'plane':
                weights = np.load(path.parent/'plane_model_raw.npz')['cell_volumes_m3']
                common = (table['n_over_n0_analytical'] >= .001) & ~table['nozzle_lip_zone'].astype(bool)
                fields = ('n_over_n0', 'ux', 'ur', 'T_over_T0')
                coordinate = table['x_over_D']
            elif kind == 'centerline':
                weights = geometry['volume_sum_m3'][geometry['export_bin_index']]
                common = table['n_over_n0_analytical'] >= .001
                fields = ('n_over_n0', 'U_sqrt_beta0', 'T_over_T0')
                coordinate = table['x_over_D']
            else:
                weights = table['solid_angle_sr']
                common = table['n_over_n0_analytical'] >= .001*np.nanmax(table['n_over_n0_analytical'])
                fields = ('n_over_n0', 'uR', 'F', 'T_over_T0')
                coordinate = table['theta_deg']
            regions = ((0, .5), (.5, 2), (2, 5), (5, 10.000001)) if kind != 'angular' else (
                (0, 15), (15, 30), (30, 60), (60, 90.000001))
            for field in fields:
                if field+'_virtual' not in table:
                    continue
                scale = None if field in ('n_over_n0', 'F') else 1.
                record = score(table, field, weights, common, scale=scale)
                record.update(kn=kn, comparison=kind, observable=field, quality=report['quality'])
                global_rows.append(record)
                for bounds in regions:
                    row = score(table, field, weights, common, scale=scale, bounds=bounds, coordinate=coordinate)
                    row.update(kn=kn, comparison=kind, observable=field,
                        region=f'{bounds[0]:g}-{min(bounds[1], 90 if kind == "angular" else 10):g}',
                        region_coordinate='theta_deg' if kind == 'angular' else 'X/D', quality=report['quality'])
                    region_rows.append(row)
    write_json(OUTPUT/'metrics_review.json', global_rows)
    write_rows(OUTPUT/'metrics_review.csv', global_rows)
    write_json(OUTPUT/'regional_review.json', region_rows)
    write_rows(OUTPUT/'regional_review.csv', region_rows)
    write_rows(OUTPUT/'validity_review.csv', diagnostic_rows)
    write_json(OUTPUT/'validity_review.json', diagnostic_rows)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        comparisons = [('plane', 'n_over_n0', 'Plane density'),
            ('centerline', 'n_over_n0', 'Centerline density'),
            ('centerline', 'U_sqrt_beta0', 'Centerline axial velocity'),
            ('centerline', 'T_over_T0', 'Centerline local-cell temperature'),
            ('angular', 'n_over_n0', 'Angular density'), ('angular', 'uR', 'Angular outward velocity'),
            ('angular', 'F', 'Angular particle flux'), ('angular', 'T_over_T0', 'Angular local-cell temperature')]
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), sharey=True)
        for ax, kn in zip(axes, (10., 1.)):
            for index, (kind, observable, label) in enumerate(comparisons):
                row = next(r for r in global_rows if r['kn'] == kn and r['comparison'] == kind and r['observable'] == observable)
                ratio = row['corrected_error']/row['baseline_error']
                ax.plot([1, ratio], [index, index], color='tab:green' if ratio < 1 else 'tab:red', lw=2)
                ax.scatter([ratio], [index], color='tab:green' if ratio < 1 else 'tab:red', marker='x' if row['invalid_model_on_common_mask'] else 'o')
                ax.text(ratio*1.06, index, f'{ratio:.2f}x', va='center', fontsize=8)
            ax.axvline(1, color='black', ls='--', lw=1)
            ax.set(xscale='log', xlim=(.5, 17), xlabel='Virtual error / Cai error (smaller is better)', title=f'Kn={kn:g}')
            ax.grid(axis='x', alpha=.2)
        axes[0].set_yticks(np.arange(len(comparisons)), [item[2] for item in comparisons])
        axes[0].invert_yaxis()
        fig.suptitle('Current aggregate errors: selected improvements coexist with broader degradation', fontsize=11)
        fig.text(.5, .01, 'x: non-realizable model contributions on common mask. These point estimates omit DSMC confidence intervals; numerical convergence is unresolved.', ha='center', fontsize=8)
        fig.tight_layout(rect=(0, .04, 1, .94))
        fig.savefig(OUTPUT/'error_ratios.png', dpi=170)
        plt.close(fig)
    except ImportError:
        pass
    print(f'Reviewed {len(diagnostic_rows)} cases, {len(global_rows)} aggregate scores and {len(region_rows)} predeclared diagnostic region scores.')


if __name__ == '__main__':
    main()
