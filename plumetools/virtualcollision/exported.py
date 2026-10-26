"""Three-way comparisons from verified, compact Cai-family result exports.

The original DSMC/Cai columns are preserved. Only source parameters, geometry,
and numerical controls enter the predictive correction. The signed residence
correction is interpolated in raw moments; exact Cai raw moments supply the
ballistic point field. This numerical smoothing is reported, never disguised as
a new physical closure or fitted to DSMC.
"""
from __future__ import annotations
from dataclasses import asdict
from pathlib import Path
import math
import platform
import time
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from plumetools.cai2012 import analytical
from plumetools.cai2012ext import comparison
from .assessment import inventory, write_json, write_rows, dsmc_reference_metrics
from .export_geometry import load_export_context, meridional_geometry, centerline_geometry, angular_geometry
from .model import NumericalOptions, VirtualCollisionModel, _cai_raw
from .kernel import VHSKernel
from .moments import derive


def _table(values):
    return {key: np.asarray(values[key]) for key in values.dtype.names}


def _pack(raw):
    return np.concatenate((raw['M0'][..., None], raw['M1'],
                           raw['M2'].reshape(raw['M2'].shape[:-2]+(9,))), axis=-1)


def _derive(raw, mass):
    return derive(raw[..., 0], raw[..., 1:4], raw[..., 4:].reshape(raw.shape[:-1]+(3, 3)), mass)


class BallisticCache:
    """Reuse the source-only, normalized exact Cai tensors across Kn cases."""
    def __init__(self):
        self.values = {}

    def evaluate(self, points, source, options):
        points = np.asarray(points)
        meridian = np.column_stack((points[:, 0], np.hypot(points[:, 1], points[:, 2])))
        coordinates, inverse = np.unique(np.round(meridian, 12), axis=0, return_inverse=True)
        signature = (source.radius_m, source.temperature_K, source.speed_ratio,
                     source.species.mass_kg, options.cai_radial, options.cai_azimuthal)
        cache = self.values.setdefault(signature, {})
        absent = np.array([i for i, q in enumerate(coordinates) if tuple(q) not in cache], dtype=int)
        if len(absent):
            p = coordinates[absent]
            threshold = max(2*math.pi*source.radius_m/options.cai_radial,
                            4*math.pi*source.radius_m/options.cai_azimuthal)
            floor = 4*math.pi*source.radius_m/analytical.MAX_AZIMUTHAL
            edges = [math.inf, threshold]
            while edges[-1]/2 > floor:
                edges.append(edges[-1]/2)
            edges.append(0.)
            for upper, lower in zip(edges[:-1], edges[1:]):
                ids = np.flatnonzero((p[:, 0] > lower) & (p[:, 0] <= upper))
                if not len(ids):
                    continue
                raw = _cai_raw(np.column_stack((p[ids], np.zeros(len(ids)))), source, options)
                normalized = _pack(raw)/source.number_density_m3
                for index, value in zip(ids, normalized):
                    cache[tuple(p[index])] = value
        result = np.stack([cache[tuple(q)] for q in coordinates])[inverse]*source.number_density_m3
        return result


def correction_grid(diameter, *, refinement=1):
    """Fixed numerical grid independent of DSMC field values."""
    x = np.array([0, .025, .05, .1, .18, .3, .45, .65, .9, 1.2, 1.6, 2.1,
                  2.7, 3.4, 4.2, 5.1, 6.1, 7.2, 8.4, 9.6, 10.8])
    r = np.array([0, .025, .05, .09, .15, .23, .34, .48, .65, .85, 1.1, 1.4,
                  1.8, 2.3, 2.9, 3.6, 4.4, 5.3, 6.3, 7.4, 8.6, 10., 11.5])
    for _ in range(int(refinement)):
        x = np.sort(np.r_[x, .5*(x[1:]+x[:-1])])
        r = np.sort(np.r_[r, .5*(r[1:]+r[:-1])])
    x, r = x*diameter, r*diameter
    xc = .5*(x[1:]+x[:-1])
    rc = np.sqrt(.5*(r[1:]**2+r[:-1]**2))
    xb = np.repeat(np.column_stack((x[:-1], x[1:])), len(rc), axis=0)
    rb = np.tile(np.column_stack((r[:-1], r[1:])), (len(xc), 1))
    return xc, rc, xb, rb


class PointPrediction:
    def __init__(self, model, grid, ballistic_cache):
        self.model, self.source, self.cache = model, model.source, ballistic_cache
        self.x, self.r, self.xb, self.rb = grid
        self.correction = model.raw_correction_bins(self.xb, self.rb)
        batches = _pack(self.correction['batch_raw_moments'])
        self.interpolators = [RegularGridInterpolator((self.x, self.r),
            batch.reshape(len(self.x), len(self.r), 13), bounds_error=True) for batch in batches]

    def evaluate(self, points):
        points = np.asarray(points)
        base = self.cache.evaluate(points, self.source, self.model.options)
        radius = np.hypot(points[:, 1], points[:, 2])
        # Below the first ring centre, use even continuation and the analytic
        # axis limits for odd components. This is an explicit symmetry boundary
        # rule for the numerical interpolant, not a positivity repair.
        query = np.column_stack((points[:, 0], np.maximum(radius, self.r[0])))
        deltas = np.stack([interpolator(query) for interpolator in self.interpolators])
        inner = radius < self.r[0]
        if np.any(inner):
            factor = radius[inner]/self.r[0]
            deltas[:, inner, 2:4] *= factor[None, :, None]
            tensor = deltas[:, inner, 4:].reshape(len(deltas), int(inner.sum()), 3, 3).copy()
            transverse = .5*(tensor[..., 1, 1]+tensor[..., 2, 2])
            f2 = factor[None, :]**2
            tensor[..., 1, 1] = transverse+f2*(tensor[..., 1, 1]-transverse)
            tensor[..., 2, 2] = transverse+f2*(tensor[..., 2, 2]-transverse)
            for j, k in ((0, 1), (1, 0), (0, 2), (2, 0)):
                tensor[..., j, k] *= factor[None, :]
            tensor[..., 1, 2] *= f2
            tensor[..., 2, 1] *= f2
            deltas[:, inner, 4:] = tensor.reshape(len(deltas), int(inner.sum()), 9)
        raw_batches = base[None, :, :]+deltas
        mean = _derive(raw_batches.mean(axis=0), self.source.species.mass_kg)
        batch_fields = [_derive(batch, self.source.species.mass_kg) for batch in raw_batches]
        return dict(fields=mean, raw_batches=raw_batches, batch_fields=batch_fields,
                    ballistic=_derive(base, self.source.species.mass_kg),
                    inner_symmetry_points=int(inner.sum()))


def _standard_error(values):
    return np.std(np.asarray(values), axis=0, ddof=1)/math.sqrt(len(values))


def _point_columns(prediction, source):
    f = prediction['fields']
    output = dict(n_over_n0_virtual=f['n']/source.number_density_m3,
        ux_virtual=f['U'][:, 0]/source.thermal_speed,
        ur_virtual=f['U'][:, 1]/source.thermal_speed,
        uphi_virtual=f['U'][:, 2]/source.thermal_speed,
        T_over_T0_virtual=f['T_trans']/source.temperature_K,
        p_over_p0_virtual=f['n']/source.number_density_m3*f['T_trans']/source.temperature_K,
        virtual_valid=f['valid'])
    for key, values in (
        ('n_over_n0', [b['n']/source.number_density_m3 for b in prediction['batch_fields']]),
        ('ux', [b['U'][:, 0]/source.thermal_speed for b in prediction['batch_fields']]),
        ('ur', [b['U'][:, 1]/source.thermal_speed for b in prediction['batch_fields']]),
        ('T_over_T0', [b['T_trans']/source.temperature_K for b in prediction['batch_fields']])):
        output[key+'_virtual_se'] = _standard_error(values)
    return output


def _profile(prediction, geometry, source, *, angular=False):
    """The exact legacy weights, with separate raw-aggregate central T."""
    volume = geometry['volumes_m3']
    if angular:
        weights = geometry['fractions']*volume[:, None]
        def total(value):
            return weights.T @ np.asarray(value)
        selection = slice(None)
        ct, st = np.cos(geometry['theta_rad']), np.sin(geometry['theta_rad'])
    else:
        def total(value):
            value = np.asarray(value)
            if value.ndim == 1:
                return np.bincount(geometry['bin_index'], weights=volume*value, minlength=geometry['n_bins'])
            return np.column_stack([total(value[:, j]) for j in range(value.shape[1])])
        selection = geometry['export_bin_index']
    vol = total(np.ones(len(volume)))
    def reduce(fields):
        n, u, t = fields['n'], fields['U'], fields['T_trans']
        molecules = total(n)
        with np.errstate(divide='ignore', invalid='ignore'):
            out = dict(n_over_n0=total(n)/vol/source.number_density_m3,
                T_over_T0=total(n*t)/molecules/source.temperature_K,
                p_over_p0=total(n*t)/vol/source.number_density_m3/source.temperature_K)
            if angular:
                ur = ct*u[:, 0]+st*u[:, 1]
                ut = -st*u[:, 0]+ct*u[:, 1]
                out.update(uR=total(n*ur)/molecules/source.thermal_speed,
                    utheta=total(n*ut)/molecules/source.thermal_speed,
                    uphi=total(n*u[:, 2])/molecules/source.thermal_speed,
                    F=total(comparison.mass_flux_per_steradian(n/source.number_density_m3,
                        ur/source.thermal_speed, geometry['radius_m'], source.radius_m, source.speed_ratio))/vol)
            else:
                out['U_sqrt_beta0'] = total(fields['M1'][:, 0])/molecules/source.thermal_speed
            # Raw trace plus mean vector retains between-cell dispersion.
            mean_vector = total(fields['M1'])/molecules[:, None]
            raw_square = total(np.trace(fields['M2'], axis1=1, axis2=2))/molecules
            out['T_over_T0_raw_aggregated'] = (raw_square-np.sum(mean_vector**2, axis=1))/(
                3*source.species.gas_constant_j_per_kg_k*source.temperature_K)
        for key in out:
            out[key] = out[key][selection]
        invalid = total((~fields['valid']).astype(float)) > 0
        out['virtual_valid'] = (~invalid)[selection] & np.isfinite(out['n_over_n0']) & (out['n_over_n0'] > 0)
        return out
    mean = reduce(prediction['fields'])
    baseline = reduce(prediction['ballistic'])
    batches = [reduce(fields) for fields in prediction['batch_fields']]
    columns = {key+'_virtual': value for key, value in mean.items() if key != 'virtual_valid'}
    columns['virtual_valid'] = mean['virtual_valid']
    for key in mean:
        if key != 'virtual_valid':
            columns[key+'_virtual_se'] = _standard_error([b[key] for b in batches])
    return columns, baseline


def _verification(original, rebuilt, pairs):
    output = {}
    for old, new in pairs:
        if old not in original or new not in rebuilt:
            continue
        a, b = np.asarray(original[old]), np.asarray(rebuilt[new])
        output[old] = dict(max_absolute=float(np.max(np.abs(a-b))),
            normalized_L2=float(np.linalg.norm(a-b)/max(np.linalg.norm(a), 1e-300)))
    return output


def _write_table(path, table):
    write_rows(path, [{key: np.asarray(value)[i].item() for key, value in table.items()}
                      for i in range(len(next(iter(table.values()))))])


def run_case(case, options, cache, *, refinement=1, plots=True, quality='exploratory'):
    start, cpu = time.perf_counter(), time.process_time()
    context = load_export_context(case)
    source = context['source']
    model = VirtualCollisionModel(source, VHSKernel(source.species), options)
    model.preprocess()
    grid = correction_grid(2*source.radius_m, refinement=refinement)
    predictor = PointPrediction(model, grid, cache)
    output = Path(case)/'results'/'virtual_collision'
    output.mkdir(parents=True, exist_ok=True)
    tables, checks, point_flags = {}, {}, {}
    for name, geometry in (('plane', meridional_geometry(case, context)),
                            ('centerline', centerline_geometry(context)),
                            ('angular', angular_geometry(context))):
        print(f'  {name}: {len(geometry["points_m"]):,} reconstructed cells', flush=True)
        prediction = predictor.evaluate(geometry['points_m'])
        original = _table(geometry['table'])
        if name == 'plane':
            columns = _point_columns(prediction, source)
            b = prediction['ballistic']
            rebuilt = dict(n=b['n']/source.number_density_m3, ux=b['U'][:, 0]/source.thermal_speed,
                           T=b['T_trans']/source.temperature_K)
            checks[name] = _verification(original, rebuilt,
                [('n_over_n0_analytical', 'n'), ('ux_analytical', 'ux'), ('T_over_T0_analytical', 'T')])
        else:
            columns, rebuilt = _profile(prediction, geometry, source, angular=name == 'angular')
            checks[name] = _verification(original, rebuilt,
                [('n_over_n0_analytical', 'n_over_n0'), ('T_over_T0_analytical', 'T_over_T0'),
                 ('U_sqrt_beta0_analytical', 'U_sqrt_beta0'), ('uR_analytical', 'uR'), ('F_analytical', 'F')])
        table = {**original, **columns}
        tables[name] = table
        _write_table(output/(name+'_three_way.csv'), table)
        fields = prediction['fields']
        np.savez_compressed(output/(name+'_model_raw.npz'), points_m=geometry['points_m'],
            cell_volumes_m3=geometry['volumes_m3'], **{key: fields[key] for key in
                ('M0', 'M1', 'M2', 'U', 'T_tensor', 'T_trans', 'P_tensor', 'valid')})
        point_flags[name] = dict(cells=len(fields['n']), invalid=int((~fields['valid']).sum()),
                                inner_symmetry_points=prediction['inner_symmetry_points'])
    np.savez_compressed(output/'correction_grid_raw.npz', x_centres_m=grid[0], r_centres_m=grid[1],
        x_bounds_m=grid[2], r_bounds_m=grid[3],
        **{'batch_delta_'+key: values for key, values in predictor.correction['batch_raw_moments'].items()})
    warning = model.diagnostics['multiple_event_warning']
    metadata = dict(case_name=Path(case).name, case_path=str(case), kn=context['cfg'].exit.knudsen,
        quality=quality+('; ONE-EVENT EXPOSURE LIMIT EXCEEDED' if warning else ''),
        estimator_note='Verified original cell geometry; finite ring-grid raw correction',
        temperature_note='Profiles preserve legacy averages of local cell T; raw-aggregate T is stored separately.',
        numerical_options=asdict(options), grid_refinement=refinement, grid_bins=len(grid[2]),
        geometry_verification=context['verification'], ballistic_export_reproduction=checks,
        point_validity=point_flags, diagnostics=model.diagnostics,
        numerical_uncertainty='independent batch SE propagated through interpolation and legacy aggregation; quadrature/domain/smoothing error additional',
        missing_data=['DSMC tensor components', 'independent instantaneous time blocks for DSMC uncertainty'],
        limitation='Comparisons exploratory: reference quality and numerical convergence not established. Raw arrays preserved; no clipping or fallback.',
        axis_rule='Below first radial node: even continuation; odd components vanish at axis; transverse tensor equality at axis.',
        wall_s=time.perf_counter()-start, core_s=time.process_time()-cpu, platform=platform.platform())
    metrics = []
    for name, table in tables.items():
        for field in ('n_over_n0', 'U_sqrt_beta0', 'uR', 'F', 'T_over_T0'):
            if field not in table or field+'_virtual' not in table:
                continue
            n = table['n_over_n0_analytical']
            floor = comparison.DENSITY_FLOOR_OVER_N0 if name != 'angular' else comparison.ANGULAR_FLOOR_FRACTION*np.nanmax(n)
            common = np.isfinite(table[field]) & (n >= floor)
            if name == 'plane':
                common &= ~table.get('nozzle_lip_zone', np.zeros(len(n), dtype=bool)).astype(bool)
                weights = meridional_geometry(case, context)['volumes_m3']
            elif name == 'centerline':
                cg = centerline_geometry(context)
                weights = cg['volume_sum_m3'][cg['export_bin_index']]
            else:
                weights = table['solid_angle_sr']
            score = dsmc_reference_metrics(table[field], table[field+'_analytical'], table[field+'_virtual'],
                weights=weights, common_mask=common, scale=None if field in ('n_over_n0', 'F') else 1.)
            score.update(comparison=name, observable=field, headline_eligible=False,
                         invalid_virtual_on_common_mask=int(np.count_nonzero(common & ~table['virtual_valid'])))
            metrics.append(score)
    write_json(output/'metrics.json', metrics)
    write_rows(output/'metrics.csv', metrics)
    if plots:
        from .export_figures import plot_case
        metadata['plots'] = [str(p) for p in plot_case(output, tables['plane'], tables['centerline'], tables['angular'], metadata)]
    write_json(output/'report.json', metadata)
    (output/'README.md').write_text(
        '# Three-method comparison\n\n'+metadata['quality']+'\n\n'+metadata['estimator_note']+'.\n\n'+
        metadata['temperature_note']+'\n\n'+metadata['limitation']+'\n\n'+
        'DSMC tensors and independent time-block uncertainty are unavailable. '
        'See report.json for mesh checks, Cai reproduction, sampling uncertainty, exposure and validity.\n', encoding='utf-8')
    (Path(case)/'results'/'THREE_METHOD_COMPARISON.md').write_text(
        '# DSMC, Cai and virtual-collision comparison\n\n'+metadata['quality']+'\n\n'
        '[Data, diagnostics and all plots](virtual_collision/README.md).\n\n'
        '![Density](virtual_collision/density_three_way.png)\n\n'
        '![Centerline](virtual_collision/centerline_density_three_way.png)\n\n'
        '![Angular flux](virtual_collision/angular_flux_three_way.png)\n', encoding='utf-8')
    return metadata


def run_exports(root, study, options, *, refinement=1, plots=True, selected_kn=None):
    root, study = Path(root).resolve(), Path(study).resolve()
    doc = inventory(root, families=(study.name,))
    cache, processed, skipped = BallisticCache(), [], []
    for row in doc['cases']:
        if row['kind'] != 'case' or row['source_type'] != 'direct-nozzle':
            continue
        if selected_kn is not None and row['nominal_kn'] not in selected_kn:
            continue
        case = root/row['path']
        required = ('meridional_plane.csv', 'centerline_extended.csv', 'angular_r10D.csv')
        missing = [name for name in required if not (case/'results'/name).is_file()]
        if missing:
            skipped.append(dict(path=str(case), kn=row['nominal_kn'], missing=missing,
                                reason='missing existing DSMC/Cai comparison exports'))
            continue
        print(f'Running {case} (Kn={row["nominal_kn"]:g})', flush=True)
        try:
            quality = ('NON-PRODUCTION DSMC REFERENCE' if row['nominal_kn'] in (.001, .0001)
                       else 'EXPLORATORY; DSMC STATIONARITY INCONCLUSIVE')
            report = run_case(case, options, cache, refinement=refinement, plots=plots, quality=quality)
            processed.append(report)
        except (ValueError, OSError, RuntimeError) as exc:
            skipped.append(dict(path=str(case), kn=row['nominal_kn'], reason=str(exc)))
            print(f'  FAILED: {exc}', flush=True)
    output = study/'results'/'virtual_collision'
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'summary.json', dict(source_revision=doc['source_revision'], processed=processed, skipped=skipped))
    write_rows(output/'cases.csv', [dict(path=r['case_path'], kn=r['kn'], quality=r['quality'],
        output=str(Path(r['case_path'])/'results'/'virtual_collision'), wall_s=r['wall_s'],
        core_s=r['core_s'], multiple_event_warning=r['diagnostics']['multiple_event_warning']) for r in processed])
    write_json(output/'missing_data.json', skipped)
    lines = ['# Extension three-method comparisons', '',
        'Existing DSMC and Cai exports plus reduced virtual-collision predictions. '
        'All comparisons remain exploratory; see each case for quality and model-validity flags.', '',
        '| Kn | Case | Comparison |', '|---:|---|---|']
    for report in processed:
        target = Path(report['case_path'])/'results'/'THREE_METHOD_COMPARISON.md'
        import os
        link = Path(os.path.relpath(target, output)).as_posix()
        lines.append(f'| {report["kn"]:g} | {Path(report["case_path"]).parent.name}/{Path(report["case_path"]).name} | [Open plots]({link}) |')
    lines += ['', 'Missing inputs:', '']
    lines += [f'- {record["path"]}: {record.get("reason", "unavailable")} '
              +', '.join(record.get('missing', [])) for record in skipped]
    lines += ['', 'DSMC tensors and independent time-block uncertainty are unavailable in these compact exports.', '']
    (output/'README.md').write_text('\n'.join(lines), encoding='utf-8')
    return processed, skipped
