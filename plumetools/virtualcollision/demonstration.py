"""Bounded predictive-only convergence, conservation and cost demonstration.

No DSMC data enter these calculations. This is not a scientific accuracy study.
"""
from __future__ import annotations
from dataclasses import asdict, replace
from pathlib import Path
import os
import platform
import time
import numpy as np

from . import DiskSource, NumericalOptions, VirtualCollisionModel
from .kernel import VHSKernel
from .assessment import load_yaml, write_json, write_rows
from plumetools.cai2012.gas import VhsSpecies, number_density_from_kn


def _peak_memory():
    try:
        import resource
        raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return raw if platform.system() == 'Darwin' else raw*1024
    except (ImportError, AttributeError):
        return None


def _save(path, result, x_bounds, r_bounds):
    arrays = {key: result[key] for key in ('M0', 'M1', 'M2', 'n', 'U', 'T_tensor',
                'T_trans', 'P_tensor', 'valid', 'tensor_eigenvalues_K')}
    arrays.update({'se_'+key: val for key, val in result['standard_error'].items()})
    np.savez_compressed(path, x_bounds_m=x_bounds, r_bounds_m=r_bounds, **arrays)


def _relative_change(ref, other):
    return float(np.linalg.norm(other-ref)/max(np.linalg.norm(ref), 1e-300))


def run_model_demonstration(cfg, output, smoke=True, plots=True):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    root = Path(cfg.get('_root', Path(__file__).resolve().parents[2]))
    path = root/cfg.get('source_case', 'cases/cai2012-kn-extension/baseCase')/'case.yaml'
    physical = load_yaml(path)
    if physical.get('model') != 'cai2012_circular_plume':
        raise ValueError('demonstration supports only configured Cai uniform disk sources')
    gas, exit_cfg = physical['gas'], physical['exit']
    species = VhsSpecies(name=gas['species_name'], **{key: gas[key] for key in
        ('mass_kg', 'diameter_m', 'omega', 't_ref_K', 'internal_degrees_of_freedom')})
    law = physical['dsmc']['binary_collision_model']
    kernel = VHSKernel(species, law=law)
    numerical = dict(cfg['smoke' if smoke else 'full'])
    kns = numerical.pop('kn', cfg.get('kn', [10.]))
    opts = NumericalOptions(**numerical)
    diameter = float(physical['nozzle']['diameter_m'])
    convention = gas['mean_free_path_convention']
    if exit_cfg['characteristic_length'] != 'diameter':
        raise ValueError('unsupported Kn characteristic length')
    x = np.array([.5, 1.5, 4.]) if smoke else np.linspace(.5, 10., 24)
    xb = diameter*np.column_stack((x-.2, x+.2))
    rb = np.tile([0., .25*diameter], (len(x), 1))
    costs, rows, reports, refinement_rows = [], [], [], []
    for kn in kns:
        density = number_density_from_kn(float(kn), diameter, species,
            float(exit_cfg['T0_K']), convention=convention)
        source = DiskSource(diameter/2, density, float(exit_cfg['T0_K']),
                            float(exit_cfg['speed_ratio']), species)
        model = VirtualCollisionModel(source, kernel, opts)
        cpu = time.process_time()
        start = time.perf_counter()
        cai = model.evaluate_bins(xb, rb, collisions=False)
        cai_time, cai_cpu = time.perf_counter()-start, time.process_time()-cpu
        cpu = time.process_time()
        start = time.perf_counter()
        diagnostics = model.preprocess()
        preprocessing, preprocessing_cpu = time.perf_counter()-start, time.process_time()-cpu
        cpu = time.process_time()
        start = time.perf_counter()
        vc = model.evaluate_bins(xb, rb)
        evaluation, evaluation_cpu = time.perf_counter()-start, time.process_time()-cpu
        # The preprocessing cache is reusable: its key includes every source,
        # collision-law and numerical control. No persisted stale cache is read.
        model.preprocess()
        tag = f'Kn{float(kn):g}'.replace('.', 'p')
        _save(output/(tag+'_cai.npz'), cai, xb, rb)
        _save(output/(tag+'_virtual.npz'), vc, xb, rb)
        write_json(output/(tag+'_diagnostics.json'), vc['diagnostics'])
        cost = dict(kn=float(kn), bins=len(x), requested_outputs='M0,M1,M2,n,U,T_tensor,T_trans,P_tensor',
            cai_evaluation_wall_s=cai_time, cai_core_s=cai_cpu,
            virtual_preprocessing_wall_s=preprocessing, virtual_preprocessing_core_s=preprocessing_cpu,
            virtual_evaluation_wall_s=evaluation, virtual_evaluation_core_s=evaluation_cpu,
            virtual_total_over_cai=(preprocessing+evaluation)/cai_time,
            peak_process_rss_bytes=_peak_memory(),
            peak_memory_scope='process lifetime high-water mark, not incremental allocation',
            cold_preprocessing_cache=True, warm_evaluation_cache=True,
            samples=opts.samples, options=asdict(opts),
            dsmc_production_core_s=None, dsmc_cost_reason='credible production resource accounting unavailable')
        costs.append(cost)
        for i, location in enumerate(x):
            rows.append(dict(kn=float(kn), x_over_D=float(location), estimator='volume-averaged cylindrical ring',
                r_low_over_D=0., r_high_over_D=.25, x_half_width_over_D=.2,
                cai_n_over_n0=cai['n'][i]/density, virtual_n_over_n0=vc['n'][i]/density,
                virtual_n_se_over_n0=vc['standard_error']['n'][i]/density,
                cai_Ux_m_s=cai['U'][i, 0], virtual_Ux_m_s=vc['U'][i, 0],
                cai_T_K=cai['T_trans'][i], virtual_T_K=vc['T_trans'][i],
                numerical_realizability=bool(vc['valid'][i]),
                single_event_warning=diagnostics['multiple_event_warning'],
                dsmc_reference='unavailable in predictive demonstration'))
        reports.append(dict(kn=float(kn), source=asdict(source), diagnostics=diagnostics,
                            valid_bins=int(vc['valid'].sum()), bins=len(x)))
        if plots:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
                for ax, key, scale, label in zip(axes, ('n', 'U', 'T_trans'),
                    (density, source.thermal_speed, source.temperature_K),
                    ('n/n0', 'Ux sqrt(beta0)', 'T/T0')):
                    baseline = cai[key][:, 0] if key == 'U' else cai[key]
                    corrected = vc[key][:, 0] if key == 'U' else vc[key]
                    se = vc['standard_error'][key][:, 0] if key == 'U' else vc['standard_error'][key]
                    ax.plot(x, baseline/scale, 'k-', label='Cai collisionless')
                    # Raw invalid values remain visible; a cross marks their status.
                    ax.errorbar(x, corrected/scale, yerr=se/scale, fmt='o-',
                                label='Virtual-collision model', color='tab:orange')
                    ax.plot(x[~vc['valid']], corrected[~vc['valid']]/scale, 'rx', ms=9)
                    ax.set(xlabel='X/D', ylabel=label)
                    ax.grid(alpha=.2)
                axes[0].legend(fontsize=7)
                fig.suptitle(f'Kn={float(kn):g}; same ring bins; predictive demonstration, no DSMC score', fontsize=9)
                fig.tight_layout()
                fig.savefig(output/(tag+'_model_only.png'), dpi=160)
                plt.close(fig)
            except ImportError:
                pass
        # A deliberately small, one-control-at-a-time assessment. Keep identical
        # seeds where possible. These changes measure uncertainty, not acceptance
        # of convergence; no physical parameters are optimized against DSMC.
        if float(kn) == float(kns[0]):
            variants = dict(samples=replace(opts, samples=2*opts.samples),
                hazard=replace(opts, hazard_steps=2*opts.hazard_steps),
                partner_quadrature=replace(opts, partner_radial=2*opts.partner_radial,
                    partner_azimuthal=2*opts.partner_azimuthal, partner_speed=2*opts.partner_speed),
                spatial_quadrature=replace(opts, bin_order=2*opts.bin_order,
                    residence_order=2*opts.residence_order),
                domain=replace(opts, domain_length_D=2*opts.domain_length_D,
                    domain_radius_D=2*opts.domain_radius_D))
            for name, refined in variants.items():
                start = time.perf_counter()
                alternative = VirtualCollisionModel(source, kernel, refined)
                alternative.preprocess()
                result = alternative.evaluate_bins(xb, rb)
                _save(output/(tag+'_refinement_'+name+'.npz'), result, xb, rb)
                write_json(output/(tag+'_refinement_'+name+'_diagnostics.json'), result['diagnostics'])
                refinement_rows.append(dict(kn=float(kn), control=name, options=asdict(refined),
                    density_relative_change=_relative_change(vc['M0'], result['M0']),
                    momentum_relative_change=_relative_change(vc['M1'], result['M1']),
                    second_moment_relative_change=_relative_change(vc['M2'], result['M2']),
                    density_se_rms_over_n0=float(np.sqrt(np.mean(result['standard_error']['n']**2))/density),
                    exposure_mean=result['diagnostics']['exposure_mean'],
                    wall_s=time.perf_counter()-start,
                    conclusion='measured change only; no convergence acceptance inferred'))
    write_rows(output/'model_profiles.csv', rows)
    write_rows(output/'cost.csv', costs)
    write_rows(output/'convergence.csv', refinement_rows)
    write_json(output/'convergence.json', refinement_rows)
    return dict(mode='smoke' if smoke else 'full', source_case=str(path),
        scope='predictive bounded demonstration; no DSMC input or fitted parameter',
        numerical_options=asdict(opts), hardware=dict(platform=platform.platform(),
            processor=platform.processor(), logical_cpus=os.cpu_count(),
            thread_environment={key: os.environ.get(key) for key in
                ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}),
        costs=costs, cases=reports, refinements=refinement_rows,
        uncertainty='batch SE is sampling uncertainty only; DSMC uncertainty unavailable',
        conclusion='Scientific effectiveness and convergence have not been established by this smoke demonstration.')
