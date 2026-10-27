#!/usr/bin/env python3
"""Bounded source-only audit of representative-event hazard discretization.

No DSMC fields or spatial transport grids are evaluated. The physical model is
not modified. At finite panel order its event CDF uses nu_mid, whereas partner
selection normalizes by nu_event. R=nu_event/nu_mid measures the discrepancy:
the implemented sampled pair density carries a factor 1/R relative to the
stated symmetric frozen-hazard pair measure. Refinement must make R approach 1
over the event-rate-weighted population before that approximation is negligible.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import sys
import time

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np
from plumetools.virtualcollision import NumericalOptions, VirtualCollisionModel
from plumetools.virtualcollision.export_geometry import load_export_context
from plumetools.virtualcollision.kernel import VHSKernel
from plumetools.virtualcollision.model import ray_exit
from plumetools.virtualcollision.assessment import write_json, source_revision


def weighted_quantiles(values, weights, quantiles):
    order = np.argsort(values)
    values, weights = np.asarray(values)[order], np.asarray(weights)[order]
    positive = weights > 0
    values, weights = values[positive], weights[positive]
    if not len(weights) or weights.sum() <= 0:
        return [None for _ in quantiles]
    # Midpoint cumulative weights avoid assigning a full atom to either edge.
    cumulative = (np.cumsum(weights)-.5*weights)/weights.sum()
    return np.interp(quantiles, cumulative, values).tolist()


def audit_one(case, samples, steps, seed):
    started, cpu = time.perf_counter(), time.process_time()
    context = load_export_context(case)
    source = context['source']
    options = NumericalOptions(samples=samples, batches=8, hazard_steps=steps,
        partner_radial=4, partner_azimuthal=16, partner_speed=16, seed=seed,
        domain_length_D=20., domain_radius_D=16., chunk_size=128)
    model = VirtualCollisionModel(source, VHSKernel(source.species), options)
    diagnostics = model.preprocess()
    ratios, event_rates, locations, panel_indices = [], [], [], []
    per_batch = samples//options.batches
    for batch_id, rays in enumerate(model.events):
        positions, velocities = source.sample_flux(per_batch, seed+batch_id*104729)
        if len(rays) != 4*per_batch:
            raise RuntimeError('Audit requires one positive-probability event proposal per sampled source ray')
        for i, (position, velocity) in enumerate(zip(positions, velocities)):
            location, _, event_rate = rays[4*i]
            duration, _ = ray_exit(position, velocity, model.length, model.radius)
            increments, edges = model._hazard(position, velocity, duration)
            event_time = float((location[0]-position[0])/velocity[0])
            panel = int(np.searchsorted(edges, event_time, side='right')-1)
            panel = max(0, min(panel, steps-1))
            nu_mid = increments[panel]/(edges[panel+1]-edges[panel])
            nu_event = model._frequency(location, velocity)[0]
            if not nu_mid > 0 or not np.isfinite(nu_event):
                raise FloatingPointError('Nonpositive proposal hazard or failed event frequency')
            ratios.append(nu_event/nu_mid)
            event_rates.append(event_rate)
            locations.append(location.tolist())
            panel_indices.append(panel)
    ratios, weights = np.array(ratios), np.array(event_rates)
    total = weights.sum()
    quantiles = [.01, .05, .25, .5, .75, .95, .99]
    return dict(case=str(case.relative_to(ROOT)), kn=float(context['cfg'].exit.knudsen),
        samples=samples, batches=options.batches, seed=seed, hazard_steps=steps,
        partner_quadrature=[4, 16, 16], domain_length_D=20., domain_radius_D=16.,
        wall_s=time.perf_counter()-started, core_s=time.process_time()-cpu,
        exposure_mean=diagnostics['exposure_mean'],
        ensemble_event_rate_over_emission=diagnostics['ensemble_event_rate_over_emission'],
        remaining_exposure_weighted_mean=diagnostics['remaining_exposure_weighted_mean'],
        event_weighted_mean_nu_event_over_nu_mid=float(weights@ratios/total),
        event_weighted_mean_nu_mid_over_nu_event=float(weights@(1/ratios)/total),
        event_weight_fraction_ratio_outside_0p9_1p1=float(weights[(ratios < .9)|(ratios > 1.1)].sum()/total),
        ratio_min=float(ratios.min()), ratio_max=float(ratios.max()),
        weighted_quantiles=dict(zip(map(str, quantiles), weighted_quantiles(ratios, weights, quantiles))),
        unweighted_quantiles=dict(zip(map(str, quantiles), np.quantile(ratios, quantiles).tolist())),
        nu_event_over_nu_mid=ratios.tolist(), normalized_event_weights=(weights/total).tolist(),
        event_locations_m=locations, panel_indices=panel_indices,
        interpretation='Diagnostic only. Source samples and event locations change with panel refinement; not a field-accuracy or DSMC validation result.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=128)
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--steps', type=int, nargs='+', default=[24, 48, 96])
    parser.add_argument('--kn', type=float, nargs='+', choices=(1., 10.), default=[10., 1.])
    parser.add_argument('--output', type=Path, default=ROOT/'cases/cai2012-kn-extension/results/virtual_collision/review/hazard_audit.json')
    args = parser.parse_args()
    if not 8 <= args.samples <= 128 or args.samples % 8:
        parser.error('This bounded audit requires samples divisible by 8, in [8,128]')
    if any(steps < 2 or steps > 96 for steps in args.steps) or len(args.steps)*len(args.kn) > 6:
        parser.error('This bounded audit permits at most six runs and hazard orders in [2,96]')
    rows = []
    for kn in args.kn:
        case = ROOT/f'cases/cai2012-kn-extension/Cases_run1_graded/Kn{kn:g}_np1x'
        for steps in args.steps:
            print(f'Hazard audit Kn={kn:g}, samples={args.samples}, panels={steps}', flush=True)
            row = audit_one(case, args.samples, steps, args.seed)
            rows.append(row)
            print(f'  R mean={row["event_weighted_mean_nu_event_over_nu_mid"]:.6g}, '
                  f'weight outside +/-10%={row["event_weight_fraction_ratio_outside_0p9_1p1"]:.3%}, '
                  f'wall={row["wall_s"]:.2f}s', flush=True)
            write_json(args.output, dict(source_revision=source_revision(ROOT), platform=platform.platform(),
                scope='source-only diagnostic; no DSMC, no field transport, no model changes',
                definition='R=nu_event/nu_mid; implemented pair measure differs by factor 1/R at finite panels',
                thread_environment={name: os.environ[name] for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')},
                runs=rows))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
