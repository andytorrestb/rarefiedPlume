"""Seven initial cases; append-only measured numerical revisions."""
from dataclasses import replace
import json
from pathlib import Path
import re

from plumetools.refinement.plan import verify_sources, write_plan, plan_identity
from .case import write_case, load_case, sizing
from .mesh import initial_mesh, refine_mesh
from .physics import source_for_kn, predict_breakdown


def label(x):
    return f'{float(x):g}'.replace('.','p').replace('-','m')


def generate_initial(config, output, *, kn_values=None, dry_run=False):
    values = config.knudsen_numbers if kn_values is None else tuple(kn_values)
    if any(k not in config.knudsen_numbers for k in values):
        raise ValueError('selected Kn is absent from the study')
    rows = []
    for kn in values:
        model = source_for_kn(config,kn)
        mesh = initial_mesh(model)
        p = Path(output)/f'Kn{label(kn)}'/f'r0_ppc{label(config.targets.min_mean_parcels)}'
        num, sampling = sizing(model,mesh,config.targets.min_mean_parcels)
        row = {'path':str(p), 'Kn':kn,'numerics':num, 'sampling':sampling,
               'breakdown':predict_breakdown(model,**vars(config.breakdown)),
               'physical_identity':mesh.identity}
        if not dry_run:
            row['record'] = write_case(p,model,mesh,target_parcels=config.targets.min_mean_parcels)
        rows.append(row)
    return rows


def revise_case(parent, plan, *, correction='needed', delta_t_s=None, dry_run=False):
    if correction not in ('needed','mesh','parcels','timestep'):
        raise ValueError('correction must be needed, mesh, parcels or timestep')
    if plan.action == 'criteria_satisfied' and delta_t_s is None:
        return {'action':'no_revision_needed','reason':'selected measured criteria satisfied'}
    if plan.action == 'unresolved':
        raise ValueError('measurements are unresolved; obtain adequate sampling before revision')
    parent = Path(parent).resolve()
    verify_sources(parent,plan)
    model, mesh, old = load_case(parent)
    weight, dt = old['numerics']['particle_weight'], old['numerics']['delta_t_s']
    changed = []
    if correction in ('needed','mesh') and plan.action == 'needs_mesh':
        mesh = refine_mesh(mesh,plan)
        weight = None
        changed.append('mesh')
    if correction in ('needed','parcels') and plan.recommended_weight is not None and plan.recommended_weight < old['numerics']['particle_weight']:
        if weight is not None:
            weight = plan.recommended_weight
        changed.append('parcels')
    if correction in ('needed','timestep') and delta_t_s is not None and delta_t_s < dt:
        dt = delta_t_s
        changed.append('timestep')
    if not changed:
        return {'action':'no_revision_needed','reason':'no selected correction is required'}
    if 'mesh' in changed:
        dt = min(dt,sizing(model,mesh,plan.target_parcels)[0]['delta_t_s'])
    revision = 1+max(int(p.name.split('_')[0][1:]) for p in parent.parent.glob('r*_ppc*')
                    if re.fullmatch(r'r\d+_ppc[0-9pem+]+',p.name) and (p/'case.json').is_file())
    child = parent.parent/f'r{revision}_ppc{label(plan.target_parcels)}'
    row = {'action':'create_revision','path':str(child),'corrections':changed,
           'parent_plan_id':plan_identity(plan),'requires_flow_rerun':True,
           'numerics':sizing(model,mesh,plan.target_parcels,weight,dt)[0]}
    if not dry_run:
        def finalize(stage):
            write_plan(plan,stage/'parent-plan')
            (stage/'revision.json').write_text(json.dumps({k:v for k,v in row.items() if k not in ('path','record')},indent=2)+'\n')
        row['record'] = write_case(child,model,mesh,target_parcels=plan.target_parcels,
                                   particle_weight=weight,delta_t_s=dt,parent=parent,_finalize=finalize)
    return row
