"""Small fixed-count patches; domain changes are explicit local experiments."""
from dataclasses import replace
import json
import math
from pathlib import Path
import re

import numpy as np

from .config import digest
from .physics import source_for_kn,_physical_hash,predict_breakdown
from .mesh import initial_mesh


def local_identity(model):
    cfg = model.config
    physical = {**cfg.physical_inputs,'exit':{**cfg.physical_inputs['exit'],'knudsen':model.source.knudsen}}
    return {'kind':'fixed_count_local','sha256':digest({
        'physical':physical,'reference_sha256':cfg.reference_sha256,'inlet':{k:float(v) for k,v in model.state(model.r_in).items()},
        'r_in_m':model.r_in,'cells':100,'targets':vars(cfg.targets),'local':vars(cfg.local),
        'time':vars(cfg.time),'controller':vars(cfg.controller),'run':vars(cfg.run)})}


def _patch(model,span,angle):
    if not math.isfinite(span+angle) or span<=64*np.spacing(model.r_in) or angle<=0:
        raise ValueError('local patch dimensions are below representable geometry precision')
    cfg = replace(model.config,domain=replace(model.config.domain,
        inlet_over_D=model.r_in/model.source.diameter_m,
        outer_over_D=(model.r_in+span)/model.source.diameter_m,angle_deg=math.degrees(angle)))
    return replace(model,config=cfg,r_out=model.r_in+span,sampling_record=None,
                   inlet_state={k:float(v) for k,v in model.state(model.r_in).items()},
                   physical_hash=_physical_hash(cfg,model.source.knudsen,model.r_in,model.r_in+span))


def contract_local_model(model,mesh,target):
    """Fit100 actual cells inside a trusted uniform upper edge bound."""
    if not math.isfinite(target) or target<=0:
        raise ValueError('local contraction requires a finite positive measured target')
    target *= .95  # Small geometric margin, never a further resolution sweep.
    span = min(model.r_out-model.r_in,100*target/math.sqrt(1+2*math.tan(mesh.angle_rad)**2))
    angle = min(mesh.angle_rad,math.atan(target/(2*(model.r_in+span))))
    child = _patch(model,span,angle)
    result = initial_mesh(child)
    if result.n_cells!=100 or result.transverse!=1 or np.max(result.geometry.max_edge_m)>target*(1+1e-10):
        raise ValueError('local patch failed its fixed-count edge construction')
    return child,result


def initial_local_model(config,kn):
    model = source_for_kn(config,kn)
    span = min(config.local.max_span_over_D*model.source.diameter_m,model.r_out-model.r_in)
    model = _patch(model,span,math.radians(config.domain.angle_deg))
    mesh = initial_mesh(model)
    lam = np.asarray(model.state(np.linspace(model.r_in,model.r_out,101))['lambda'])
    if not np.all(np.isfinite(lam)&(lam>0)):
        raise ValueError('local initial patch has no valid predicted MFP coverage')
    model,_ = contract_local_model(model,mesh,config.targets.edge_over_mfp*float(np.min(lam)))
    return model


def initial_weight(model,mesh,target_parcels):
    target = max(target_parcels,model.config.local.initial_parcels)
    return float(np.min(mesh.density_estimate*mesh.geometry.volumes_m3)/target)


def reference_breakdown(model):
    outer = model.config.local.reference_outer_over_D*model.source.diameter_m
    result = predict_breakdown(replace(model,r_out=outer),**vars(model.config.breakdown))
    return {**result,'scope':'analytical reference, not measured by the local patch',
            'reference_outer_over_D':model.config.local.reference_outer_over_D,
            'outside_patch':bool(result.get('r_m',0)>model.r_out) if 'r_m' in result else None}


def revise_local_case(parent,plan,*,correction='needed',delta_t_s=None,dry_run=False):
    from .case import load_case,write_case,sizing
    from plumetools.refinement.plan import verify_sources,write_plan,plan_identity
    parent = Path(parent).resolve()
    verify_sources(parent,plan)
    model,mesh,old = load_case(parent)
    if not model.config.local.enabled or plan.action=='unresolved':
        raise ValueError('local contraction requires an enabled local case and trusted measurements')
    weight,dt = old['numerics']['particle_weight'],old['numerics']['delta_t_s']
    changed = []
    if correction in ('needed','mesh') and plan.action=='needs_mesh':
        if (len(plan.target_edge_m)!=100 or not np.all(np.isfinite(plan.target_edge_m)&(plan.target_edge_m>0))
                or not np.all(np.isfinite(plan.density_estimate)&(plan.density_estimate>0))):
            raise ValueError('local contraction requires complete measured coverage')
        old_centres = mesh.geometry.centres_m[:,0]
        model,mesh = contract_local_model(model,mesh,float(np.min(plan.target_edge_m)))
        density = np.interp(mesh.geometry.centres_m[:,0],old_centres,plan.density_estimate)
        mesh = replace(mesh,layer_density=density,density_kind='parent-measured')
        weight = initial_weight(model,mesh,plan.target_parcels)
        dt = min(dt,sizing(model,mesh,plan.target_parcels,weight)[0]['delta_t_s'])
        changed.append('domain')
    if correction in ('needed','parcels') and plan.recommended_weight is not None and plan.recommended_weight<weight:
        weight = plan.recommended_weight
        changed.append('parcels')
    if correction in ('needed','timestep') and delta_t_s is not None and delta_t_s<dt:
        dt = delta_t_s
        changed.append('timestep')
    if not changed:
        return {'action':'no_revision_needed','reason':'no selected local correction is required'}
    revisions = [int(p.name.split('_')[0][1:]) for p in parent.parent.glob('r*_ppc*')
                 if re.fullmatch(r'r\d+_ppc[0-9pem+]+',p.name) and (p/'case.json').is_file()]
    from .generation import label
    child = parent.parent/f'r{max(revisions)+1}_ppc{label(plan.target_parcels)}'
    num,sampling = sizing(model,mesh,plan.target_parcels,weight,dt)
    row = {'action':'create_revision','path':str(child),'corrections':changed,'numerics':num,
           'sampling':sampling,'physical_identity':mesh.identity,
           'parent_plan_id':plan_identity(plan),'requires_flow_rerun':True}
    if not dry_run:
        def finalize(stage):
            write_plan(plan,stage/'parent-plan')
            (stage/'revision.json').write_text(json.dumps({k:v for k,v in row.items() if k!='record'},indent=2)+'\n')
        row['record'] = write_case(child,model,mesh,target_parcels=plan.target_parcels,
                                   particle_weight=weight,delta_t_s=dt,parent=parent,_finalize=finalize)
    return row
