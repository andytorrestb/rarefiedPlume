"""Small fixed-count patches; domain changes are explicit local experiments."""
from dataclasses import replace
import math

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
    angle = min(mesh.angle_rad,math.atan(target/(2*model.r_out)))
    span = min(model.r_out-model.r_in,100*target/math.sqrt(1+2*math.tan(angle)**2))
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
