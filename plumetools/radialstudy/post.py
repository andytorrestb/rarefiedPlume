"""Stream measured windows through moments, then assess actual cells and profiles."""
from __future__ import annotations

from dataclasses import asdict
import itertools
import json
import math
from pathlib import Path

import numpy as np

from plumetools.rarefaction.core import Context,translational_temperature,mean_free_path
from plumetools.rarefaction.settings import DiagnosticSettings
from plumetools.cai2012ext.overlap import normalised_l2,verdict
from plumetools.refinement.assess import assess
from plumetools.refinement.geometry import read_geometry,field_values,file_sha256
from plumetools.refinement.model import Measurements
from plumetools.refinement.plan import plan_refinement,write_plan,json_safe,plan_identity
from plumetools.refinement.report import write_assessment
from .case import load_case
from .config import digest
from .run import read_windows,collect_case


def serial(value):
    if isinstance(value,np.ndarray):
        return serial(value.tolist())
    if isinstance(value,dict):
        return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [serial(v) for v in value]
    return json_safe(value)


def primitive(raw,species):
    rho = np.asarray(raw['rho'])
    momentum = np.asarray(raw['momentum'])
    energy = np.asarray(raw['energy'])
    with np.errstate(divide='ignore',invalid='ignore'):
        u = momentum/rho[:,None]
    return {'n':rho/species.mass_kg,'U':u,'T':translational_temperature(rho,momentum,energy,species.gas_constant_j_per_kg_k)}


def project_moments(source_bounds,raw,target_bounds):
    """Conservative frustum-overlap projection; preserve density/momentum/energy."""
    source_bounds,target_bounds = np.asarray(source_bounds),np.asarray(target_bounds)
    if np.any(np.diff(source_bounds)<=0) or np.any(np.diff(target_bounds)<=0):
        raise ValueError('projection boundaries must increase')
    result = {k:np.zeros((len(target_bounds)-1,)+np.asarray(v).shape[1:]) for k,v in raw.items()}
    coverage = np.zeros(len(target_bounds)-1)
    for j,(a,b) in enumerate(zip(target_bounds[:-1],target_bounds[1:])):
        first = max(0,np.searchsorted(source_bounds,a,side='right')-1)
        last = min(len(source_bounds)-1,np.searchsorted(source_bounds,b,side='left'))
        denom = b**3-a**3
        for i in range(first,last):
            lo,hi = max(a,source_bounds[i]),min(b,source_bounds[i+1])
            if hi>lo:
                w = (hi**3-lo**3)/denom
                coverage[j] += w
                for k,v in raw.items():
                    result[k][j] += np.asarray(v)[i]*w
        if coverage[j] < 1-1e-10:
            for value in result.values():
                value[j] = np.nan
    return result,coverage


def _read_window(case,w,n,species):
    roles = w['field_roles']
    folder = case/w['time']
    density = field_values(folder/roles['number_density'],n,[0,-3,0,0,0,0,0])
    temperature = field_values(folder/roles['temperature'],n,[0,0,0,1,0,0,0])
    velocity = field_values(folder/roles['velocity'],n,[0,1,-1,0,0,0,0],vector=True)
    parcels = field_values(folder/roles['mean_parcels'],n,[0,0,0,0,0,0,0])
    rho = density*species.mass_kg
    raw = {'rho':rho,'momentum':rho[:,None]*velocity,
           'energy':.5*rho*(3*species.gas_constant_j_per_kg_k*temperature+np.sum(velocity**2,axis=1)),
           'parcels':parcels}
    # Empty/sentinel cells remain unresolved, even if other windows populate them.
    raw['bad'] = ((~np.isfinite(density)) | (density<=0) | (~np.isfinite(temperature))
                  | (temperature<=0) | (~np.all(np.isfinite(velocity),axis=1))
                  | (~np.isfinite(parcels)) | (parcels<=0)).astype(float)
    for stem,dims in (('mfp',[0,1,0,0,0,0,0]),('mct',[0,0,1,0,0,0,0])):
        path = folder/f'{stem}_{species.name}'
        if path.is_file():
            v = field_values(path,n,dims)
            raw[stem] = np.where(np.isfinite(v)&(v>0)&(v<1e20),v,np.nan)
    return raw


def _profile(raw,mesh,geometry):
    nt = mesh.transverse**2
    v = geometry.volumes_m3.reshape(-1,nt)
    result = {}
    for k in ('rho','momentum','energy'):
        data = np.asarray(raw[k])
        reshaped = data.reshape((len(v),nt)+data.shape[1:])
        weights = v if data.ndim==1 else v[:,:,None]
        denom = v.sum(axis=1) if data.ndim==1 else v.sum(axis=1)[:,None]
        result[k] = np.sum(reshaped*weights,axis=1)/denom
    return result


def _statistics(blocks,species,volumes,expected):
    if len(blocks) != expected:
        return {'complete':False,'stationary':False,'independence_supported':False,'reason':'missing complete blocks'},{}
    states = [primitive(b,species) for b in blocks]
    se,drift,correlated = {},False,False
    for key in ('n','Ux','T'):
        data = np.stack([st['U'][:,0] if key=='Ux' else st[key] for st in states])
        se[key] = np.std(data,axis=0,ddof=1)/math.sqrt(len(data))
        if not np.all(np.isfinite(data)):
            return {'complete':True,'stationary':False,'independence_supported':False,'reason':'invalid block profile'},se
        mid = len(data)//2
        left,right = data[:mid],data[mid:]
        diff = np.abs(left.mean(axis=0)-right.mean(axis=0))
        error = np.sqrt(np.var(left,axis=0,ddof=1)/len(left)+np.var(right,axis=0,ddof=1)/len(right)) if min(len(left),len(right))>1 else np.full(data.shape[1],np.inf)
        scale = np.maximum(np.abs(data.mean(axis=0)),1e-30)
        drift |= bool(np.any(diff > np.maximum(3*error,.01*scale)))
        # A block series may look steady yet retain strong temporal correlation.
        series = data@volumes/np.sum(volumes)
        if len(series)>=6 and np.ptp(series)>1e-12*max(abs(np.mean(series)),1e-30):
            if np.std(series[:-1])>0 and np.std(series[1:])>0:
                correlated |= abs(float(np.corrcoef(series[:-1],series[1:])[0,1])) > .8
    return {'complete':True,'stationary':not drift and not correlated,
            'independence_supported':not correlated,'drift_detected':drift,
            'strong_lag_one_correlation':correlated,'blocks':len(blocks),
            'test':'half-series 3SE with 1% floor; |lag1| <= 0.8; diagnostic, not proof'},se


def measured_breakdown(r,n,u,lam,cbar,valid,*,threshold=.05,density_se=None):
    """General signed Bird diagnostic using a numerical radial density gradient."""
    r,n,u,lam,cbar = map(np.asarray,(r,n,u,lam,cbar))
    if len(r)<3 or not np.all(valid) or not np.all(np.isfinite(n)&(n>0)):
        return {'status':'UNAVAILABLE','reason':'insufficient valid radial gradient coverage'}
    # numpy's nonuniform three-point stencil is appropriate for this 1-D profile;
    # the shared 3-D least-squares fit intentionally rejects its rank-one stencil.
    gradient = np.gradient(np.log(n),r,edge_order=2)
    b = -lam*u/cbar*gradient
    result = {'basis':'measured_radial_gradient','definition':'general','B':b.tolist(),'threshold':threshold}
    if not np.all(np.isfinite(b)):
        return {**result,'status':'UNAVAILABLE','reason':'nonfinite gradient'}
    if density_se is not None:
        noise = lam*np.abs(u)/cbar*np.asarray(density_se)/n/np.gradient(r)
        if np.any(~np.isfinite(noise)) or np.any(noise>.5*threshold):
            return {**result,'status':'UNAVAILABLE','reason':'gradient noise exceeds half the threshold'}
    if b[0]>=threshold:
        return {**result,'status':'ALREADY_NONCONTINUUM_AT_INLET'}
    crosses = np.flatnonzero((b[:-1]<threshold)&(b[1:]>=threshold))
    if not len(crosses):
        return {**result,'status':'NOT_REACHED_IN_DOMAIN'}
    i = int(crosses[0])
    return {**result,'status':'CROSSING_IN_DOMAIN','r_m':float(np.interp(threshold,b[i:i+2],r[i:i+2]))}


def _sizing_statistics(windows,species,settings):
    result = {'basis':'adaptive_sizing','complete':False,'stationary':False,
              'independence_supported':False,'relative_tolerance':settings.stability_relative,
              'floors':{'n':settings.density_floor,'U':settings.velocity_floor,'T':settings.temperature_floor},
              'reason':'need more complete post-startup windows'}
    if len(windows)<settings.min_windows:
        return result,{}
    if any(np.any(w['bad']) for w in windows):
        return {**result,'reason':'missing or invalid per-cell measurements'},{}
    mid = len(windows)//2
    states = [primitive({key:np.mean([w[key] for w in group],axis=0)
                         for key in ('rho','momentum','energy')},species)
              for group in (windows[:mid],windows[mid:])]
    changes = {}
    for key,floor in (('n',settings.density_floor),('U',settings.velocity_floor),('T',settings.temperature_floor)):
        a,b = (st[key] for st in states)
        magnitude = (lambda x:np.linalg.norm(x,axis=1)) if key=='U' else np.abs
        ratio = magnitude(a-b)/np.maximum(np.maximum(magnitude(a),magnitude(b)),floor)
        changes[key] = float(np.max(ratio))
    stable = all(math.isfinite(x) and x<=settings.stability_relative for x in changes.values())
    return {**result,'stationary':stable,'blocks':2,'max_relative_change':changes,
            'reason':'stable sizing measurements' if stable else 'per-cell block means still changing'},{}


def _compute(case, *, incremental=None):
    case = Path(case).resolve()
    model,mesh,record = load_case(case)
    if incremental is None:
        incremental = model.config.controller.enabled
    windows = read_windows(case)
    if not windows:
        raise ValueError('no measured post-startup windows')
    observed_windows = len(windows)
    if incremental:
        windows = windows[-model.config.controller.min_windows:]
    geometry = read_geometry(case,time='0')
    if len(geometry.volumes_m3) != mesh.n_cells:
        raise ValueError('actual mesh differs from radial metadata')
    source_files = {p:h for w in windows for p,h in w['source_files'].items()}
    source_files.update(record['input_files'])
    source_files['case.json'] = file_sha256(case/'case.json')
    for relative in ('0/C','0/V','run-state.json','sampling-windows.json'):
        source_files[relative] = file_sha256(case/relative)
    total,blocks,weight_total = None,[],0.
    recent = []
    complete = True
    for block,group in itertools.groupby(windows,key=lambda w:w['block']):
        accum,steps,count = None,0.,0
        for w in group:
            raw = _read_window(case,w,mesh.n_cells,model.source.species)
            if incremental:
                recent.append(raw)
            if accum is None:
                accum = {k:np.zeros_like(v) for k,v in raw.items()}
            if set(raw) != set(accum):
                raise ValueError('inconsistent native cross-check fields across windows')
            for k,v in raw.items():
                accum[k] += v*w['steps']
            steps += w['steps']; count += 1
        if total is None:
            total = {k:np.zeros_like(v) for k,v in accum.items()}
        for k,v in accum.items():
            total[k] += v
        weight_total += steps
        if count == record['sampling']['windows_per_block']:
            blocks.append(_profile({k:v/steps for k,v in accum.items()},mesh,geometry))
        else:
            complete = False
    raw = {k:v/weight_total for k,v in total.items()}
    states = primitive(raw,model.source.species)
    bad = raw['bad']>0
    states['n'][bad] = np.nan; states['T'][bad] = np.nan
    profile_raw = _profile(raw,mesh,geometry)
    profile = primitive(profile_raw,model.source.species)
    volume = geometry.volumes_m3.reshape(-1,mesh.transverse**2).sum(axis=1)
    sampling,se = (_sizing_statistics(recent,model.source.species,model.config.controller) if incremental else
                   _statistics(blocks,model.source.species,volume,record['sampling']['blocks']))
    expected = record['sampling']['blocks']*record['sampling']['windows_per_block']
    if not incremental:
        sampling['complete'] &= complete and len(windows)==expected
        sampling['stationary'] &= sampling['complete']
    sampling['windows'] = len(windows)
    sampling['observed_windows'] = observed_windows
    window = {'start_s':windows[0]['start_s'],'end_s':windows[-1]['end_s'],
              'basis':'mean_per_cell','source':'verified','stationary':sampling['stationary'],
              'steps':int(weight_total),'native_basis':'disjoint_native_windows'}
    provenance = {'source_files':source_files,'windows':windows,'sampling':sampling,
                  'particle_weight':record['numerics']['particle_weight'],'case_id':record['case_id']}
    measurement = Measurements(states['n'],states['T'],raw['parcels'],window,digest(source_files),provenance)
    assessed = assess(geometry,measurement,model.source.species,convention=model.source.convention,targets=model.config.targets)
    plan = plan_refinement(geometry,measurement,assessed,particle_weight=record['numerics']['particle_weight'])
    ctx = Context({'number_density':states['n'],'temperature':states['T'],'parcels_per_step':raw['parcels']},
                  species=model.source.species,settings=DiagnosticSettings(),convention=model.source.convention)
    cbar = np.sqrt(8*model.source.gas_constant*states['T']/math.pi)
    tau = assessed.lambda_m/cbar
    speed = np.linalg.norm(states['U'],axis=1)+3*np.sqrt(model.source.gas_constant*states['T'])
    dt = record['numerics']['delta_t_s']
    collision_ratio = dt/tau
    transport_ratio = dt*speed/mesh.transport_length_m
    physical = np.isfinite(collision_ratio)&np.isfinite(transport_ratio)&(tau>0)&(speed>0)
    failures = physical & ((collision_ratio>model.config.time.collision_ratio*(1+1e-12))
                           | (transport_ratio>model.config.time.transport_ratio*(1+1e-12)))
    limit = np.minimum(model.config.time.collision_ratio*tau,
                       model.config.time.transport_ratio*mesh.transport_length_m/speed)
    recommended = float(np.min(limit)) if np.all(physical) and np.any(failures) else None
    time_result = {'status':('unresolved' if not np.all(physical) else 'needs_timestep' if np.any(failures) else 'criteria_satisfied'),
                   'failed_cells':int(failures.sum()),'unresolved_cells':int((~physical).sum()),
                   'max_collision_ratio':float(np.max(collision_ratio)),'max_transport_ratio':float(np.max(transport_ratio)),
                   'collision_limit':model.config.time.collision_ratio,'transport_limit':model.config.time.transport_ratio}
    status = assessed.summary['action']
    if status=='criteria_satisfied' and time_result['status']!='criteria_satisfied':
        status = time_result['status']
    if not sampling['stationary']:
        status = 'unresolved'
    centres = mesh.geometry.centres_m[:,0].reshape(-1,mesh.transverse**2).mean(axis=1)
    predicted = model.state(centres)
    lam_profile = mean_free_path(profile['n'],profile['T'],model.source.species,convention=model.source.convention)
    ppc = np.average(raw['parcels'].reshape(-1,mesh.transverse**2),axis=1)
    profile_dict = {**profile_raw,'boundaries_m':mesh.boundaries,'r_m':centres,'n':profile['n'],
        'Ux':profile['U'][:,0],'T':profile['T'],'se':se,'blocks':blocks,
        'n_over_ne':profile['n']/model.source.n_e,'Ux_over_Ue':profile['U'][:,0]/model.source.U_e,
        'T_over_Te':profile['T']/model.source.T_e,'n_predicted':predicted['n'],
        'U_predicted':predicted['U'],'T_predicted':predicted['T'],'lambda_predicted_m':predicted['lambda'],
        'lambda_measured_m':lam_profile,'mean_parcels':ppc,
        'max_edge_m':geometry.max_edge_m.reshape(-1,mesh.transverse**2).max(axis=1),
        'max_edge_over_measured_mfp':(geometry.max_edge_m/assessed.lambda_m).reshape(-1,mesh.transverse**2).max(axis=1),
        'max_edge_over_predicted_mfp':(geometry.max_edge_m/model.state(geometry.centres_m[:,0])['lambda']).reshape(-1,mesh.transverse**2).max(axis=1),
        'max_collision_ratio':collision_ratio.reshape(-1,mesh.transverse**2).max(axis=1),
        'max_transport_ratio':transport_ratio.reshape(-1,mesh.transverse**2).max(axis=1)}
    valid_profile = assessed.valid.reshape(-1,mesh.transverse**2).all(axis=1)
    measured_b = measured_breakdown(centres,profile['n'],profile['U'][:,0],lam_profile,
          np.sqrt(8*model.source.gas_constant*profile['T']/math.pi),valid_profile,
          threshold=model.config.breakdown.threshold,density_se=se.get('n'))
    result = {'path':str(case),'case_id':record['case_id'],'Kn':model.source.knudsen,
              'diameter_m':model.source.diameter_m,'status':status,
              'physical_identity':record['physical_identity'],'sampling_identity':record['sampling'],
              'numerics':record['numerics'],'sampling':sampling,'assessment':assessed.summary,
              'measured_total_mean_parcels':float(np.sum(raw['parcels'])),
              'timestep':time_result,'recommended_delta_t_s':recommended,
              'predicted_breakdown':record['breakdown'],'measured_breakdown':measured_b,
              'reference_breakdown':record.get('reference_breakdown'),
              'native_cross_checks':{k:raw[k] for k in ('mfp','mct') if k in raw},
              'diagnostic_strict_count_valid_cells':int(np.count_nonzero(ctx['valid'])),
              'comparison_required':False,'profile':profile_dict,'plan_id':plan_identity(plan)}
    return result,(geometry,measurement,assessed,plan)


def assess_case(case):
    case = Path(case)
    try:
        result,_ = _compute(case)
        return serial(result)
    except (ValueError,OSError) as exc:
        model,mesh,r = load_case(case)
        missing = not (case/'run-state.json').exists()
        return {'path':str(case),'case_id':r['case_id'],'Kn':model.source.knudsen,
                'status':'missing_measurements' if missing else 'unresolved','reason':str(exc),
                'numerics':r['numerics'],'physical_identity':r['physical_identity'],
                'sampling_identity':r['sampling'],'predicted_breakdown':r['breakdown'],
                'reference_breakdown':r.get('reference_breakdown'),
                'comparison_required':False}


def case_plan(case):
    result,data = _compute(case)
    return data[-1],result['recommended_delta_t_s']


def checkpoint_assessment(case):
    """Assess a partial run for sizing, without certifying flow convergence."""
    try:
        result,data = _compute(case,incremental=True)
        return serial(result),data[-1] if result['status']!='unresolved' else None
    except (ValueError,OSError) as exc:
        return {'status':'unresolved','reason':str(exc)},None


def comparison_metrics(coarse,fine,species,*,targets=None,trusted=False):
    boundaries = np.asarray(coarse['boundaries_m'])
    keys = ('rho','momentum','energy')
    projected,coverage = project_moments(fine['boundaries_m'],{k:fine[k] for k in keys},boundaries)
    a,b = primitive(coarse,species),primitive(projected,species)
    volumes = np.diff(boundaries**3)
    metrics = {}
    # Project each block's moments, then recompute T; never project mean T.
    block_a = [primitive(raw,species) for raw in coarse.get('blocks',[])]
    block_b = [primitive(project_moments(fine['boundaries_m'],raw,boundaries)[0],species)
               for raw in fine.get('blocks',[])]
    for key in ('n','Ux','T'):
        values = lambda st:st['U'][:,0] if key=='Ux' else st[key]
        se = lambda blocks: np.std(np.stack([values(x) for x in blocks]),axis=0,ddof=1)/math.sqrt(len(blocks)) if len(blocks)>=2 else None
        metric = normalised_l2(values(a),values(b),volumes,se(block_a) if trusted else None,se(block_b) if trusted else None)
        target = (targets or {}).get(key)
        if target is not None:
            metric.update(target=target,verdict=verdict(metric,target) if trusted and np.all(coverage>=1-1e-10) else 'inconclusive')
        metrics[key] = metric
    metrics['covered_volume_fraction'] = float(np.sum(volumes*np.minimum(coverage,1))/volumes.sum())
    return metrics


def analyse_study(root,*,kn_values=None):
    root = Path(root)
    paths = []
    for p in sorted(root.rglob('case.json')):
        if '.stage-' in str(p):
            continue
        statepath = p.parent.parent/'controller-state.json'
        if statepath.exists() and json.loads(statepath.read_text()).get('current')!=p.parent.name:
            continue
        paths.append(p.parent)
    cases = [assess_case(p) for p in paths]
    if kn_values is not None:
        cases = [c for c in cases if c['Kn'] in kn_values]
    comparisons = []
    for a,b in itertools.combinations(cases,2):
        if a['physical_identity']!=b['physical_identity'] or a['sampling_identity']!=b['sampling_identity']:
            continue
        if 'profile' not in a or 'profile' not in b:
            continue
        same_mesh = a['numerics']['mesh_id']==b['numerics']['mesh_id']
        same_target = a['numerics']['target_parcels']==b['numerics']['target_parcels']
        if not same_mesh and not same_target:
            continue
        # Higher cell/population count supplies the reference side.
        rank = lambda c:(c['numerics']['n_cells'],1/c['numerics']['particle_weight'])
        coarse,fine = sorted((a,b),key=rank)
        model,_,_ = load_case(coarse['path'])
        tol = model.config.analysis
        targets = {k:v for k,v in zip(('n','Ux','T'),(tol.density_tolerance,tol.velocity_tolerance,tol.temperature_tolerance)) if v is not None}
        trusted = all(c['status']=='criteria_satisfied' and c['sampling']['independence_supported'] for c in (a,b))
        comparisons.append({'coarse':coarse['path'],'reference':fine['path'],
            'kind':'parcels' if same_mesh else 'mesh','trusted':trusted,
            'metrics':comparison_metrics(coarse['profile'],fine['profile'],model.source.species,targets=targets,trusted=trusted)})
    return serial({'cases':cases,'comparisons':comparisons,'scope':'radial surrogate; streamtube totals are not full-plume resource estimates'})


def post_case(case,*,plots=True):
    from .report import write_case_report
    result,data = _compute(case)
    case = Path(case)
    output = case/'postProcessing/radialstudy'/result['plan_id'][:16]
    write_assessment(*data[:3],output/'assessment')
    write_plan(data[-1],output/'plan')
    result = serial(result)
    write_case_report(output,result,plots=plots)
    return {'status':result['status'],'output':str(output),'assessment':result['assessment'],
            'timestep':result['timestep'],'sampling':result['sampling'],'recommended_delta_t_s':result['recommended_delta_t_s']}


def analyse(root):
    from .report import write_study_report
    analysis = analyse_study(root)
    output = Path(root).parent/'analysis'
    paths = write_study_report(output,analysis)
    return {'cases':len(analysis['cases']),'comparisons':len(analysis['comparisons']),'files':list(map(str,paths))}
