"""Serial numerical sizing with durable summaries and clean, owned retries."""
from __future__ import annotations

from dataclasses import replace
from datetime import datetime,timezone
import json
import math
import os
from pathlib import Path
import re
import time
import uuid

from plumetools.foamio.checkpoints import plan_restart
from plumetools.foamio.times import time_directories
from plumetools.refinement.geometry import file_sha256,mesh_fingerprint
from plumetools.refinement.plan import read_plan,write_plan,plan_identity,json_safe
from .case import load_case
from .generation import revise_case
from .post import checkpoint_assessment
from .run import _save,_has_output,_restart_state
from .segments import run_segment,recover_control,verify_segment_inputs
from .locking import exclusive_lock


class AdaptiveStop(RuntimeError):
    """An honest non-passing outcome, with the evidence retained."""


def _attempt(root,name):
    if not re.fullmatch(r'r\d+_ppc[0-9pem+]+',name):
        raise ValueError('invalid controller attempt path')
    path = root/name
    if path.is_symlink() or path.resolve().parent!=root.resolve():
        raise ValueError('controller attempt must be directly inside its owned Kn directory')
    return path


def _statepath(root):
    return root/'controller-state.json'


def _log(root,event,**values):
    path = root.parent/'adaptive-study.log'
    stamp = datetime.now(timezone.utc).isoformat(timespec='seconds')
    data = json_safe(values)
    line = stamp+' '+event+' '+' '.join(f'{k}={json.dumps(v,separators=(",",":"),allow_nan=False)}' for k,v in data.items())+'\n'
    with exclusive_lock(root.parent/'.adaptive-log.lock',blocking=True):
        with path.open('a',encoding='utf-8',newline='\n') as stream:
            stream.write(line)
            stream.flush()
            os.fsync(stream.fileno())


def _owner(case,state):
    marker = case/'controller-owner.json'
    expected = {'owner_id':state['owner_id'],'case_id':load_case(case)[2]['case_id']}
    if marker.exists() and json.loads(marker.read_text())!=expected:
        raise ValueError('attempt is owned by a different controller')
    _save(marker,expected)


def _read_state(root,record):
    state = json.loads(_statepath(root).read_text())
    if state.get('schema_version')!=1 or state.get('physical_identity')!=record['physical_identity']:
        raise ValueError('controller physical identity changed; generate a fresh study tree')
    for name,item in state['attempts'].items():
        path = _attempt(root,name)
        _,_,r = load_case(path)
        verify_segment_inputs(path,r)
        if r['case_id']!=item['case_id']:
            raise ValueError('controller attempt inputs changed')
        marker = path/'controller-owner.json'
        if not marker.is_file() or json.loads(marker.read_text())!={'owner_id':state['owner_id'],'case_id':r['case_id']}:
            raise ValueError('missing or changed controller ownership record')
    if state['current'] not in state['attempts']:
        raise ValueError('controller has no recognized current attempt')
    return state


def _checkpoint(case,record):
    if not (case/'run-state.json').exists():
        if _has_output(case):
            raise ValueError('output has no native run provenance; refusing a fresh initialization')
        return 0.
    _,restart = _restart_state(case,record,mesh_fingerprint(case))
    return float(restart.resume)


def _verify_accepted(case,state,record):
    certificate = state['certificate']
    if mesh_fingerprint(case)!=certificate['mesh_sha256']:
        raise ValueError('accepted mesh changed')
    for name,expected in certificate['source_files'].items():
        path = (case/name).resolve()
        if not path.is_relative_to(case.resolve()) or not path.is_file() or file_sha256(path)!=expected:
            raise ValueError(f'accepted measurement or input {name} changed')
    latest = _checkpoint(case,record)
    if not math.isclose(latest,state['accepted_time_s'],rel_tol=0,abs_tol=record['numerics']['delta_t_s']/2):
        raise ValueError('accepted checkpoint changed')


def adaptive_preview(case,*,resume=False):
    """Read-only queue entry, including recovery and early-acceptance checks."""
    case = Path(case).resolve()
    model,_,record = load_case(case)
    verify_segment_inputs(case,record)
    if not model.config.controller.enabled:
        raise ValueError('checkpoint controller is not enabled for this case')
    root = case.parent
    state = None
    if _statepath(root).exists():
        if not resume:
            raise ValueError('existing controller state; use Allrestart or --resume')
        state = _read_state(root,record)
        case = _attempt(root,state['current'])
        model,_,record = load_case(case)
        if state['outcome']=='accepted':
            _verify_accepted(case,state,record)
    elif _has_output(case):
        raise ValueError('existing unowned run; generate a fresh Cases-adaptive tree before using the controller')
    action = ('skip' if state and state['outcome']=='accepted' else
              'retry' if state and state.get('pending_retry') else 'resume' if state else 'fresh')
    start = _checkpoint(case,record) if action not in ('skip','retry') else 0.
    next_end = min(record['sampling']['end_s'],
                   (round(start/record['sampling']['cadence_s'])+model.config.controller.check_every_windows)*record['sampling']['cadence_s'])
    gaps = []
    if action!='skip':
        if record['numerics']['predicted_total_parcels']>model.config.run.max_predicted_parcels:
            gaps.append('predicted parcel budget exceeded')
        charged = state['total_steps'] if state else 0
        steps = max(0,round((next_end-start)/record['numerics']['delta_t_s']))
        if charged+steps>model.config.time.max_steps:
            gaps.append('cumulative step budget exceeded')
    return {'path':str(case),'Kn':record['source']['exit']['Kn'],'action':action,
            'outcome':state['outcome'] if state else 'untouched',
            'retries':len(state['attempts'])-1 if state else 0,
            'n_cells':record['numerics']['n_cells'],'half_angle_deg':model.config.domain.angle_deg,
            'outer_over_D':model.config.domain.outer_over_D,'next_end_s':next_end,
            'budget_gaps':gaps,'numerics':record['numerics']}


def clean_attempt(case,root,owner_id):
    """Delete only native outputs of a durably superseded owned attempt."""
    case,root = Path(case),Path(root)
    supplied = case.resolve()
    case = _attempt(root,case.name)
    if supplied!=case.resolve():
        raise ValueError('cleanup source must be inside its owned Kn directory')
    marker = case/'controller-owner.json'
    if not marker.is_file() or json.loads(marker.read_text()).get('owner_id')!=owner_id:
        raise ValueError('cleanup requires a controller-owned attempt')
    state = json.loads(_statepath(root).read_text())
    pending = state.get('pending_retry')
    if (not pending or pending.get('phase')!='cleanup' or pending['parent']!=case.name
            or state['owner_id']!=owner_id or state['outcome']=='accepted'):
        raise ValueError('cleanup requires a committed owned retry transition')
    child = _attempt(root,pending['child'])
    _,_,record = load_case(child)
    verify_segment_inputs(child,record)
    _output_files(case)  # Validate all output containers, including unknown links.
    inventory = json.loads((case/'controller-outputs.json').read_text())
    targets = []
    for name,expected in inventory.items():
        path = case/name
        if path.is_symlink() or not path.resolve().is_relative_to(case.resolve()):
            raise ValueError(f'cleanup refuses linked or outside output: {path}')
        if path.exists():
            if file_sha256(path)!=expected:
                raise ValueError(f'owned native output changed before cleanup: {path}')
            targets.append(path)
    directories = set()
    for path in targets:
        directories.update(p for p in path.parents if p!=case and p.is_relative_to(case))
        path.unlink()
    for directory in sorted(directories,key=lambda p:len(p.parts),reverse=True):
        if directory.exists() and not any(directory.iterdir()):
            directory.rmdir()


def _output_files(case):
    containers = [case/'constant/polyMesh',case/'0/C',case/'0/V',case/'0/lagrangian',
                  case/'postProcessing',case/'fieldMeasurements',case/'boundaries']
    containers += [case/name for _,name in time_directories(case)]
    containers += [p for p in case.glob('log.*') if p.name in (
        'log.blockMesh','log.checkMesh','log.writeCellCentres','log.writeCellVolumes','log.dsmcInitialise')
        or re.fullmatch(r'log\.dsmcFoam\.stage-\d+',p.name)]
    files = []
    for path in containers:
        if path.is_symlink() or getattr(path,'is_junction',lambda:False)() or not path.resolve().is_relative_to(case.resolve()):
            raise ValueError(f'cleanup refuses linked or outside output: {path}')
        entries = [path]+list(path.rglob('*')) if path.is_dir() else [path]
        for entry in entries:
            if (entry.is_symlink() or getattr(entry,'is_junction',lambda:False)()
                    or not entry.resolve().is_relative_to(case.resolve())):
                raise ValueError(f'cleanup refuses linked or outside output: {entry}')
            if entry.is_file():
                files.append(entry)
    return files


def _capture_outputs(case,before):
    path = case/'controller-outputs.json'
    owned = json.loads(path.read_text()) if path.exists() else {}
    for output in _output_files(case):
        name = output.relative_to(case).as_posix()
        sha = file_sha256(output)
        if before.get(name)!=sha:
            owned[name] = sha
    _save(path,owned)


def _stop(root,state,outcome,reason):
    state['outcome'] = outcome
    state['reason'] = reason
    _save(_statepath(root),state)
    _log(root,outcome,Kn=state['Kn'],attempt=state['current'],reason=reason)
    raise AdaptiveStop(reason)


def _finish_retry(root,state):
    pending = state['pending_retry']
    parent = _attempt(root,pending['parent'])
    child = _attempt(root,pending['child'])
    plan = read_plan(parent/pending['plan'])
    if not child.exists():
        row = revise_case(parent,plan,delta_t_s=pending['delta_t_s'])
        if Path(row['path'])!=child:
            raise ValueError('retry destination changed during publication; preserve both attempts')
    _,_,record = load_case(child)
    revision = json.loads((child/'revision.json').read_text())
    if (record['parent']['case_id']!=state['attempts'][parent.name]['case_id']
            or record['physical_identity']!=state['physical_identity']
            or revision['parent_plan_id']!=plan_identity(plan)):
        raise ValueError('prepared retry does not match the committed correction')
    verify_segment_inputs(child,record)
    _owner(child,state)
    state['attempts'][child.name] = {'case_id':record['case_id'],'status':'active'}
    pending['phase'] = 'cleanup'
    _save(_statepath(root),state)
    _log(root,'retry',Kn=state['Kn'],attempt=parent.name,next_attempt=child.name,
         reason=pending['reason'],old_settings=pending['old_settings'],new_settings=record['numerics'])
    clean_attempt(parent,root,state['owner_id'])
    _log(root,'cleanup',Kn=state['Kn'],attempt=parent.name)
    state['attempts'][parent.name]['status'] = 'superseded'
    state.update(current=child.name,pending_retry=None,last_checked_s=0.,outcome='running')
    _save(_statepath(root),state)


def run_adaptive(case,*,resume=False,allow_budget_exceedance=False,progress=None):
    """Run/retry one Kn case; the persistent study log survives every retry."""
    case = Path(case).resolve()
    legacy_lock = case.parent.parent/'.adaptive-study.lock'
    if legacy_lock.exists():
        with exclusive_lock(legacy_lock):
            pass  # Refuse to race a still-running older serial controller.
    with exclusive_lock(case.parent/'.adaptive-controller.lock'):
        return _run_adaptive(case,resume=resume,allow_budget_exceedance=allow_budget_exceedance,progress=progress)


def _run_adaptive(case,*,resume=False,allow_budget_exceedance=False,progress=None):
    case = Path(case).resolve()
    preview = adaptive_preview(case,resume=resume)
    if preview['action']=='skip':
        return {'status':'already_finished','path':preview['path']}
    root = case.parent
    model,_,record = load_case(case)
    if _statepath(root).exists():
        state = _read_state(root,record)
    else:
        state = {'schema_version':1,'owner_id':uuid.uuid4().hex,'Kn':preview['Kn'],
                 'physical_identity':record['physical_identity'],'current':case.name,
                 'attempts':{case.name:{'case_id':record['case_id'],'status':'active'}},
                 'total_steps':0,'wall_seconds':0.,'last_checked_s':0.,'outcome':'running','pending_retry':None}
        _owner(case,state)
        _save(_statepath(root),state)
        _log(root,'start',Kn=state['Kn'],attempt=case.name,numerics=record['numerics'])
    try:
        while True:
            if state.get('pending_retry'):
                _finish_retry(root,state)
            active = _attempt(root,state['current'])
            recover_control(active)
            model,_,record = load_case(active)
            verify_segment_inputs(active,record)
            current = _checkpoint(active,record)
            cadence,limit = record['sampling']['cadence_s'],record['sampling']['end_s']
            if current>record['sampling']['start_s']+cadence*1e-8 and current>state['last_checked_s']+cadence*1e-8:
                result,plan = checkpoint_assessment(active)
                assessment,timestep = result.get('assessment',{}),result.get('timestep',{})
                summary = {'status':result['status'],'sampling':result.get('sampling',{}),
                           'assessment':assessment,'timestep':timestep,'reason':result.get('reason')}
                _save(active/'controller-assessment.json',json_safe(summary))
                _log(root,'assessment',Kn=state['Kn'],attempt=active.name,time_s=current,
                     wall_seconds=state['wall_seconds'],steps=state['total_steps'],cells=record['numerics']['n_cells'],
                     particle_weight=record['numerics']['particle_weight'],delta_t_s=record['numerics']['delta_t_s'],
                     min_mean_parcels=assessment.get('min_mean_parcels'),max_edge_over_mfp=assessment.get('max_edge_over_mfp'),
                     mesh_failed=assessment.get('mesh_failed_cells'),occupancy_failed=assessment.get('occupancy_failed_cells'),
                     timestep_failed=timestep.get('failed_cells'),unresolved=assessment.get('unresolved_cells'),
                     action=result['status'],reason=result.get('reason',result.get('sampling',{}).get('reason')))
                if progress:
                    progress(f"Kn{state['Kn']:g} {active.name} t={current:g}: {result['status']}")
                if result['status']=='criteria_satisfied' and plan is not None:
                    state.update(outcome='accepted',accepted_time_s=current,last_checked_s=current,
                                 certificate={'mesh_sha256':plan.mesh_sha256,'source_files':plan.provenance['source_files']})
                    _save(_statepath(root),state)
                    _log(root,'accepted',Kn=state['Kn'],attempt=active.name,time_s=current)
                    return {'status':'accepted','path':str(active),'retries':len(state['attempts'])-1,'steps':state['total_steps']}
                if plan is not None:
                    if len(state['attempts'])-1>=model.config.controller.max_retries:
                        _stop(root,state,'budget_limited','automatic retry limit reached')
                    if plan.action=='needs_parcels':
                        target = max(plan.target_parcels,model.config.controller.retry_parcels)
                        plan = replace(plan,recommended_weight=record['numerics']['particle_weight']*assessment['min_mean_parcels']/target)
                    try:
                        proposed = revise_case(active,plan,delta_t_s=result.get('recommended_delta_t_s'),dry_run=True)
                    except ValueError as exc:
                        if 'budget' in str(exc):
                            _stop(root,state,'budget_limited',str(exc))
                        raise
                    if proposed['action']!='create_revision':
                        _stop(root,state,'unresolved','measured failure has no executable correction')
                    num = proposed['numerics']
                    if not allow_budget_exceedance and num['predicted_total_parcels']>model.config.run.max_predicted_parcels:
                        _stop(root,state,'budget_limited','corrected parcel budget exceeded')
                    if not allow_budget_exceedance and state['total_steps']+round(cadence/num['delta_t_s'])>model.config.time.max_steps:
                        _stop(root,state,'budget_limited','cumulative step budget cannot start corrected attempt')
                    folder = f'controller-plans/check-{round(current/cadence)}'
                    write_plan(plan,active/folder)
                    state['pending_retry'] = {'parent':active.name,'child':Path(proposed['path']).name,
                        'plan':folder,'delta_t_s':result.get('recommended_delta_t_s'),'phase':'intent',
                        'reason':result['status'],'old_settings':record['numerics']}
                    state['last_checked_s'] = current
                    _save(_statepath(root),state)
                    continue
                state['last_checked_s'] = current
                _save(_statepath(root),state)
            if current>=limit-record['numerics']['delta_t_s']/2:
                _stop(root,state,'unresolved','maximum sampling limit reached without trusted numerical acceptance')
            end = min(limit,(round(current/cadence)+model.config.controller.check_every_windows)*cadence)
            steps = round((end-current)/record['numerics']['delta_t_s'])
            if not allow_budget_exceedance:
                if state['total_steps']+steps>model.config.time.max_steps:
                    _stop(root,state,'budget_limited','cumulative step budget exceeded before next segment')
                if record['numerics']['predicted_total_parcels']>model.config.run.max_predicted_parcels:
                    _stop(root,state,'budget_limited','predicted parcel budget exceeded')
            state['total_steps'] += steps  # Conservative charge survives partial writes/process loss.
            state['outcome'] = 'running'
            _save(_statepath(root),state)
            started = time.monotonic()
            before = {p.relative_to(active).as_posix():file_sha256(p) for p in _output_files(active)}
            try:
                run_segment(active,end)
            finally:
                _capture_outputs(active,before)
                state['wall_seconds'] += time.monotonic()-started
                _save(_statepath(root),state)
            _log(root,'checkpoint',Kn=state['Kn'],attempt=active.name,time_s=end,steps=state['total_steps'])
    except AdaptiveStop:
        raise
    except BaseException as exc:
        state['outcome'] = 'failed'
        _save(_statepath(root),state)
        _log(root,'solver_failure',Kn=state['Kn'],attempt=state['current'],reason=str(exc))
        raise
