"""Native write-boundary stops with narrowly audited controlDict overrides."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import time

from .case import load_case
from .locking import exclusive_lock
from .run import (_save,_native_run,_has_output,_restart_state,mesh_case,
                  require_environment,verify_inputs)
from plumetools.refinement.geometry import file_sha256


def _canonical(case,record):
    manifest = case/'controller-control.json'
    data = json.loads(manifest.read_text())
    canonical = data['canonical'].encode('utf-8')
    expected = record['input_files']['system/controlDict']
    if (hashlib.sha256(canonical).hexdigest()!=expected or data['canonical_sha256']!=expected
            or file_sha256(case/'system/controlDict') not in (expected,data['active_sha256'])):
        raise ValueError('controlDict changed outside the recorded controller override; preserve it and inspect controller-control.json')
    return canonical


def verify_segment_inputs(case,record):
    """Read-only validation, including a pending exact controller override."""
    case = Path(case)
    if not (case/'controller-control.json').exists():
        return verify_inputs(case,record)
    _canonical(case,record)
    verify_inputs(case,{**record,'input_files':{k:v for k,v in record['input_files'].items()
                                              if k!='system/controlDict'}})


def _recover_control(case):
    case = Path(case)
    manifest = case/'controller-control.json'
    if manifest.exists():
        _,_,record = load_case(case)
        verify_segment_inputs(case,record)
        canonical = _canonical(case,record)
        path = case/'system/controlDict'
        temporary = path.with_suffix('.controller.tmp')
        temporary.write_bytes(canonical)
        temporary.replace(path)
        manifest.unlink()


def recover_control(case):
    case = Path(case)
    with exclusive_lock(case/'controller.lock'):
        _recover_control(case)


def run_segment(case,end_s,*,runner=None):
    case = Path(case)
    with exclusive_lock(case/'controller.lock'):
        return _run_segment(case,end_s,runner=runner)


def _run_segment(case,end_s,*,runner=None):
    """Initialize once, then continue unchanged native clouds to exact writes."""
    case = Path(case).resolve()
    _recover_control(case)
    _,_,record = load_case(case)
    verify_inputs(case,record)
    cadence = record['sampling']['cadence_s']
    index = end_s/cadence
    if (not math.isfinite(end_s) or end_s<=0 or
            end_s>record['sampling']['end_s']+record['numerics']['delta_t_s']/2 or
            not math.isclose(index,round(index),abs_tol=1e-8,rel_tol=0)):
        raise ValueError('segment end must be a native write boundary within the sampling limit')
    statepath = case/'run-state.json'
    fresh = not _has_output(case)
    if not fresh and not statepath.exists():
        raise ValueError('existing initialization/output has no original run-state.json; refusing reinitialization')
    environment = require_environment()
    runner = runner or _native_run
    geometry = mesh_case(case,runner=runner)
    state = {'schema_version':1,'case_id':record['case_id'],'mesh_sha256':geometry['mesh_sha256'],
             'sampling':record['sampling'],'environment':environment,'segments':[]}
    start = 0.
    if not fresh:
        state,restart = _restart_state(case,record,geometry['mesh_sha256'])
        start = float(restart.resume)
    if end_s<=start+record['numerics']['delta_t_s']/2:
        raise ValueError('segment requires a checkpoint later than the existing one')
    segment = {'start_s':start,'end_s':end_s,'status':'started',
               'steps':round((end_s-start)/record['numerics']['delta_t_s'])}
    state['segments'].append(segment)
    _save(statepath,state)
    control = case/'system/controlDict'
    original = control.read_bytes()
    active,count = re.subn(rb'(?m)^endTime\s+[^;]+;',
                          f'endTime {end_s:.17g};'.encode(),original)
    if count!=1:
        raise ValueError('expected one canonical endTime for staged native execution')
    _save(case/'controller-control.json',{'canonical':original.decode('utf-8'),
          'canonical_sha256':hashlib.sha256(original).hexdigest(),
          'active_sha256':hashlib.sha256(active).hexdigest()})
    try:
        temporary = control.with_suffix('.controller.tmp')
        temporary.write_bytes(active)
        temporary.replace(control)
        if fresh:
            runner(['dsmcInitialise+','-case',str(case)],case,case/'log.dsmcInitialise')
        log = case/f'log.dsmcFoam.stage-{len(state["segments"])-1}'
        solver_started = time.monotonic()
        try:
            runner(['dsmcFoam+','-case',str(case)],case,log)
        finally:
            segment['solver_wall_seconds'] = time.monotonic()-solver_started
        _recover_control(case)
        _,checkpoint = _restart_state(case,record,geometry['mesh_sha256'])
        if not math.isclose(float(checkpoint.resume),end_s,abs_tol=record['numerics']['delta_t_s']/2,rel_tol=0):
            raise ValueError(f'solver did not write the requested complete checkpoint at {end_s:g}; see {log}')
        segment['status'] = 'completed'
    except BaseException:
        segment['status'] = 'failed'
        raise
    finally:
        _save(statepath,state)
        _recover_control(case)
    return dict(segment)
