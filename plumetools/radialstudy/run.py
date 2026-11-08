"""Serial plus lifecycle; immutable inputs and auditable disjoint write windows."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil

import numpy as np

from plumetools.foamio.checkpoints import plan_restart
from plumetools.foamio.times import time_directories
from plumetools.refinement.geometry import file_sha256,mesh_fingerprint,read_geometry
from plumetools.refinement.openfoam import _run,_require_mesh_ok
from .case import load_case


def require_environment():
    if os.environ.get('WM_PROJECT_VERSION') not in ('v2512','2512'):
        raise ValueError('source OpenFOAM v2512/etc/bashrc before running; ensure plus binaries and libraries are on PATH/LD_LIBRARY_PATH')
    names = ('dsmcFoam+','dsmcInitialise+','blockMesh','checkMesh','postProcess')
    paths = {n:shutil.which(n) for n in names}
    if any(p is None for p in paths.values()):
        raise ValueError(f'missing native tools: {[n for n,p in paths.items() if p is None]}; source OpenFOAM then add the installed hyStrath plus binaries')
    return paths


def verify_inputs(case, record):
    changed = [p for p,h in record['input_files'].items()
               if not (case/p).is_file() or file_sha256(case/p) != h]
    if changed:
        raise ValueError(f'case inputs changed: {changed}; create a new case instead of modifying an existing run')


def _save(path, data):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
    tmp.replace(path)


def _native_run(command,cwd,log):
    timeout = None if str(command[0]) == 'dsmcFoam+' else 300
    return _run(command,cwd,log,timeout=timeout)


def mesh_case(case, *, runner=None):
    case = Path(case).resolve()
    model, mesh, record = load_case(case)
    verify_inputs(case,record)
    runner = runner or _native_run
    require_environment()
    if not (case/'constant/polyMesh/points').exists():
        if time_directories(case):
            raise ValueError('cannot build a new mesh over existing solver output')
        runner(['blockMesh','-case',str(case)],case,case/'log.blockMesh')
        runner(['checkMesh','-case',str(case)],case,case/'log.checkMesh')
        _require_mesh_ok(case/'log.checkMesh')
        for name in ('CellCentres','CellVolumes'):
            runner(['postProcess','-case',str(case),'-func','write'+name,'-time','0'],case,case/f'log.write{name}')
    g = read_geometry(case,time='0')
    if (len(g.volumes_m3) != mesh.n_cells
        or not np.allclose(g.centres_m,mesh.geometry.centres_m,rtol=1e-8,atol=1e-12)
        or not np.allclose(g.volumes_m3,mesh.geometry.volumes_m3,rtol=1e-8,atol=0)
        or not np.allclose(g.max_edge_m,mesh.geometry.max_edge_m,rtol=1e-8,atol=1e-12)):
        raise ValueError('actual mesh geometry differs from the immutable radial plan')
    return {'status':'mesh_verified','n_cells':mesh.n_cells,'mesh_sha256':g.mesh_sha256,
            'total_volume_m3':float(g.volumes_m3.sum())}


def run_case(case, *, resume=False, allow_budget_exceedance=False, runner=None):
    case = Path(case).resolve()
    model,mesh,record = load_case(case)
    verify_inputs(case,record)
    statepath = case/'run-state.json'
    if not resume and (statepath.exists() or time_directories(case) or (case/'0/lagrangian').exists()):
        raise ValueError('fresh run refuses existing output/initialization; use --resume with a complete checkpoint or generate a new revision')
    if record['numerics']['budget_gaps'] and not allow_budget_exceedance:
        raise ValueError('resource budget exceeded: '+'; '.join(record['numerics']['budget_gaps'])+'; inspect the plan, then explicitly allow budget exceedance if intended')
    runner = runner or _native_run
    environment = require_environment()
    geometry = mesh_case(case,runner=runner)
    state = {'schema_version':1,'case_id':record['case_id'],'mesh_sha256':geometry['mesh_sha256'],
             'sampling':record['sampling'],'environment':environment,'segments':[]}
    start = 0.
    if resume:
        if not statepath.is_file():
            raise ValueError('resume requires the original run-state.json and checkpoint; no reinitialization was performed')
        state = json.loads(statepath.read_text())
        if (state['case_id'] != record['case_id'] or state['mesh_sha256'] != geometry['mesh_sha256']
                or state.get('sampling',record['sampling']) != record['sampling']):
            raise ValueError('restart mesh, inputs or sampling identity changed')
        restart = plan_restart(case,1)
        if restart.refusal:
            raise ValueError(restart.refusal)
        if restart.resume is None:
            raise ValueError('no complete checkpoint for resume; no cloud was initialized or discarded')
        if restart.reached_end:
            return {'status':'already_finished','checkpoint':restart.resume}
        start = float(restart.resume)
        cadence = record['sampling']['cadence_s']
        if not math.isclose(start/cadence,round(start/cadence),abs_tol=1e-8,rel_tol=0):
            raise ValueError('checkpoint is not a native reset-at-write window boundary; independent sampling cannot resume')
    segment = {'start_s':start,'end_s':record['sampling']['end_s'],'status':'started'}
    state['segments'].append(segment)
    _save(statepath,state)
    try:
        if not resume:
            runner(['dsmcInitialise+','-case',str(case)],case,case/'log.dsmcInitialise')
        suffix = '' if not resume else f'.resume-{len(state["segments"])-1}'
        runner(['dsmcFoam+','-case',str(case)],case,case/f'log.dsmcFoam{suffix}')
        segment['status'] = 'completed'
    except Exception:
        segment['status'] = 'failed'
        _save(statepath,state)
        raise
    _save(statepath,state)
    return {'status':'completed','start_s':start,'end_s':record['sampling']['end_s'],
            'planned_steps':record['numerics']['planned_steps'],'state_file':'run-state.json'}


def validate_windows(windows):
    previous = -math.inf
    for w in windows:
        start,end = w['start_s'],w['end_s']
        if w.get('basis') != 'disjoint_native_windows':
            raise ValueError('independent blocks require disjoint native windows')
        if not math.isfinite(start+end) or start >= end or start < previous-1e-12:
            raise ValueError('invalid or overlapping sampling windows')
        previous = end


def collect_case(case):
    case = Path(case).resolve()
    model,mesh,record = load_case(case)
    verify_inputs(case,record)
    if not (case/'run-state.json').is_file():
        raise ValueError('missing native run provenance; run this case with Allrun first')
    state = json.loads((case/'run-state.json').read_text())
    sha = mesh_fingerprint(case)
    if state['case_id'] != record['case_id'] or state['mesh_sha256'] != sha:
        raise ValueError('mesh or physical/run inputs changed since sampling')
    sampling = record['sampling']
    cadence = sampling['cadence_s']
    species = model.reference.gas.species_name
    mandatory = [f'{x}_{species}' for x in ('rhoN','Ttra','dsmcNMean')]
    optional = [f'{x}_{species}' for x in ('mfp','mct','dsmcN','rhoM')]
    windows = []
    for end,name in time_directories(case):
        idx = round(end/cadence)
        start = (idx-1)*cadence
        if end <= sampling['start_s']+cadence*1e-8:
            continue
        if end > sampling['end_s']+cadence*1e-8:
            raise ValueError('output lies beyond declared sampling interval')
        if not math.isclose(end,idx*cadence,rel_tol=0,abs_tol=cadence*1e-8):
            raise ValueError(f'output {name} is not a disjoint native write boundary')
        folder = case/name
        velocity_name = next((f'{prefix}_{species}' for prefix in ('U','UMean')
                              if (folder/f'{prefix}_{species}').is_file()),None)
        fields = mandatory+([velocity_name] if velocity_name else [])
        missing = [f for f in mandatory if not (folder/f).is_file()]
        if velocity_name is None:
            missing.append(f'U_{species} (or UMean_{species})')
        if missing:
            raise ValueError(f'incomplete measurement at {name}: missing {missing}')
        files = fields+[f for f in optional if (folder/f).is_file()]
        w = {'start_s':start,'end_s':end,'time':name,'steps':record['numerics']['steps_per_window'],
             'basis':'disjoint_native_windows','source':'verified','occupancy_basis':'mean_per_cell',
             'mesh_sha256':sha,'physical_identity':record['physical_identity'],
             'case_id':record['case_id'],'block':int(round((start-sampling['start_s'])/cadence))//sampling['windows_per_block'],
             'field_roles':{'number_density':f'rhoN_{species}','temperature':f'Ttra_{species}',
                            'velocity':velocity_name,'mean_parcels':f'dsmcNMean_{species}'},
             'field_names':files,'source_files':{f'{name}/{f}':file_sha256(folder/f) for f in files}}
        windows.append(w)
    validate_windows(windows)
    _save(case/'sampling-windows.json',windows)
    return windows


def read_windows(case):
    case = Path(case).resolve()
    if not (case/'sampling-windows.json').exists():
        return collect_case(case)
    windows = json.loads((case/'sampling-windows.json').read_text())
    validate_windows(windows)
    _,_,record = load_case(case)
    verify_inputs(case,record)
    sha = mesh_fingerprint(case)
    for w in windows:
        if w['mesh_sha256'] != sha or w['case_id'] != record['case_id']:
            raise ValueError('sampling window mesh or case changed')
        for name,h in w['source_files'].items():
            if not (case/name).is_file() or file_sha256(case/name) != h:
                raise ValueError(f'sampled field {name} changed; recollect and reassess')
    cutoff = record['sampling']['start_s']+record['sampling']['cadence_s']*1e-8
    present = {name for value,name in time_directories(case) if value>cutoff}
    if present != {w['time'] for w in windows}:
        return collect_case(case)
    return windows
