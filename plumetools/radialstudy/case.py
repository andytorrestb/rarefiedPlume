"""Native plus dictionaries and atomic, idempotent case creation."""
from __future__ import annotations

from dataclasses import asdict, replace
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
import yaml
from scipy.integrate import simpson

from plumetools.refinement.geometry import file_sha256
from plumetools.refinement.plan import _publish_stage
from plumetools.refinement.openfoam import _stage_system
from .config import digest, load_study
from .physics import source_for_kn, require_same_physics, predict_breakdown, _physical_hash, _restore_source
from .mesh import render_block_mesh, write_radial_mesh, read_radial_mesh


def dictionary(name, content):
    return f'FoamFile {{ version 2.0; format ascii; class dictionary; object {name}; }}\n{content}\n'


def _write_text(path,content):
    with Path(path).open('w',encoding='utf-8',newline='\n') as stream:
        stream.write(content)


def sizing(model, mesh, target_parcels=20., particle_weight=None, delta_t_s=None):
    if not math.isfinite(target_parcels) or target_parcels <= 0:
        raise ValueError('target_parcels must be finite and positive')
    molecules = mesh.density_estimate * mesh.geometry.volumes_m3
    weight = float(np.min(molecules)/target_parcels) if particle_weight is None else float(particle_weight)
    if not math.isfinite(weight) or weight <= 0:
        raise ValueError('particle_weight must be finite and positive')
    cfg, src = model.config, model.source
    # Physical durations are independent of the mesh, weight and timestep.
    frozen = model.sampling_record
    if frozen is None:
        radii = np.geomspace(model.r_in,model.r_out,4097)
        transit = float(simpson(1/model.state(radii)['U'],x=radii))
        cadence = cfg.time.sampling_transits*transit/(cfg.time.blocks*cfg.time.windows_per_block)
        startup = math.ceil(cfg.time.startup_transits*transit/cadence-1e-12)*cadence
        sampling = cfg.time.sampling_transits*transit
        sampling_record = {'transit_s': transit, 'requested_startup_s': cfg.time.startup_transits*transit,
                           'start_s': startup, 'end_s': startup+sampling, 'cadence_s': cadence,
                           'blocks': cfg.time.blocks, 'windows_per_block': cfg.time.windows_per_block,
                           'reset_at_output': True, 'basis': 'disjoint_native_windows'}
    else:
        transit,cadence = frozen['transit_s'],frozen['cadence_s']
        if (cfg.time.blocks!=frozen['blocks'] or cfg.time.windows_per_block!=frozen['windows_per_block']
            or not math.isclose(cfg.time.startup_transits*transit,frozen['requested_startup_s'],rel_tol=1e-12)
            or not math.isclose(cfg.time.sampling_transits*transit,frozen['end_s']-frozen['start_s'],rel_tol=1e-12)):
            raise ValueError('case changes frozen physical sampling settings')
        sampling_record = dict(frozen)
    sample = model.state(mesh.geometry.centres_m[:,0])
    tau = sample['lambda']/sample['c_bar']
    speed = sample['U']+3*np.sqrt(src.gas_constant*sample['T'])
    limit = float(min(np.min(cfg.time.collision_ratio*tau),
                      np.min(cfg.time.transport_ratio*mesh.transport_length_m/speed)))
    requested = limit if delta_t_s is None else float(delta_t_s)
    if not math.isfinite(requested) or requested <= 0:
        raise ValueError('delta_t_s must be finite and positive')
    steps_window = math.ceil(cadence/requested)
    dt = cadence/steps_window
    steps = round(sampling_record['end_s']/dt)
    pop = molecules/weight
    gaps = []
    if steps > cfg.time.max_steps:
        gaps.append(f'{steps} steps exceed budget {cfg.time.max_steps}')
    if np.sum(pop) > cfg.run.max_predicted_parcels:
        gaps.append(f'{np.sum(pop):.6g} predicted parcels exceed budget {cfg.run.max_predicted_parcels}')
    numeric = {'mesh_id': digest({'boundaries': mesh.boundaries.tolist(), 'transverse': mesh.transverse,
                                 'angle_rad': mesh.angle_rad}),
               'n_cells': mesh.n_cells, 'target_parcels': target_parcels, 'particle_weight': weight,
               'density_estimate_kind': mesh.density_kind, 'predicted_min_parcels': float(pop.min()),
               'predicted_total_parcels': float(pop.sum()), 'weight_controlling_cell': int(np.argmin(molecules)),
               'delta_t_s': dt, 'predicted_delta_t_limit_s': limit, 'steps_per_window': steps_window,
               'planned_steps': steps, 'budget_gaps': gaps}
    return numeric, sampling_record


def _write_inputs(root, model, mesh, numerical, sampling, seed):
    src, gas = model.source, model.reference.gas
    species = gas.species_name
    st = {k: float(v) for k,v in model.state(model.r_in).items()}
    _stage_system(root)
    (root/'0').mkdir(exist_ok=True)
    (root/'constant').mkdir(exist_ok=True)
    content = f'''nEquivalentParticles {numerical['particle_weight']:.17g};
seedNumber {seed};
BinaryCollisionModel VariableHardSphere;
VariableHardSphereCoeffs {{ Tref {gas.t_ref_K:.17g}; }}
collisionPartnerSelectionModel noTimeCounter;
coordinateSystem dsmcCartesian;
typeIdList ({species});
moleculeProperties {{ {species} {{ mass {gas.mass_kg:.17g}; diameter {gas.diameter_m:.17g};
omega {gas.omega:.17g}; alpha 1; rotationalDegreesOfFreedom 0; }} }}'''
    _write_text(root/'constant/dsmcProperties',dictionary('dsmcProperties',content))
    _write_text(root/'constant/dynamicMeshDict',dictionary('dynamicMeshDict','dynamicFvMesh staticFvMesh;'))
    _write_text(root/'system/blockMeshDict',render_block_mesh(mesh))
    control = f'''application dsmcFoam+;
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {sampling['end_s']:.17g};
deltaT {numerical['delta_t_s']:.17g};
writeControl timeStep;
writeInterval {numerical['steps_per_window']};
purgeWrite 0;
writeFormat ascii;
writePrecision 17;
writeCompression off;
timeFormat general;
timePrecision 17;
runTimeModifiable false;
adjustTimeStep no;
nTerminalOutputs 100;'''
    _write_text(root/'system/controlDict',dictionary('controlDict',control))
    init = f'''configurations (configuration {{ type dsmcMeshFill;
numberDensities {{ {species} 0; }}
translationalTemperature {st['T']:.17g}; rotationalTemperature 0;
vibrationalTemperature 0; electronicTemperature 0; velocity (0 0 0); }});'''
    _write_text(root/'system/dsmcInitialiseDict',dictionary('dsmcInitialiseDict',init))
    patch = 'dsmcPatchBoundaries (\n'
    for name in ('inlet','outlet'):
        patch += f'boundary {{ patchBoundaryProperties {{ patchName {name}; }} boundaryModel dsmcDeletionPatch; dsmcDeletionPatchProperties {{ allSpecies yes; }} }}\n'
    patch += 'boundary { patchBoundaryProperties { patchName sides; } boundaryModel dsmcSpecularWallPatch; dsmcSpecularWallPatchProperties {} }\n);\ndsmcCyclicBoundaries ();\n'
    patch += f'''dsmcGeneralBoundaries (boundary {{ generalBoundaryProperties {{ patchName inlet; }}
boundaryModel dsmcFreeStreamInflowPatch;
dsmcFreeStreamInflowPatchProperties {{ typeIds ({species}); numberDensities {{ {species} {st['n']:.17g}; }}
translationalTemperature {st['T']:.17g}; velocity ({st['U']:.17g} 0 0); }} }});'''
    _write_text(root/'system/boundariesDict',dictionary('boundariesDict',patch))
    field = f'''dsmcFields (field {{ fieldModel dsmcVolFields;
timeProperties {{ timeOption write; resetAtOutput on; }}
dsmcVolFieldsProperties {{ fieldName {species}; typeIds ({species}); measureMeanFreePath true;
mfpReferenceTemperature {gas.t_ref_K:.17g}; averagingAcrossManyRuns false; }} }});'''
    _write_text(root/'system/fieldPropertiesDict',dictionary('fieldPropertiesDict',field))
    _write_text(root/'system/controllersDict',dictionary('controllersDict','dsmcStateControllers ();\ndsmcFluxControllers ();'))
    _write_text(root/'system/loadBalanceDict',dictionary('loadBalanceDict','enableBalancing false;\nmaximumAllowableImbalance 0.1;'))
    _write_text(root/'system/chemReactDict',dictionary('chemReactDict','reactions ();'))
    write_radial_mesh(mesh, root)
    physical = {**model.config.physical_inputs, 'exit': asdict(model.reference.exit)}
    _write_text(root/'reference.yaml',yaml.safe_dump(physical,sort_keys=False))
    study = {'model': 'cai_exit_radial_resolution', 'reference_case': 'reference.yaml',
             'knudsen_numbers': [src.knudsen], **model.config.settings()}
    study['targets']['min_mean_parcels'] = numerical['target_parcels']
    _write_text(root/'study.yaml',yaml.safe_dump(study,sort_keys=False))
    for script, cmd in (('Allmesh','mesh'),('Allrun','run'),('Allpost','post')):
        p = root/script
        _write_text(p,f'#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "$0")"\nexec "${{PYTHON:-python3}}" -m plumetools.radialstudy {cmd} "$PWD" "$@"\n')
        p.chmod(0o755)


def write_case(destination, model, mesh, *, target_parcels=20., particle_weight=None,
               delta_t_s=None, parent=None):
    destination = Path(destination).resolve()
    require_same_physics(mesh.identity, model.identity(angle_rad=mesh.angle_rad))
    num, sampling = sizing(model,mesh,target_parcels,particle_weight,delta_t_s)
    seed = 1+(model.config.run.seed+int(digest({'physical':mesh.identity,'numerics':num}),16)) % 2147483646
    parent_record = None
    if parent is not None:
        parent = Path(parent).resolve()
        old = json.loads((parent/'case.json').read_text())
        require_same_physics(old['physical_identity'],mesh.identity)
        if old['sampling'] != sampling:
            raise ValueError('numerical revision changes physical sampling times')
        parent_record = {'relative_path': os.path.relpath(parent,destination).replace('\\','/'), 'case_id': old['case_id']}
    record = {'schema_version':1, 'model':'cai_exit_radial_resolution',
              'physical_identity':mesh.identity, 'reference_sha256':model.config.reference_sha256,
              'physical_inputs':{**model.config.physical_inputs,'exit':asdict(model.reference.exit)},
              'source':model.source.as_dict(), 'inlet':{k:float(v) for k,v in model.state(model.r_in).items()},
              'numerics':num, 'sampling':sampling, 'seed':seed,
              'seed_policy':{'base_seed':model.config.run.seed,'method':'sha256 physical identity and numerical settings; serial'},
              'parent':parent_record, 'generator':'plumetools.radialstudy/1',
              'solver':{'application':'dsmcFoam+','initialiser':'dsmcInitialise+','openfoam':'v2512','execution':'serial'},
              'breakdown':predict_breakdown(model,**asdict(model.config.breakdown))}
    record['case_id'] = digest(record)
    if destination.exists():
        existing_path = destination/'case.json'
        existing = json.loads(existing_path.read_text()) if existing_path.is_file() else {}
        if existing.get('case_id') != record['case_id']:
            keys = [k for k in record if record.get(k) != existing.get(k)]
            raise ValueError(f'case conflict at {destination}: changed {keys}')
        changed = [p for p,h in existing['input_files'].items() if not (destination/p).is_file() or file_sha256(destination/p) != h]
        if changed:
            raise ValueError(f'case input conflict: {changed}')
        return existing
    destination.parent.mkdir(parents=True,exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=destination.name+'.stage-',dir=destination.parent))
    try:
        _write_inputs(stage,model,mesh,num,sampling,seed)
        record['input_files'] = {p.relative_to(stage).as_posix():file_sha256(p) for p in sorted(stage.rglob('*')) if p.is_file()}
        _write_text(stage/'case.json',json.dumps(record,sort_keys=True,indent=2)+'\n')
        _publish_stage(stage,destination)
    except Exception:
        shutil.rmtree(stage)
        raise
    return record


def load_case(case):
    case = Path(case)
    record = json.loads((case/'case.json').read_text())
    if digest({k:v for k,v in record.items() if k not in ('case_id','input_files')}) != record['case_id']:
        raise ValueError('case record no longer matches its case identity')
    cfg = replace(load_study(case/'study.yaml'),reference_sha256=record['reference_sha256'])
    model = source_for_kn(cfg,record['source']['exit']['Kn'])
    mesh = read_radial_mesh(case)
    identity = record['physical_identity']
    if not (math.isclose(model.r_in,identity['r_in_m'],rel_tol=1e-12)
            and math.isclose(model.r_out,identity['r_out_m'],rel_tol=1e-12)):
        raise ValueError('case domain differs from its frozen physical identity')
    # Rehash physical inputs with the recorded domain, rather than a fresh
    # optimizer's last-bit differences on another Python/SciPy platform.
    fingerprint = _physical_hash(cfg,model.source.knudsen,identity['r_in_m'],identity['r_out_m'])
    model = replace(model,source=_restore_source(model.source,record['source']),
                    r_in=identity['r_in_m'],r_out=identity['r_out_m'],physical_hash=fingerprint,
                    inlet_state=record['inlet'],sampling_record=record['sampling'])
    require_same_physics(record['physical_identity'],model.identity(angle_rad=mesh.angle_rad))
    return model, mesh, record
