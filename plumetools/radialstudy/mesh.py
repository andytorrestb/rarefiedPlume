"""Conformal square streamtubes with actual edge bounds and immutable domains."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import cached_property
import json
import math
from pathlib import Path
import shutil
import tempfile

import numpy as np

from plumetools.refinement.model import Geometry, RefinementResult
from plumetools.refinement.plan import verify_sources, _publish_stage
from plumetools.refinement.geometry import mesh_fingerprint
from plumetools.refinement.openfoam import _run, _stage_system, _require_mesh_ok
from .config import digest

HEADER = 'FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }\n'


@dataclass(frozen=True)
class RadialMesh:
    boundaries: np.ndarray
    transverse: int
    angle_rad: float
    identity: dict
    layer_density: np.ndarray
    density_kind: str = 'predicted'
    max_cells: int = 1000000

    def __post_init__(self):
        object.__setattr__(self, 'boundaries', np.asarray(self.boundaries, dtype=float))
        object.__setattr__(self, 'layer_density', np.asarray(self.layer_density, dtype=float))
        if (self.boundaries.ndim != 1 or len(self.boundaries) < 2
                or not np.all(np.isfinite(self.boundaries)) or np.any(np.diff(self.boundaries) <= 0)
                or self.boundaries[0] <= 0 or self.transverse < 1
                or not 0 < self.angle_rad < math.pi/4
                or self.layer_density.shape != (len(self.boundaries)-1,)
                or not np.all(np.isfinite(self.layer_density) & (self.layer_density > 0))):
            raise ValueError('invalid radial mesh geometry or density estimates')
        if self.n_cells > self.max_cells:
            raise ValueError(f'mesh requires {self.n_cells} cells; budget is {self.max_cells}; no coarsening applied')

    @property
    def n_cells(self):
        return (len(self.boundaries)-1) * self.transverse**2

    @cached_property
    def vertices(self):
        t, n = math.tan(self.angle_rad), self.transverse
        a, b = np.repeat(self.boundaries[:-1], n*n), np.repeat(self.boundaries[1:], n*n)
        j = np.tile(np.arange(n), (len(self.boundaries)-1)*n)
        k = np.tile(np.repeat(np.arange(n), n), len(self.boundaries)-1)
        y0, y1 = -1+2*j/n, -1+2*(j+1)/n
        z0, z1 = -1+2*k/n, -1+2*(k+1)/n
        return np.stack([np.column_stack((r, r*t*y, r*t*z)) for r,y,z in (
            (a,y0,z0),(b,y0,z0),(b,y1,z0),(a,y1,z0),
            (a,y0,z1),(b,y0,z1),(b,y1,z1),(a,y1,z1))], axis=1)

    @cached_property
    def geometry(self):
        a, b = self.boundaries[:-1], self.boundaries[1:]
        n, t = self.transverse, math.tan(self.angle_rad)
        # Stable polynomial differences, important for dense near-inlet cells.
        xc = .75 * (a+b)*(a*a+b*b)/(a*a+a*b+b*b)
        v = 4*t*t/(3*n*n) * (b-a)*(a*a+a*b+b*b)
        j = np.tile(np.arange(n), (len(a))*n)
        k = np.tile(np.repeat(np.arange(n), n), len(a))
        x = np.repeat(xc, n*n)
        c = np.column_stack((x, x*t*(-1+(2*j+1)/n), x*t*(-1+(2*k+1)/n)))
        ends = ([1,2,3,0,5,6,7,4,4,5,6,7], [0,1,2,3,4,5,6,7,0,1,2,3])
        edge = np.linalg.norm(self.vertices[:,ends[0]]-self.vertices[:,ends[1]], axis=2).max(axis=1)
        return Geometry(c, np.repeat(v, n*n), edge, digest(self.metadata()))

    @property
    def transport_length_m(self):
        t = math.tan(self.angle_rad)
        length = np.minimum(np.diff(self.boundaries), 2*self.boundaries[:-1]*t/self.transverse)
        return np.repeat(length / math.sqrt(1+2*t*t), self.transverse**2)

    @property
    def density_estimate(self):
        return np.repeat(self.layer_density, self.transverse**2)

    def metadata(self):
        return {'schema_version': 1, 'boundaries_m': self.boundaries.tolist(),
                'transverse': self.transverse, 'angle_rad': self.angle_rad,
                'identity': self.identity, 'layer_density': self.layer_density.tolist(),
                'density_kind': self.density_kind, 'max_cells': self.max_cells}


def initial_mesh(model):
    cfg, d = model.config, model.source.diameter_m
    if cfg.mesh.initial_cells is not None:
        bounds = np.linspace(model.r_in,model.r_out,cfg.mesh.initial_cells+1)
        angle = math.radians(cfg.domain.angle_deg if cfg.domain.angle_deg is not None else .5)
        centres = .75*(bounds[:-1]+bounds[1:])*(bounds[:-1]**2+bounds[1:]**2)/(bounds[:-1]**2+bounds[:-1]*bounds[1:]+bounds[1:]**2)
        return RadialMesh(bounds,1,angle,model.identity(angle_rad=angle),
                          model.state(centres)['n'],max_cells=cfg.mesh.max_cells)
    samples = np.geomspace(model.r_in, model.r_out, 8193)
    lam = np.asarray(model.state(samples)['lambda'])
    if not np.all(np.isfinite(lam) & (lam > 0)) or np.any(np.diff(lam) < 0):
        raise ValueError('predicted MFP must be positive and monotone for radial meshing')
    target = np.minimum(cfg.mesh.predicted_edge_over_mfp*lam, cfg.mesh.geometric_cap_over_D*d)
    angle = (math.radians(cfg.domain.angle_deg) if cfg.domain.angle_deg is not None
             else math.atan(min(math.tan(math.radians(.5)), .25*float(np.min(target/samples)))))
    factor = math.sqrt(1+2*math.tan(angle)**2)
    bounds = [model.r_in]
    while bounds[-1] < model.r_out:
        idx = max(0, np.searchsorted(samples, bounds[-1], side='right')-1)
        nxt = min(model.r_out, bounds[-1]+target[idx]/factor*(1-1e-10))
        if nxt <= bounds[-1]:
            raise ValueError('radial step below coordinate precision')
        bounds.append(nxt)
        if len(bounds)-1 > cfg.mesh.max_cells:
            raise ValueError('radial cell budget exceeded; no coarsening applied')
    bounds = np.array(bounds)
    # Shared transverse subdivision is required for conforming block faces.
    bound = np.minimum(cfg.mesh.predicted_edge_over_mfp*model.state(bounds[:-1])['lambda'],
                       cfg.mesh.geometric_cap_over_D*d)
    nt = max(1, math.ceil(float(np.max(2*bounds[1:]*math.tan(angle)/bound))))
    centres = .75*(bounds[:-1]+bounds[1:])*(bounds[:-1]**2+bounds[1:]**2)/(bounds[:-1]**2+bounds[:-1]*bounds[1:]+bounds[1:]**2)
    return RadialMesh(bounds, nt, angle, model.identity(angle_rad=angle),
                      model.state(centres)['n'], max_cells=cfg.mesh.max_cells)


def refine_mesh(parent, plan):
    if plan.action == 'criteria_satisfied':
        return parent
    if plan.action != 'needs_mesh':
        raise ValueError(f'no trusted mesh correction: {plan.action}')
    if len(plan.target_edge_m) != parent.n_cells:
        raise ValueError('source mesh does not match radial metadata')
    bounds = np.asarray(plan.target_edge_m).reshape(-1, parent.transverse**2).min(axis=1)
    density = np.asarray(plan.density_estimate).reshape(-1, parent.transverse**2).min(axis=1)
    if not np.all(np.isfinite(bounds) & (bounds > 0)) or not np.all(np.isfinite(density) & (density > 0)):
        raise ValueError('incomplete measured coverage; no mesh estimate')
    t = math.tan(parent.angle_rad)
    factor = math.sqrt(1+2*t*t)
    nt = max(parent.transverse, math.ceil(float(np.max(2*parent.boundaries[1:]*t/bounds))))
    cuts = np.maximum(1, np.ceil(np.diff(parent.boundaries)*factor/bounds*(1-1e-12)).astype(int))
    count = int(np.sum(cuts))*nt*nt
    if count > parent.max_cells:
        raise ValueError(f'refinement requires {count} cells; budget {parent.max_cells}')
    new = np.concatenate([np.linspace(a,b,int(n)+1)[:-1] for a,b,n in
                          zip(parent.boundaries[:-1], parent.boundaries[1:], cuts)] + [parent.boundaries[-1:]])
    return RadialMesh(new, nt, parent.angle_rad, parent.identity, np.repeat(density, cuts),
                      'parent-measured', parent.max_cells)


def render_block_mesh(mesh):
    t, n = math.tan(mesh.angle_rad), mesh.transverse
    lines = [HEADER, 'convertToMeters 1;\nvertices\n(']
    for r in mesh.boundaries:
        for y,z in ((-1,-1),(1,-1),(1,1),(-1,1)):
            lines.append(f'({r:.17g} {y*r*t:.17g} {z*r*t:.17g})')
    lines += [');\nblocks\n(']
    for i in range(len(mesh.boundaries)-1):
        a,b = 4*i,4*(i+1)
        ids = (a,b,b+1,a+1,a+3,b+3,b+2,a+2)
        lines.append(f'hex ({" ".join(map(str,ids))}) (1 {n} {n}) simpleGrading (1 1 1)')
    lines += [');\nedges ();\nboundary\n(', 'inlet { type patch; faces ((0 3 2 1)); }']
    b = 4*(len(mesh.boundaries)-1)
    lines.append(f'outlet {{ type patch; faces (({b} {b+1} {b+2} {b+3})); }}')
    lines.append('sides { type wall; faces (')
    for i in range(len(mesh.boundaries)-1):
        a,b = 4*i,4*(i+1)
        for ids in ((a,a+1,b+1,b),(a+1,a+2,b+2,b+1),(a+2,a+3,b+3,b+2),(a+3,a,b,b+3)):
            lines.append('('+' '.join(map(str,ids))+')')
    lines += ['); }', ');\nmergePatchPairs ();\n']
    return '\n'.join(lines)


def write_radial_mesh(mesh, artifact):
    artifact = Path(artifact)
    artifact.mkdir(parents=True, exist_ok=True)
    (artifact/'radial-mesh.json').write_text(json.dumps(mesh.metadata(), sort_keys=True, indent=2)+'\n')


def read_radial_mesh(artifact):
    raw = json.loads((Path(artifact)/'radial-mesh.json').read_text())
    if raw['schema_version'] != 1:
        raise ValueError('unsupported radial mesh schema')
    return RadialMesh(raw['boundaries_m'], raw['transverse'], raw['angle_rad'], raw['identity'],
                      raw['layer_density'], raw['density_kind'], raw['max_cells'])


def materialize(mesh, destination):
    destination = Path(destination)
    if destination.exists():
        raise ValueError(f'destination exists: {destination}')
    destination.mkdir(parents=True)
    _stage_system(destination)
    (destination/'system/blockMeshDict').write_text(render_block_mesh(mesh))
    write_radial_mesh(mesh, destination)
    (destination/'0').mkdir(exist_ok=True)
    _run(['blockMesh', '-case', str(destination)], destination, destination/'log.blockMesh')
    _run(['checkMesh', '-case', str(destination)], destination, destination/'log.checkMesh')
    _require_mesh_ok(destination/'log.checkMesh')
    _run(['postProcess', '-case', str(destination), '-func', 'writeCellCentres', '-time', '0'], destination, destination/'log.C')
    _run(['postProcess', '-case', str(destination), '-func', 'writeCellVolumes', '-time', '0'], destination, destination/'log.V')


class RadialRefiner:
    def refine(self, source, plan, destination):
        if plan.action == 'criteria_satisfied':
            return RefinementResult(None, 'radial', plan.mesh_sha256, len(plan.source_geometry.volumes_m3), 'criteria_satisfied', False, ())
        source, destination = Path(source).resolve(), Path(destination).resolve()
        verify_sources(source, plan)
        if destination.exists() or destination.is_relative_to(source):
            raise ValueError('refinement needs a fresh destination outside the source')
        parent = read_radial_mesh(source)
        if not np.allclose(plan.source_geometry.centres_m, parent.geometry.centres_m, rtol=1e-8, atol=1e-12):
            raise ValueError('radial source geometry differs from plan')
        child = refine_mesh(parent, plan)
        destination.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=destination.name+'.stage-', dir=destination.parent))
        # materialize demands a fresh directory; its fixed child avoids any deletion.
        built = stage/'mesh'
        try:
            materialize(child, built)
            from plumetools.refinement.geometry import read_geometry
            g = read_geometry(built, time='0')
            index = np.searchsorted(parent.boundaries, g.centres_m[:,0], side='right')-1
            target = np.asarray(plan.target_edge_m).reshape(-1,parent.transverse**2).min(axis=1)[index]
            if np.any(g.max_edge_m > target*(1+1e-9)):
                raise ValueError('built radial mesh misses inherited target')
            metadata = {'backend': 'radial', 'target_status': 'inherited_targets_met',
                        'requires_flow_rerun': True, 'density_kind': 'parent-measured',
                        'parent_cells': parent.n_cells, 'new_cells': child.n_cells,
                        'parent_transverse': parent.transverse, 'new_transverse': child.transverse,
                        'source_sha256': plan.source_sha256, 'mesh_sha256': g.mesh_sha256}
            (built/'refinement-result.json').write_text(json.dumps(metadata, indent=2)+'\n')
            _publish_stage(built,destination)
            stage.rmdir()
            return RefinementResult(destination, 'radial', g.mesh_sha256, child.n_cells, 'inherited_targets_met', True, ())
        except Exception:
            # Leave failed staging/logs available to diagnose; source is untouched.
            raise
