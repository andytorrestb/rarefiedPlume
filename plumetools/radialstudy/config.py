"""Strict study controls, separate from the inherited physical reference."""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
from pathlib import Path

import yaml

from plumetools.cai2012.config import load_case_config
from plumetools.refinement.model import Targets


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                     separators=(',', ':')).encode()).hexdigest()


@dataclass(frozen=True)
class Domain:
    inlet_over_D: float | None = None
    outer_over_D: float = 25.
    angle_deg: float | None = None


@dataclass(frozen=True)
class Mesh:
    predicted_edge_over_mfp: float = 1.
    geometric_cap_over_D: float = .05
    max_cells: int = 1000000


@dataclass(frozen=True)
class Time:
    startup_transits: float = 2.
    sampling_transits: float = 6.
    blocks: int = 6
    windows_per_block: int = 2
    collision_ratio: float = .2
    transport_ratio: float = .2
    max_steps: int = 10000000


@dataclass(frozen=True)
class Run:
    cores: int = 1
    seed: int = 1729
    max_predicted_parcels: int = 10000000


@dataclass(frozen=True)
class Breakdown:
    definition: str = 'general'
    threshold: float = .05


@dataclass(frozen=True)
class Analysis:
    density_tolerance: float | None = None
    velocity_tolerance: float | None = None
    temperature_tolerance: float | None = None


@dataclass(frozen=True)
class RadialStudyConfig:
    path: Path
    reference: object
    reference_sha256: str
    physical_inputs: dict
    knudsen_numbers: tuple
    domain: Domain
    mesh: Mesh
    targets: Targets
    time: Time
    run: Run
    breakdown: Breakdown
    analysis: Analysis

    def settings(self):
        return {key: asdict(getattr(self, key)) for key in
                ('domain', 'mesh', 'targets', 'time', 'run', 'breakdown', 'analysis')}


def _section(cls, raw, name):
    if not isinstance(raw, dict):
        raise ValueError(f'{name} must be a mapping')
    unknown = set(raw) - {f.name for f in fields(cls)}
    if unknown:
        raise ValueError(f'{name}: unknown keys {sorted(unknown)}')
    return cls(**raw)


def load_study(path: Path) -> RadialStudyConfig:
    path = Path(path).resolve()
    raw = yaml.safe_load(path.read_text(encoding='utf-8-sig'))
    allowed = {'model', 'reference_case', 'knudsen_numbers', 'domain', 'mesh',
               'targets', 'time', 'run', 'breakdown', 'analysis'}
    if not isinstance(raw, dict) or set(raw) - allowed:
        raise ValueError('study: expected mapping with known study keys')
    if raw.get('model') != 'cai_exit_radial_resolution':
        raise ValueError('study.model must be cai_exit_radial_resolution')
    refpath = (path.parent / raw.get('reference_case', 'reference.yaml')).resolve()
    reference = load_case_config(refpath)
    if reference.exit.characteristic_length != 'diameter':
        raise ValueError('radial study requires exit Kn based on diameter')
    if reference.geometry.symmetry_mode != 'none' or reference.gas.internal_degrees_of_freedom != 0:
        raise ValueError('reference must be a full-domain monatomic gas')
    physical = {'model': 'cai2012_circular_plume', **{k: asdict(getattr(reference, k))
                for k in ('nozzle', 'exit', 'gas')}}
    cfg = RadialStudyConfig(path, reference, hashlib.sha256(refpath.read_bytes()).hexdigest(),
        physical, tuple(float(k) for k in raw.get('knudsen_numbers', [100, 10, 1, .1, .01, .001, .0001])),
        *[_section(cls, raw.get(key, {}), key) for key, cls in (
            ('domain', Domain), ('mesh', Mesh), ('targets', Targets), ('time', Time),
            ('run', Run), ('breakdown', Breakdown), ('analysis', Analysis))])
    if not cfg.knudsen_numbers or len(set(cfg.knudsen_numbers)) != len(cfg.knudsen_numbers):
        raise ValueError('knudsen_numbers must be nonempty and unique')
    positive = [*cfg.knudsen_numbers, cfg.domain.outer_over_D, cfg.mesh.predicted_edge_over_mfp,
                cfg.mesh.geometric_cap_over_D, cfg.time.startup_transits, cfg.time.sampling_transits,
                cfg.time.collision_ratio, cfg.time.transport_ratio, cfg.breakdown.threshold]
    positive += [x for x in (cfg.domain.inlet_over_D, cfg.domain.angle_deg,
                  *asdict(cfg.analysis).values()) if x is not None]
    if any(not math.isfinite(x) or x <= 0 for x in positive):
        raise ValueError('study values must be finite and positive')
    for x in (cfg.mesh.max_cells, cfg.time.blocks, cfg.time.windows_per_block,
              cfg.time.max_steps, cfg.run.cores, cfg.run.seed, cfg.run.max_predicted_parcels):
        if isinstance(x, bool) or not isinstance(x, int) or x <= 0:
            raise ValueError('counts, seed and budgets must be positive integers')
    if cfg.time.blocks < 2 or cfg.time.windows_per_block < 2:
        raise ValueError('sampling requires at least two blocks and two windows per block')
    if cfg.run.cores != 1:
        raise ValueError('this runner supports serial execution (cores: 1)')
    if cfg.domain.angle_deg is not None and cfg.domain.angle_deg >= 45:
        raise ValueError('angle_deg must be below 45 degrees')
    if cfg.breakdown.definition not in ('general', 'radial_approx'):
        raise ValueError('breakdown.definition must be general or radial_approx')
    return cfg
