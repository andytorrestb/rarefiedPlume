"""``unigas.yaml``: the uniGasFoam settings of one cai2012-unigas case.

The cai2012 ``case.yaml`` schema rejects new sections, so the solver settings
live beside it, as interface mode keeps ``interface.yaml``. Every key has the
value the approved design fixes (spec §4); a case file only needs the keys it
changes, and the family's template states them all explicitly.

Two of the owner's decisions are refusals rather than options:

* decision 4 -- DSMC cells respect the collision-time limit
  (``adaptive.dsmc_cells_time_step_limited`` may not be false);
* decision 5 -- cell weighting is off (``cell_weighted`` may not be true).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path

import yaml

FILENAME = "unigas.yaml"

USP_MODEL = "unifiedStochasticParticleSBGK"


class UniGasConfigError(ValueError):
    """A ``unigas.yaml`` the family must not run."""


@dataclass(frozen=True)
class CollisionSettings:
    model: str = "hybrid"
    bgk_model: str = USP_MODEL
    dsmc_partner_model: str = "noTimeCounter"
    dsmc_model: str = "variableHardSphere"
    macro_interpolation: bool = True
    theta: float = 1.0e-3


@dataclass(frozen=True)
class AdaptiveSettings:
    adaptation_interval: int = 20
    max_sub_cell_size_mfp_ratio: float = 0.5
    max_time_step_mct_ratio: float = 0.2
    max_courant_number: float = 0.5
    dsmc_cells_time_step_limited: bool = True
    max_sub_cell_levels: int = 10


@dataclass(frozen=True)
class DecompositionSettings:
    breakdown_max: float = 0.05
    theta: float = 0.2
    smoothing_passes: int = 5
    decomposition_interval: int = 100


@dataclass(frozen=True)
class WriteSettings:
    startup_interval_transits: float = 1.0
    keep_lagrangian: bool = False


@dataclass(frozen=True)
class FieldSettings:
    mean_field: str = "Ar"
    window_field: str = "ArWindow"


@dataclass(frozen=True)
class VerdictSettings:
    max_effective_spacing_over_mfp: float | None = None   # None: checks.max_cell_over_mfp
    max_time_step_mct_ratio: float = 0.2
    min_particles_per_dsmc_cell: float | None = None       # None: checks.min_particles_per_cell
    min_particles_per_sub_cell: float = 2.0
    occupancy_floor_n_over_n0: float = 1.0e-3


@dataclass(frozen=True)
class UniGasSettings:
    cell_weighted: bool = False
    collision: CollisionSettings = field(default_factory=CollisionSettings)
    adaptive: AdaptiveSettings = field(default_factory=AdaptiveSettings)
    decomposition: DecompositionSettings = field(default_factory=DecompositionSettings)
    writes: WriteSettings = field(default_factory=WriteSettings)
    fields: FieldSettings = field(default_factory=FieldSettings)
    verdict: VerdictSettings = field(default_factory=VerdictSettings)


_SECTIONS = {
    "collision": CollisionSettings,
    "adaptive": AdaptiveSettings,
    "decomposition": DecompositionSettings,
    "writes": WriteSettings,
    "fields": FieldSettings,
    "verdict": VerdictSettings,
}


def _section(cls, data, name: str, where: str):
    if data is None:
        return cls()
    if not isinstance(data, dict):
        raise UniGasConfigError(f"{where}: '{name}' must be a mapping")
    known = {f.name: f for f in dataclasses.fields(cls)}
    unknown = sorted(set(data) - set(known))
    if unknown:
        raise UniGasConfigError(
            f"{where}: unknown key(s) {', '.join(f'{name}.{k}' for k in unknown)}")
    values = {}
    for key, value in data.items():
        default = getattr(cls(), key)
        values[key] = _coerce(value, default, f"{name}.{key}", where)
    return cls(**values)


def _coerce(value, default, key: str, where: str):
    if isinstance(default, bool):
        if not isinstance(value, bool):
            raise UniGasConfigError(f"{where}: {key} must be true or false")
        return value
    if isinstance(default, int):
        if isinstance(value, bool) or not isinstance(value, int):
            raise UniGasConfigError(f"{where}: {key} must be an integer")
        return value
    if isinstance(default, float) or default is None:
        if value is None and default is None:
            return None
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise UniGasConfigError(f"{where}: {key} must be a number")
        return float(value)
    if isinstance(default, str):
        if not isinstance(value, str) or not value:
            raise UniGasConfigError(f"{where}: {key} must be a non-empty string")
        return value
    raise UniGasConfigError(f"{where}: {key} has an unsupported type")


def _positive(value, key: str, where: str) -> None:
    if value is not None and not value > 0:
        raise UniGasConfigError(f"{where}: {key} must be positive, got {value}")


def _validate(s: UniGasSettings, where: str) -> None:
    if s.cell_weighted:
        raise UniGasConfigError(
            f"{where}: cell_weighted must be false -- the owner's decision 5 "
            f"(spec 'Decisions made by the owner'): one particle weight, as in the "
            f"dsmcFoam suite")
    if s.collision.model != "hybrid":
        raise UniGasConfigError(
            f"{where}: collision.model must be 'hybrid' (USP-DSMC), got "
            f"'{s.collision.model}'")
    if s.collision.bgk_model != USP_MODEL:
        raise UniGasConfigError(
            f"{where}: collision.bgk_model must be '{USP_MODEL}', got "
            f"'{s.collision.bgk_model}'")
    if not s.adaptive.dsmc_cells_time_step_limited:
        raise UniGasConfigError(
            f"{where}: adaptive.dsmc_cells_time_step_limited must be true -- the "
            f"owner's decision 4: DSMC cells respect the collision-time limit")
    if s.adaptive.max_sub_cell_levels < 1:
        raise UniGasConfigError(f"{where}: adaptive.max_sub_cell_levels must be >= 1")
    for name in ("adaptation_interval", "max_sub_cell_size_mfp_ratio",
                 "max_time_step_mct_ratio", "max_courant_number"):
        _positive(getattr(s.adaptive, name), f"adaptive.{name}", where)
    for name in ("breakdown_max", "theta", "decomposition_interval"):
        _positive(getattr(s.decomposition, name), f"decomposition.{name}", where)
    if s.decomposition.smoothing_passes < 0:
        raise UniGasConfigError(f"{where}: decomposition.smoothing_passes must be >= 0")
    _positive(s.collision.theta, "collision.theta", where)
    _positive(s.writes.startup_interval_transits, "writes.startup_interval_transits", where)
    for name in ("max_effective_spacing_over_mfp", "max_time_step_mct_ratio",
                 "min_particles_per_dsmc_cell", "min_particles_per_sub_cell",
                 "occupancy_floor_n_over_n0"):
        _positive(getattr(s.verdict, name), f"verdict.{name}", where)
    if s.fields.mean_field == s.fields.window_field:
        raise UniGasConfigError(
            f"{where}: fields.mean_field and fields.window_field must differ "
            f"(both '{s.fields.mean_field}')")


def settings_from_dict(data: dict, where: str) -> UniGasSettings:
    """Build and validate settings from a parsed ``unigas.yaml`` mapping."""
    data = {} if data is None else data
    if not isinstance(data, dict):
        raise UniGasConfigError(f"{where}: must be a mapping")
    unknown = sorted(set(data) - set(_SECTIONS) - {"cell_weighted"})
    if unknown:
        raise UniGasConfigError(f"{where}: unknown key(s) {', '.join(unknown)}")
    cell_weighted = data.get("cell_weighted", False)
    if not isinstance(cell_weighted, bool):
        raise UniGasConfigError(f"{where}: cell_weighted must be true or false")
    settings = UniGasSettings(
        cell_weighted=cell_weighted,
        **{name: _section(cls, data.get(name), name, where)
           for name, cls in _SECTIONS.items()})
    _validate(settings, where)
    return settings


def load_settings(case_dir: Path) -> UniGasSettings:
    """Read and validate ``<case_dir>/unigas.yaml``."""
    path = Path(case_dir) / FILENAME
    if not path.is_file():
        raise UniGasConfigError(f"no {path}: every cai2012-unigas case has a {FILENAME}")
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise UniGasConfigError(f"{path}: not valid YAML: {exc}") from exc
    return settings_from_dict(data, str(path))


def as_dict(settings: UniGasSettings) -> dict:
    """The settings as plain data, for provenance."""
    return dataclasses.asdict(settings)
