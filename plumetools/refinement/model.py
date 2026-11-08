"""SI-unit contracts shared by refinement readers, decisions and backends."""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Optional, Protocol

import numpy as np


@dataclass(frozen=True)
class Geometry:
    """Actual centres, volumes and longest cell edges on an identified mesh."""
    centres_m: np.ndarray
    volumes_m3: np.ndarray
    max_edge_m: np.ndarray
    mesh_sha256: str

    def __post_init__(self):
        centres = np.asarray(self.centres_m, dtype=float)
        volumes = np.asarray(self.volumes_m3, dtype=float)
        edges = np.asarray(self.max_edge_m, dtype=float)
        n = len(volumes) if volumes.ndim == 1 else 0
        if not n or centres.shape != (n, 3) or edges.shape != (n,):
            raise ValueError("cell geometry requires centres (n,3), volumes (n), edges (n)")
        if not np.all(np.isfinite(centres)):
            raise ValueError("cell centres must be finite")
        if not np.all(np.isfinite(volumes) & (volumes > 0)) or not math.isfinite(float(volumes.sum())):
            raise ValueError("cell volumes and their sum must be finite and positive")
        if not np.all(np.isfinite(edges) & (edges > 0)):
            raise ValueError("cell edge lengths must be finite and positive")
        if not str(self.mesh_sha256):
            raise ValueError("cell geometry requires a mesh identity")
        for name, value in (("centres_m", centres), ("volumes_m3", volumes), ("max_edge_m", edges)):
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class Measurements:
    """Sampled state/counts; retain invalid values for unresolved-cell reporting."""
    number_density: np.ndarray
    temperature_K: np.ndarray
    mean_parcels: np.ndarray
    window: dict
    source_sha256: str
    provenance: dict = field(default_factory=dict)

    def __post_init__(self):
        n = np.asarray(self.number_density, dtype=float)
        if n.ndim != 1:
            raise ValueError("number_density must have one entry per cell")
        for name in ("number_density", "temperature_K", "mean_parcels"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.shape != n.shape:
                raise ValueError(f"{name} must have {n.size} cell entries")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class Targets:
    """Measured acceptance targets, separate from an initial mesh prediction."""
    edge_over_mfp: float = 0.8
    min_mean_parcels: float = 20.0

    def __post_init__(self):
        for name in ("edge_over_mfp", "min_mean_parcels"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class Assessment:
    """Statuses: -1 unresolved, 0 failing, 1 passing; summary covers every cell."""
    targets: Targets
    lambda_m: np.ndarray
    target_edge_m: np.ndarray
    mesh_status: np.ndarray
    occupancy_status: np.ndarray
    valid: np.ndarray
    summary: dict


@dataclass(frozen=True)
class RefinementPlan:
    """Parent-bound edge targets and an optional global parcel-weight estimate."""
    action: str
    mesh_sha256: str
    source_sha256: str
    source_geometry: Geometry
    target_edge_m: np.ndarray
    density_estimate: np.ndarray
    selected_cells: np.ndarray
    target_parcels: float
    recommended_weight: Optional[float]
    provenance: dict


@dataclass(frozen=True)
class RefinementResult:
    """Fresh mesh artifact; inherited targets are never child measurements."""
    destination: Optional[Path]
    backend: str
    mesh_sha256: str
    n_cells: int
    target_status: str
    requires_flow_rerun: bool
    logs: list


class MeshRefiner(Protocol):
    """A materializer leaves source results and solver initialization alone."""
    def refine(self, source: Path, plan: RefinementPlan, destination: Path) -> RefinementResult:
        """Materialize a new mesh, or return a no-op for a compliant plan."""
        ...
