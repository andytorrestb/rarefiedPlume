"""Structured configuration for one hybrid case: load, overlay a profile, validate.

A case is described by one YAML file (``hybrid.yaml``). It has a ``base``
section and named ``profiles`` that overlay it -- ``smoke``, ``verification``,
``production`` by convention. Selecting a profile deep-merges it over ``base``;
the merged document is what every step reads and what the provenance records.

Nothing here derives physics: :func:`build_source` hands the physical inputs to
:mod:`plumetools.hybridplume.physics`, and every derived number comes from there.

Values are labelled where they come from (paper, repository, assumption) in the
YAML itself; :data:`PAPER_VALUES` lists the ones Cai & Wang (2012) state, so the
report can show which inputs a configuration has moved away from them.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import yaml

from .physics import Gas, SourceState

SCHEMA = "hybridplume/1"

#: What Cai & Wang (2012) state. Everything else is a repository value or an
#: assumption, and is labelled as such in the case YAML.
PAPER_VALUES = {
    "source.diameter_m": 0.2,
    "source.speed_ratio": 2.0,
    "gas.name": "Ar",
    "knudsen_cases": [100.0, 0.1, 0.01],
}

BACKENDS = ("dsmcFoam",)
INTERFACE_MODES = ("automatic", "manual", "none", "fitted-ellipse", "fitted-ellipsoid")
MANUAL_SURFACES = ("sphere", "cylinder", "box")


class ConfigError(ValueError):
    """A configuration that cannot describe a valid case."""


def deep_merge(base: dict, overlay: dict) -> dict:
    out = copy.deepcopy(base)
    for key, value in (overlay or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _require(doc: dict, dotted: str):
    node = doc
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            raise ConfigError(f"missing required key {dotted!r}")
        node = node[part]
    return node


@dataclass(frozen=True)
class CaseConfig:
    """The merged document plus where it came from."""

    doc: dict
    path: Path | None
    profile: str

    def get(self, dotted: str, default=None):
        node = self.doc
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def require(self, dotted: str):
        return _require(self.doc, dotted)

    @property
    def digest(self) -> str:
        """SHA-256 of the merged document, canonically serialised."""
        text = json.dumps(self.doc, sort_keys=True, default=str)
        return hashlib.sha256(text.encode()).hexdigest()


def load(path: Path, profile: str | None = None, overrides: dict | None = None) -> CaseConfig:
    """Load ``hybrid.yaml``, overlay ``profile`` then ``overrides``, validate."""
    path = Path(path)
    raw = yaml.safe_load(path.read_text())
    return from_document(raw, profile=profile, overrides=overrides, path=path)


def from_document(raw: dict, *, profile: str | None = None,
                  overrides: dict | None = None, path: Path | None = None) -> CaseConfig:
    if not isinstance(raw, dict) or raw.get("schema") != SCHEMA:
        raise ConfigError(f"expected schema {SCHEMA!r}, got {raw.get('schema') if isinstance(raw, dict) else raw!r}")
    base = raw.get("base")
    if not isinstance(base, dict):
        raise ConfigError("no 'base' section")
    profiles = raw.get("profiles", {}) or {}
    name = profile or raw.get("default_profile") or "base"
    if name != "base":
        if name not in profiles:
            raise ConfigError(f"unknown profile {name!r}; known: {sorted(profiles)}")
        doc = deep_merge(base, profiles[name])
    else:
        doc = copy.deepcopy(base)
    if overrides:
        doc = deep_merge(doc, overrides)
    cfg = CaseConfig(doc=doc, path=path, profile=name)
    validate(cfg)
    return cfg


def validate(cfg: CaseConfig) -> None:
    """Refuse configurations that would build a case other than the one described."""
    for key in ("source.diameter_m", "source.speed_ratio", "source.temperature_K",
                "source.knudsen", "gas.mass_kg", "gas.diameter_m", "gas.t_ref_K",
                "gas.omega", "domain.x_max_over_D", "domain.half_width_over_D",
                "mesh.core_cell_over_D", "dsmc.backend"):
        cfg.require(key)
    if cfg.get("dsmc.backend") not in BACKENDS:
        raise ConfigError(f"dsmc.backend {cfg.get('dsmc.backend')!r} not in {BACKENDS}")
    if cfg.get("interface.mode", "automatic") not in INTERFACE_MODES:
        raise ConfigError(f"interface.mode must be one of {INTERFACE_MODES}")
    if cfg.get("interface.mode") == "manual":
        kind = cfg.get("interface.manual.type")
        if kind not in MANUAL_SURFACES:
            raise ConfigError(f"interface.manual.type must be one of {MANUAL_SURFACES}")
    geometry = cfg.get("mesh.geometry", "cartesian")
    if geometry not in ("cartesian", "planar"):
        raise ConfigError("mesh.geometry must be cartesian or planar")
    fitted = cfg.get("interface.mode") == "fitted-ellipse"
    if fitted != (geometry == "planar"):
        raise ConfigError("planar geometry requires the fitted-ellipse interface mode")
    if cfg.get("interface.mode") == "fitted-ellipsoid" and geometry != "cartesian":
        raise ConfigError("fitted-ellipsoid requires cartesian 3D geometry")
    if geometry == "planar":
        depth = float(cfg.get("mesh.depth_over_D", 0))
        clearance = float(cfg.get("interface.fit.clearance_cells", 0.25))
        if not math.isfinite(depth) or depth <= 0:
            raise ConfigError("mesh.depth_over_D must be positive and finite")
        if not math.isfinite(clearance) or clearance < 0:
            raise ConfigError("interface.fit.clearance_cells must be finite and nonnegative")
    # Axisymmetry and radial weighting remain unsupported.
    for key in ("mesh.axisymmetric", "mesh.wedge", "dsmc.axisymmetric",
                "dsmc.radial_weighting", "domain.symmetry"):
        if cfg.get(key):
            raise ConfigError(f"{key} is unsupported: use explicit planar geometry for a slit prototype; "
                              "axisymmetry, wedges, symmetry planes and radial weights are unavailable")
    half = float(cfg.require("domain.half_width_over_D"))
    if half <= 0.5:
        raise ConfigError("domain.half_width_over_D must exceed 0.5 so the source disk fits")
    core_half = cfg.get("mesh.core_half_over_D")
    if core_half is not None and float(core_half) <= 0.5:
        raise ConfigError("mesh.core_half_over_D must exceed 0.5: every source face "
                          "must lie in the uniform core")
    bg = float(cfg.get("cfd.background_pressure_ratio", 1e-4))
    if not 0.0 < bg < 0.1:
        raise ConfigError("cfd.background_pressure_ratio must be in (0, 0.1): it is a "
                          "numerical floor, not an ambient gas")


def build_gas(cfg: CaseConfig) -> Gas:
    g = cfg.require("gas")
    return Gas(name=str(g.get("name", "Ar")), mass_kg=float(g["mass_kg"]),
               diameter_m=float(g["diameter_m"]), t_ref_K=float(g["t_ref_K"]),
               omega=float(g["omega"]), prandtl=float(g.get("prandtl", 2.0 / 3.0)),
               provenance=dict(g.get("provenance", {}) or {}))


def build_source(cfg: CaseConfig) -> SourceState:
    s = cfg.require("source")
    return SourceState(gas=build_gas(cfg), diameter_m=float(s["diameter_m"]),
                       speed_ratio=float(s["speed_ratio"]),
                       temperature_K=float(s["temperature_K"]),
                       knudsen=float(s["knudsen"]),
                       mean_free_path_convention=str(
                           s.get("mean_free_path_convention", "vhs")),
                       depth_m=(float(cfg.require("mesh.depth_over_D")) * float(s["diameter_m"])
                                if cfg.get("mesh.geometry") == "planar" else None))


def departures_from_paper(cfg: CaseConfig) -> dict:
    """Paper-stated inputs this configuration does not use, for the report."""
    out = {}
    if cfg.get("mesh.geometry") == "planar":
        out["source.geometry"] = {"value": "planar slit", "paper": "circular disk"}
    for key, value in PAPER_VALUES.items():
        if key == "knudsen_cases":
            kn = float(cfg.require("source.knudsen"))
            if kn not in value:
                out["source.knudsen"] = {"value": kn, "paper_cases": value,
                                         "label": "extension (not a Cai case)"}
            continue
        have = cfg.get(key)
        if have != value:
            out[key] = {"value": have, "paper": value}
    return out
