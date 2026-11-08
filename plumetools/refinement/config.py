"""Standalone gas, field-role and sampling settings; no case-family dependency."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re

import yaml

from plumetools.nozzleplume.gas import VhsSpecies
from .model import Targets

ROLES = {"number_density", "temperature", "mean_parcels", "mass_density",
         "momentum_density", "kinetic_energy_density", "cell_centres", "cell_volumes"}
GAS_KEYS = {"species_name", "molecular_mass_kg", "molecular_diameter_m", "omega",
            "t_ref_K", "mean_free_path_convention", "internal_degrees_of_freedom"}


def _mapping(value, name):
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    return dict(value)


def _keys(raw, allowed, name):
    unknown = sorted(set(raw) - set(allowed))
    if unknown:
        raise ValueError(f"{name}: unknown keys {unknown}")


@dataclass(frozen=True)
class RefinementConfig:
    """Resolved portable reader settings; ``fields`` maps roles to filenames."""
    species: object
    convention: str
    profile: str
    fields: dict
    sampling: dict
    targets: Targets
    path: Path

    def as_dict(self):
        """Semantic settings without machine-specific source paths."""
        return {"gas": asdict(self.species), "convention": self.convention,
                "profile": self.profile, "fields": dict(self.fields),
                "sampling": dict(self.sampling), "targets": asdict(self.targets)}


def load_settings(path: Path) -> RefinementConfig:
    """Read a strict standalone YAML file for a measured DSMC case."""
    path = Path(path)
    try:
        raw = _mapping(yaml.safe_load(path.read_text(encoding="utf-8-sig")), str(path))
    except (OSError, yaml.YAMLError) as exc:
        raise ValueError(f"{path}: {exc}") from None
    _keys(raw, {"gas", "fields", "sampling", "targets"}, str(path))
    gas = _mapping(raw.get("gas"), "gas")
    _keys(gas, GAS_KEYS, "gas")
    required = {"molecular_mass_kg", "molecular_diameter_m"}
    if not required <= gas.keys():
        raise ValueError(f"gas: missing {sorted(required - gas.keys())}")
    convention = str(gas.get("mean_free_path_convention", "vhs"))
    if convention not in ("vhs", "hard_sphere"):
        raise ValueError("gas.mean_free_path_convention must be vhs or hard_sphere")
    if isinstance(gas.get("species_name", "gas"), (list, dict)):
        raise ValueError("single-species gas required; mixture MFP is not supported")
    species = VhsSpecies(name=str(gas.get("species_name", "gas")),
                         mass_kg=float(gas["molecular_mass_kg"]),
                         diameter_m=float(gas["molecular_diameter_m"]),
                         omega=float(gas.get("omega", 0.5)),
                         t_ref_K=float(gas.get("t_ref_K", 273)),
                         internal_degrees_of_freedom=int(gas.get("internal_degrees_of_freedom", 0)))
    fields = _mapping(raw.get("fields"), "fields")
    _keys(fields, ROLES | {"profile", "suffix"}, "fields")
    profile = str(fields.pop("profile", "standard"))
    suffix = str(fields.pop("suffix", species.name))
    if profile == "standard":
        defaults = {"number_density": "rhoNMean", "mass_density": "rhoMMean",
                    "momentum_density": "momentumMean", "kinetic_energy_density": "linearKEMean",
                    "mean_parcels": "dsmcRhoNMean"}
    elif profile == "hystrath":
        defaults = {"number_density": f"rhoN_{suffix}", "temperature": f"Ttra_{suffix}",
                    "mean_parcels": f"dsmcNMean_{suffix}"}
    elif profile == "custom":
        defaults = {}
    else:
        raise ValueError("fields.profile must be standard, hystrath or custom")
    resolved = {"cell_centres": "C", "cell_volumes": "V", **defaults, **fields}
    need = {"number_density", "mean_parcels"}
    if "temperature" not in resolved:
        need |= {"mass_density", "momentum_density", "kinetic_energy_density"}
    if not need <= resolved.keys():
        raise ValueError(f"fields: missing roles {sorted(need - resolved.keys())}")
    for role, name in resolved.items():
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.:-]*", name):
            raise ValueError(f"fields.{role}: expected a field filename, got {name!r}")
    sampling = _mapping(raw.get("sampling"), "sampling")
    _keys(sampling, {"start_s", "end_s", "basis", "source", "stationary", "n_steps"}, "sampling")
    targets = _mapping(raw.get("targets"), "targets")
    _keys(targets, {"edge_over_mfp", "min_mean_parcels"}, "targets")
    return RefinementConfig(species, convention, profile, resolved, sampling, Targets(**targets), path)
