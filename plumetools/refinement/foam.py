"""Read normalized standard, hyStrath or explicitly mapped DSMC fields."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Optional

from plumetools.foamio.times import latest_time
from plumetools.mesh._parse import strip_comments
from plumetools.rarefaction.core import translational_temperature
from plumetools.rarefaction.foam import read_window, weight_check
from .config import RefinementConfig
from .geometry import field_values, file_sha256, mesh_files, read_geometry
from .model import Geometry, Measurements
from .assess import usable_window

FIELD_DIMS = {
    "number_density": (0, -3, 0, 0, 0, 0, 0), "temperature": (0, 0, 0, 1, 0, 0, 0),
    "mean_parcels": (0, 0, 0, 0, 0, 0, 0), "mass_density": (1, -3, 0, 0, 0, 0, 0),
    "momentum_density": (1, -2, -1, 0, 0, 0, 0), "kinetic_energy_density": (1, -1, -2, 0, 0, 0, 0),
}


@dataclass(frozen=True)
class Snapshot:
    """OpenFOAM I/O wrapper; the normalized API does not require this type."""
    geometry: Geometry
    measured: Measurements
    species: object
    convention: str
    particle_weight: Optional[float]
    time: str
    provenance: dict


def read_snapshot(case: Path, config: RefinementConfig, *, time: Optional[str] = None) -> Snapshot:
    """Read one reconstructed static case, without writing into it."""
    case = Path(case)
    time = time or latest_time(case)
    if not (case / time).is_dir():
        raise ValueError(f"missing time {time} in {case}")
    geometry = read_geometry(case, time=time, field_names=config.fields)
    n_cells = len(geometry.volumes_m3)
    roles = ["number_density", "mean_parcels"]
    roles += (["temperature"] if "temperature" in config.fields else
              ["mass_density", "momentum_density", "kinetic_energy_density"])
    arrays, sources = {}, mesh_files(case)
    for role in roles:
        path = case / time / config.fields[role]
        arrays[role] = field_values(path, n_cells, FIELD_DIMS[role], vector=(role == "momentum_density"))
        sources[path.relative_to(case).as_posix()] = file_sha256(path)
    for role in ("cell_centres", "cell_volumes"):
        path = case / time / config.fields[role]
        sources[path.relative_to(case).as_posix()] = file_sha256(path)
    if "temperature" not in arrays:
        arrays["temperature"] = translational_temperature(
            arrays["mass_density"], arrays["momentum_density"], arrays["kinetic_energy_density"],
            config.species.gas_constant_j_per_kg_k)
    window = dict(config.sampling)
    if config.profile == "standard":
        names = [config.fields[role] for role in roles]
        bases = list(dict.fromkeys(name[:-4] if name.endswith("Mean") else name for name in names))
        known = [name for name in bases if read_window(case,time,[name])[0] is not None]
        native, _ = read_window(case, time, known) if known else (None, "")
        state = Path(time) / "uniform/functionObjects/functionObjectProperties"
        if (case/state).is_file():
            sources[state.as_posix()] = file_sha256(case/state)
        if native is not None:
            for key, actual in (("start_s",native.start_s),("end_s",native.end_s),("n_steps",native.n_steps)):
                if key in window and not math.isclose(float(window[key]),actual,rel_tol=1e-9,abs_tol=1e-12):
                    raise ValueError(f"declared sampling window {key} disagrees with available native field metadata")
        if native is not None and len(known)==len(bases):
            window = {**window,"start_s": native.start_s, "end_s": native.end_s,
                      "basis": "mean_per_cell", "source": "verified", "n_steps": native.n_steps}
        elif usable_window(window):
            # A declaration covers every consumed role; partial native metadata
            # may corroborate it but can never verify the entire snapshot.
            window = {**window,"source":"declared","native_metadata_missing_fields":sorted(set(bases)-set(known))}
        else:
            window = {"reason":"complete averaging metadata or an explicit sampling declaration is required",
                      "native_metadata_missing_fields":sorted(set(bases)-set(known))}
    if window and "end_s" in window and abs(float(window["end_s"]) - float(time)) > 1e-9 * max(1., abs(float(time))):
        raise ValueError("sampling window end does not match the selected field time")
    props = case / "constant/dsmcProperties"
    weight, weight_report = None, {}
    if props.is_file():
        text = strip_comments(props.read_text(encoding="utf-8", errors="replace"))
        ids = re.search(r"\btypeIdList\s*\(([^)]*)\)", text)
        if ids and len(ids.group(1).split()) > 1:
            raise ValueError("single-species MFP reader cannot assess a mixture case")
        match = re.search(r"^\s*nEquivalentParticles\s+([^;\s]+)\s*;", text, re.M)
        variable = re.search(r"\bcoordinateSystem\s+(?!dsmcCartesian\b)\w+", text)
        if match and not variable:
            candidate = float(match.group(1))
            if not math.isfinite(candidate) or candidate <= 0:
                raise ValueError("nEquivalentParticles must be finite and positive")
            weight_report = weight_check(arrays["number_density"], geometry.volumes_m3,
                                         arrays["mean_parcels"], candidate)
            if weight_report.get("consistent"):
                weight = candidate
        sources[props.relative_to(case).as_posix()] = file_sha256(props)
    try:
        relative_config = config.path.resolve().relative_to(case.resolve()).as_posix()
    except ValueError:
        relative_config = None
    if relative_config and config.path.is_file():
        sources[relative_config] = file_sha256(config.path)
    provenance = {"source_files": sources, "settings": config.as_dict(), "time": time,
                  "window": window, "particle_weight_check": weight_report}
    digest = hashlib.sha256(json.dumps(provenance, sort_keys=True, allow_nan=False).encode()).hexdigest()
    measured = Measurements(arrays["number_density"], arrays["temperature"], arrays["mean_parcels"],
                            window, digest, provenance)
    return Snapshot(geometry, measured, config.species, config.convention, weight, time, provenance)
