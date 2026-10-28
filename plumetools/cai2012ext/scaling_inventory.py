"""Read-only inventory and opt-in registration for the whole-field experiment.

Registration is deliberately separate from production generation: it authorizes
analysis of explicitly named, study-local results and writes input specifications
only.  It does not bypass any legacy generation, meshing or execution guard.
Neither a profile CSV nor an image is a whole-field estimator with cell volumes.
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
import yaml

from plumetools.cai2012 import checks, inflow, mesh, post
from plumetools.cai2012.config import load_case_config, validate
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012.study import load_study
from plumetools.cai2012ext import feasibility, interface, policy, provenance


TARGET_DECADES = (100.0, 10.0, 1.0, 0.1, 0.01, 0.001, 0.0001)
REGISTRATION_SCOPE = "cai2012-kn-extension-whole-field-scaling-v1"
DEFAULT_MAX_BYTES = 512 * 1024 ** 2
BLOCK_METADATA = ("independent_blocks", "temporal_correlation_method", "block_duration_s",
                  "integrated_autocorrelation_time_s")
SUITABILITY_CRITERIA = {
    "pre_run": "Every plumetools.cai2012ext.feasibility criterion passes.",
    "mesh": "Measured bounding-box cell/mfp maximum <= the configured diagnostic target.",
    "occupancy": "The measured plume region has no cells below checks.min_particles_per_cell.",
    "stationarity": "Allconverge assessment verdict is stable-window-demonstrated.",
    "particle_weight": "The measured occupancy weight_check is consistent.",
    "coverage": "Actual matched whole-domain density and cell volumes are available.",
    "source": "No acknowledged unresolved source-boundary limitation applies to this direct-exit case.",
    "unknown": "Absent, inconclusive or unverifiable evidence does not pass.",
    "uncertainty": "Occupancy and cumulative running means are not independent blocks or confidence intervals.",
}


class InventoryError(ValueError):
    """A registration or field source is incompatible with this experiment."""


def _read(path):
    path = Path(path)
    if not path.is_file():
        return {}
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        raise InventoryError(f"Cannot read {path}: {exc}") from exc


def _local(root, value):
    path = (Path(root) / str(value)).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise InventoryError(f"Registration path must remain inside this study: {value}")
    return path


def read_registration(study_root, path=None):
    """Read an explicit analysis opt-in; no registration is implicit by default."""
    if path is None:
        return {"cases": [], "prepared_decades": []}
    root = Path(study_root).resolve()
    source = _local(root, path)
    document = _read(source)
    if document.get("scope") != REGISTRATION_SCOPE:
        raise InventoryError(f"{source}: expected scope {REGISTRATION_SCOPE!r}")
    if document.get("purpose") != "analysis-and-specification-only":
        raise InventoryError("Registration never authorizes production generation or execution")
    if document.get("enabled") is not True:
        raise InventoryError("Registration requires enabled: true as an explicit opt-in")
    names, paths = set(), set()
    for entry in document.get("cases") or []:
        if not {"name", "path", "Kn"} <= entry.keys():
            raise InventoryError("Each registration needs name, path and Kn")
        _local(root, entry["path"])
        resolved = str(_local(root, entry["path"]))
        if entry["name"] in names or resolved in paths:
            raise InventoryError("Duplicate case name or path in scaling registration")
        names.add(entry["name"])
        paths.add(resolved)
        if entry.get("field_npz"):
            _local(root, entry["field_npz"])
        if not any(policy.same_knudsen(entry["Kn"], k) for k in TARGET_DECADES):
            raise InventoryError(f"{entry['name']}: Kn is outside the seven-decade experiment")
        if policy.is_excluded(entry["Kn"]) and document.get("allow_paper_decades") is not True:
            raise InventoryError("Paper-decade analysis requires allow_paper_decades: true")
    for kn in document.get("prepared_decades") or []:
        if not any(policy.same_knudsen(kn, k) for k in TARGET_DECADES):
            raise InventoryError(f"Cannot prepare Kn={kn}: outside target decades")
        if policy.is_excluded(kn) and document.get("allow_paper_decades") is not True:
            raise InventoryError("Paper-decade specifications require allow_paper_decades: true")
    return document


def _replace_config(cfg, overrides):
    for section, values in (overrides or {}).items():
        if section not in policy.SECTIONS or not isinstance(values, dict):
            raise InventoryError(f"Invalid specification override section {section!r}")
        try:
            cfg = dataclasses.replace(cfg, **{
                section: dataclasses.replace(getattr(cfg, section), **values)})
        except TypeError as exc:
            raise InventoryError(str(exc)) from exc
    validate(cfg)
    return cfg


def _assess_config(cfg):
    """Use existing bounded analytical planners, never write mesh dictionaries."""
    geom = from_config(cfg)
    state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, state)
    run = inflow.derive_run_settings(cfg, state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    particles = mesh.estimate_particles(plan, geom, state, run.n_equivalent_particles)
    check = checks.run_checks(cfg, geom, state, plan, run,
                             nozzle=nozzle, particles=particles)
    assessment = feasibility.assess(cfg, geom, state, plan, run, nozzle,
                                    particles, check)
    return provenance.exit_conditions(cfg, geom, state), assessment


def _stratum(cfg):
    # Numerical mesh and time-step choices may differ; the source problem may not.
    value = {"D_m": cfg.nozzle.diameter_m, "S0": cfg.exit.speed_ratio,
             "T0_K": cfg.exit.T0_K, "gas": provenance.gas_parameters(cfg),
             "geometry": dataclasses.asdict(cfg.geometry),
             "characteristic_length": cfg.exit.characteristic_length}
    encoded = json.dumps(value, sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _npz_description(path, max_bytes=DEFAULT_MAX_BYTES):
    # ZIP metadata bounds decompression before numpy reads any payload.
    with zipfile.ZipFile(path) as archive:
        expanded = sum(item.file_size for item in archive.infolist())
        if expanded > max_bytes:
            raise InventoryError(f"{path}: {expanded} uncompressed bytes exceed limit {max_bytes}")
        fields = sorted(Path(item.filename).stem for item in archive.infolist())
    required = {"points", "density", "volumes"}
    if not required <= set(fields):
        raise InventoryError(f"{path}: whole-field NPZ needs {sorted(required)}")
    return fields, expanded


def _field_source(case_dir, registration, root):
    times = []
    if case_dir.is_dir():
        for path in case_dir.iterdir():
            if path.is_dir():
                try:
                    if float(path.name) > 0:
                        times.append(path.name)
                except ValueError:
                    pass
    times.sort(key=float)
    partial = sorted(p.name for p in (case_dir / "results").glob("*.csv"))
    field = registration.get("field_npz")
    if field:
        path = _local(root, field)
        if path.is_file():
            names, size = _npz_description(path)
            return {"kind": "npz", "path": str(path), "bytes_uncompressed": size}, names, times, partial
    geometry = next((t for t in reversed(times)
                     if (case_dir / t / "C").is_file() and (case_dir / t / "V").is_file()), None)
    selected = registration.get("time")
    candidates = [str(selected)] if selected is not None else list(reversed(times))
    for time in candidates:
        directory = case_dir / time
        if (directory / "rhoNMean").is_file() and geometry:
            names = sorted(p.name for p in directory.iterdir() if p.is_file())
            return {"kind": "openfoam", "path": str(case_dir), "time": time,
                    "mesh_time": geometry}, names, times, partial
    return None, [], times, partial


def reference_suitability(assessment, diagnostics, occupancy, stationarity, fit_ready):
    """Conservative, explicit quality verdict; unknown evidence never passes."""
    reasons = []
    if assessment.get("production_ready") is not True:
        reasons.extend(assessment.get("unresolved") or ["Pre-run feasibility is unknown"])
    resolution = (diagnostics.get("resolution") or {}).get("bounding_box_max") or {}
    target = (diagnostics.get("resolution") or {}).get("target_cell_over_mean_free_path")
    maximum = resolution.get("max")
    if target is None or maximum is None:
        reasons.append("Measured mesh/mfp diagnostics unavailable")
    elif float(maximum) > float(target):
        reasons.append("Measured maximum cell/mfp exceeds the diagnostic target")
    plume = next((r for r in occupancy.get("regions", []) if r.get("name") == "plume"), {})
    below = plume.get("fraction_below_floor")
    if below is None:
        reasons.append("Measured plume occupancy unavailable")
    elif float(below) > 0:
        reasons.append(f"Plume cells below configured parcel floor: {float(below):.6g}")
    if (occupancy.get("weight_check") or {}).get("consistent") is not True:
        reasons.append("Particle weight consistency unavailable or failed")
    verdict = (stationarity.get("assessment") or {}).get("verdict")
    if verdict != "stable-window-demonstrated":
        reasons.append(f"Stationarity not demonstrated ({verdict or 'unknown'})")
    if not fit_ready:
        reasons.append("Matched whole-field density and cell volumes unavailable")
    return {"suitable": not reasons, "status": "adequate" if not reasons else "questionable-or-unknown",
            "reasons": reasons, "criteria": SUITABILITY_CRITERIA}


def _case_row(root, entry, fallback_cfg=None):
    case_dir = _local(root, entry["path"])
    row = {"name": entry["name"], "path": str(case_dir), "Kn": float(entry["Kn"]),
           "registration": entry.get("registration", "study.yaml"), "warnings": []}
    try:
        exists = (case_dir / "case.yaml").is_file()
        cfg = load_case_config(case_dir) if exists else fallback_cfg
    except (OSError, ValueError) as exc:
        exists, cfg = False, None
        row["warnings"].append(f"Configuration unavailable: {exc}")
    row["configuration_exists"] = exists
    if cfg is None:
        row.update({"verified_kn": None, "compatible": False, "stratum": None,
                    "suitable": False, "fit_ready": False, "result_status": "missing-configuration",
                    "source_type": entry.get("declared_source_type", "unknown"), "field_source": None,
                    "available_fields": [], "times": []})
        return row
    source_type = "direct-exit"
    if exists:
        try:
            if interface.load_interface(case_dir) is not None:
                source_type = "downstream-interface"
        except (OSError, ValueError) as exc:
            source_type = "unknown"
            row["warnings"].append(f"Source configuration unavailable: {exc}")
    if source_type == "direct-exit":
        exit_values, assessment = _assess_config(cfg)
    else:
        exit_values = provenance.exit_conditions(cfg, from_config(cfg), inflow.from_config(cfg))
        assessment = {"status": "separate-source-problem", "production_ready": False,
                      "unresolved": ["Direct-exit feasibility criteria are not applicable to interface-source fields"]}
    row.update({"exit_state": exit_values, "verified_kn": exit_values["Kn_recomputed"],
                "gas": provenance.gas_parameters(cfg), "geometry": dataclasses.asdict(cfg.geometry),
                "stratum": _stratum(cfg), "speed_ratio": float(cfg.exit.speed_ratio),
                "feasibility": assessment, "configuration_sha256": provenance.sha256(case_dir / "case.yaml")})
    known_limitations = list(cfg.meta.get("known_limitations") or [])
    known_limitations.extend(entry.get("known_source_limitations") or [])
    known_limitations.extend(policy.source_limitations(cfg))
    if source_type == "direct-exit" and float(cfg.exit.knudsen) <= .001 and not known_limitations:
        known_limitations.append(
            "Known low-Kn limitation: direct-exit inflow surface is not the intended collisional source surface; "
            "source-boundary correction is deferred. Results remain exploratory.")
    row["known_source_limitations"] = list(dict.fromkeys(known_limitations))
    row["knudsen_suite"] = cfg.meta.get("knudsen_suite", "legacy")
    matches = policy.same_knudsen(entry["Kn"], cfg.exit.knudsen)
    row["source_type"] = source_type
    row["compatible"] = matches and source_type == "direct-exit"
    if not matches:
        row["warnings"].append("Registered and configured Kn disagree; excluded")
    if source_type != "direct-exit":
        row["warnings"].append("Interface injection changes the source problem; excluded")
    result = case_dir / "results"
    evidence = []
    for filename in ("diagnostics_summary.yaml", "occupancy.yaml", "stationarity.yaml"):
        try:
            evidence.append(_read(result / filename))
        except (InventoryError, OSError) as exc:
            row["warnings"].append(str(exc))
            evidence.append({})
    diagnostics, occupancy, stationarity = evidence
    row.update({"diagnostics": diagnostics, "occupancy": occupancy, "stationarity": stationarity})
    try:
        field, fields, times, partial = _field_source(case_dir, entry, root)
    except (InventoryError, OSError, zipfile.BadZipFile) as exc:
        row["warnings"].append(str(exc))
        field, fields, times, partial = None, [], [], []
    row.update({"field_source": field, "available_fields": fields, "times": times,
                "partial_exports": partial, "fit_ready": field is not None and row["compatible"],
                "result_status": "whole-field-available" if field else "partial-exports-only" if partial else "missing-results"})
    quality = reference_suitability(assessment, diagnostics, occupancy, stationarity, row["fit_ready"])
    if row["known_source_limitations"]:
        quality["reasons"].extend(row["known_source_limitations"])
        quality.update(suitable=False, status="questionable-or-unknown")
    row.update({"suitable": quality["suitable"] and row["compatible"], "suitability": quality})
    row["warnings"].extend(quality["reasons"])
    marker = case_dir / "NON-PRODUCTION"
    try:
        if marker.is_file():
            row["warnings"].append(marker.read_text(encoding="utf-8", errors="replace").strip())
    except OSError as exc:
        row["warnings"].append(f"Cannot read NON-PRODUCTION marker: {exc}")
    block_available = field is not None and "density_blocks" in fields
    row["uncertainty"] = {
        "status": "blocks-present-require-correlation-analysis" if block_available else "statistical-unavailable",
        "time_blocks_available": block_available,
        "note": "Cumulative averaged snapshots and occupancy are not independent uncertainty estimates.",
        "sampling": diagnostics.get("sampling") or {},
        "declared_block_metadata": {key: entry.get("uncertainty", {}).get(key)
                                    for key in BLOCK_METADATA if key in entry.get("uncertainty", {})},
    }
    return row


def _configuration_specs(root, registration, study=None):
    study = study or load_study(root / "study.yaml")
    policy.check_study(study)
    template = load_case_config(root / study.base_case)
    specs = []
    for case in study.expand():
        cfg = dataclasses.replace(template, exit=dataclasses.replace(template.exit, knudsen=float(case.knudsen)))
        cfg = _replace_config(cfg, case.overrides)
        if study.meta.get("knudsen_suite"):
            cfg = dataclasses.replace(cfg, meta=dict(cfg.meta, knudsen_suite=study.meta["knudsen_suite"]))
        specs.append((dict(name=case.name, Kn=float(case.knudsen), path=f"{study.cases_dir}/{case.name}"), cfg))
    existing = [float(entry["Kn"]) for entry, _ in specs]
    for kn in registration.get("prepared_decades") or []:
        if any(policy.same_knudsen(kn, k) for k in existing):
            continue
        name = "Kn" + f"{float(kn):g}".replace(".", "p") + "_np1x"
        cfg = dataclasses.replace(template, exit=dataclasses.replace(template.exit, knudsen=float(kn)))
        specs.append((dict(name=name, Kn=float(kn), path=f"ScalingPrepared/{name}",
                          registration="explicit-analysis-specification"), cfg))
    return specs


def _external_equivalents(root):
    """Read configuration metadata only; never read or register external results."""
    found = []
    cases_root = root.parent
    for family in ("cai2012", "cai2012-health"):
        directory = cases_root / family
        study_path = directory / "study.yaml"
        if not study_path.is_file():
            continue
        document = _read(study_path)
        template_path = directory / str(document.get("base_case", "baseCase")) / "case.yaml"
        template = _read(template_path)
        values = [c.get("Kn") for c in document.get("kn_cases") or [] if c.get("Kn") is not None]
        if not values and template.get("exit", {}).get("knudsen") is not None:
            values = [template["exit"]["knudsen"]]
        for kn in sorted(set(float(k) for k in values)):
            if any(policy.same_knudsen(kn, target) for target in TARGET_DECADES):
                found.append({"study": family, "Kn": kn, "configuration": str(study_path),
                              "source": str(template_path), "included": False,
                              "status": "configured-equivalent-provenance-only; results not inspected"})
        # Input files only: do not recursively inspect result trees or large fields.
        patterns = ("Cases/*/case.yaml", "Cases/*/*/case.yaml")
        for pattern in patterns:
            for case_path in sorted(directory.glob(pattern)):
                candidate = _read(case_path)
                kn = candidate.get("exit", {}).get("knudsen")
                if kn is not None and any(policy.same_knudsen(kn, k) for k in TARGET_DECADES):
                    found.append({"study": family, "Kn": float(kn), "configuration": str(case_path),
                                  "included": False, "status": "existing-case-config-only; results not inspected"})
    return found


def inventory(study_root, registration_path=None):
    """Inventory local registrations, actual availability and all seven target Kn."""
    root = Path(study_root).resolve()
    registration = read_registration(root, registration_path)
    study = load_study(root / "study.yaml")
    specs = _configuration_specs(root, registration, study)
    manifest_path = _local(root, study.manifest_name)
    manifest = _read(manifest_path)
    problems = policy.check_manifest(manifest)
    if problems:
        raise InventoryError(f"{study.manifest_name} violates its production policy: " + "; ".join(problems))
    if manifest and policy.suite_name(manifest) != policy.suite_name(study):
        raise InventoryError(f"{study.manifest_name} profile disagrees with the active study")
    listed = {str(Path(e["path"]).as_posix()): e for e in manifest.get("cases") or []}
    rows = []
    paths = set()
    for entry, cfg in specs:
        merged = dict(entry)
        if entry["path"] in listed:
            merged.update(listed[entry["path"]])
            merged["registration"] = f"{study.manifest_name} + study.yaml"
        rows.append(_case_row(root, merged, cfg))
        paths.add(str(_local(root, entry["path"])))
    # Explicit registration is the only route from an archived or newly supplied
    # case into this analysis. Merely finding a directory is not registration.
    for entry in registration.get("cases") or []:
        merged = dict(entry, registration=str(registration_path))
        resolved = str(_local(root, entry["path"]))
        if resolved in paths:
            index = next(i for i, row in enumerate(rows) if row["path"] == resolved)
            rows[index] = _case_row(root, merged)
        else:
            rows.append(_case_row(root, merged))
            paths.add(resolved)
    # Generated direct and interface manifests register additional local cases;
    # interface-source cases remain visible and explicitly incompatible.
    interface_manifest = _read(root / "interface-manifest.yaml")
    additional = [(e, "direct-exit") for e in manifest.get("cases") or []]
    additional += [(e, "downstream-interface") for e in interface_manifest.get("cases") or []]
    for entry, declared_source in additional:
        if not entry.get("path") or entry.get("Kn") is None:
            continue
        resolved = str(_local(root, entry["path"]))
        if resolved not in paths:
            merged = dict(entry, registration="local-generated-manifest", declared_source_type=declared_source)
            rows.append(_case_row(root, merged))
            paths.add(resolved)
    groups = {}
    for row in rows:
        if row["compatible"]:
            groups.setdefault(row["stratum"], []).append(row["name"])
    missing = [k for k in TARGET_DECADES if not any(r["fit_ready"] and policy.same_knudsen(k, r["Kn"]) for r in rows)]
    prepared = []
    for entry, cfg in specs:
        row = next(r for r in rows if r["path"] == str(_local(root, entry["path"])))
        if not row["fit_ready"]:
            prepared.append({"name": entry["name"], "Kn": entry["Kn"],
                             "configuration": dataclasses.asdict(cfg), "exit_state": row.get("exit_state"),
                             "feasibility": row.get("feasibility"),
                             "known_source_limitations": row.get("known_source_limitations", []),
                             "status": "analysis-input-specification-copy; use the generated study case for execution"})
    unregistered = []
    for path in sorted(root.glob("Cases*/*/case.yaml")):
        if str(path.parent.resolve()) not in paths:
            unregistered.append(str(path.parent.resolve()))
    code = provenance.git_provenance(root.parent.parent)
    return {"schema_version": 1, "study_root": str(root), "revision": code["commit"],
            "code": code, "registration": None if registration_path is None else str(registration_path),
            "active_manifest": str(manifest_path), "cases_dir": study.cases_dir,
            "knudsen_suite": study.meta.get("knudsen_suite", "legacy"),
            "target_decades": list(TARGET_DECADES), "cases": rows,
            "missing_results": missing, "prepared_specs": prepared,
            "compatibility_strata": groups, "unregistered_local_cases": unregistered,
            "external_equivalents": _external_equivalents(root),
            "suitability_criteria": SUITABILITY_CRITERIA,
            "production_policy": "The study profile controls generation and case guards. Scaling registration adds analysis inputs only."}


def prepare_missing_specs(study_root, output_dir, registration_path=None):
    """Write small, reviewable YAML inputs; no meshes, dictionaries, or scripts."""
    document = inventory(study_root, registration_path)
    destination = Path(output_dir).resolve()
    written = []
    for spec in document["prepared_specs"]:
        directory = destination / spec["name"]
        directory.mkdir(parents=True, exist_ok=True)
        cfg = copy.deepcopy(spec["configuration"])
        cfg["meta"].update({"scaling_specification_only": True,
                            "production_execution_authorized": False})
        for filename, value in (("case.yaml", cfg), ("specification.yaml", {k: v for k, v in spec.items() if k != "configuration"})):
            path = directory / filename
            text = yaml.safe_dump(value, sort_keys=False)
            if path.is_file() and path.read_text(encoding="utf-8") != text:
                raise InventoryError(f"Refusing to overwrite changed specification {path}")
            path.write_text(text, encoding="utf-8", newline="\n")
            written.append(path)
    return written


def load_case(row, max_bytes=DEFAULT_MAX_BYTES):
    """Load registered matched arrays, or return None for missing/partial results.

    Portable NPZ: points=(X,rho/D), density=n/ne, volumes in m^3;
    optional temperature=T/Te, occupied=meaningful-estimator mask,
    domain_valid=fluid-domain mask and velocity=(ux,uy,uz) in m/s.  NPZ files
    never contain Python objects and are checked against an uncompressed cap.
    """
    if not row.get("fit_ready") or not row.get("compatible"):
        return None
    from plumetools.cai2012ext.scaling_fit import ScalingCase

    source = row["field_source"]
    metadata = {"inventory": row, "uncertainty": row.get("uncertainty"),
                "stratum": row.get("stratum"), "synthetic": False}
    metadata.update((row.get("uncertainty") or {}).get("declared_block_metadata") or {})
    if source["kind"] == "npz":
        path = Path(source["path"])
        _npz_description(path, max_bytes)
        with np.load(path, allow_pickle=False) as values:
            points = np.asarray(values["points"], dtype=float)
            density = np.asarray(values["density"], dtype=float)
            volumes = np.asarray(values["volumes"], dtype=float)
            temperature = np.asarray(values["temperature"], dtype=float) if "temperature" in values else None
            occupied = np.asarray(values["occupied"], dtype=bool) if "occupied" in values else None
            for key in ("domain_valid", "velocity", "velocity_sqrt_beta0", "parcel_count",
                        "cylindrical_velocity_sqrt_beta0", "cartesian_centres_m",
                        "density_blocks", "block_times", "block_durations"):
                if key in values:
                    metadata[key] = np.asarray(values[key])
            for key in BLOCK_METADATA:
                if key in values:
                    if values[key].size != 1:
                        raise InventoryError(f"{path}: {key} must be scalar metadata")
                    metadata[key] = values[key].item()
        if points.ndim != 2 or points.shape[1] != 2 or density.shape != (len(points),) or volumes.shape != density.shape:
            raise InventoryError(f"{path}: mismatched points/density/volumes arrays")
        if occupied is None and "parcel_count" in metadata:
            occupied = metadata["parcel_count"] > 0
        if "velocity" in metadata and "velocity_sqrt_beta0" not in metadata:
            metadata["velocity_sqrt_beta0"] = metadata["velocity"] * np.sqrt(row["exit_state"]["beta0_s2_per_m2"])
    elif source["kind"] == "openfoam":
        directory = Path(source["path"])
        time_dir = directory / source["time"]
        geom_dir = directory / source["mesh_time"]
        # The legacy read_case deliberately permits equal weights if V is absent;
        # this experiment requires true volumes and checks them before calling it.
        required = [geom_dir / "C", geom_dir / "V", time_dir / "rhoNMean"]
        if not all(p.is_file() for p in required):
            raise InventoryError("Whole-field geometry or density disappeared since inventory")
        size = sum(p.stat().st_size for p in time_dir.iterdir()
                   if p.is_file() and p.name in set(post.REQUIRED_FIELDS) | {"dsmcRhoNMean"})
        size += sum(p.stat().st_size for p in required[:2])
        if size * 4 > max_bytes:
            raise InventoryError(f"OpenFOAM input size exceeds conservative array cap {max_bytes}; use HPC configuration")
        centres = np.asarray(post.read_internal_field(geom_dir / "C"), dtype=float)
        volumes = np.asarray(post.read_cell_values(geom_dir / "V", len(centres)), dtype=float)
        state = row["exit_state"]
        density = np.asarray(post.read_cell_values(time_dir / "rhoNMean", len(centres)), dtype=float)
        density = density / state["n0_per_m3"]
        points = np.column_stack((centres[:, 0], np.hypot(centres[:, 1], centres[:, 2]))) / state["D_m"]
        temperature = None
        occupied = density > 0
        parcel_path = time_dir / "dsmcRhoNMean"
        if parcel_path.is_file():
            count = np.asarray(post.read_cell_values(parcel_path, len(centres)), dtype=float)
            metadata["parcel_count"] = count
            occupied = np.isfinite(count) & (count > 0)
        if all((time_dir / name).is_file() for name in post.REQUIRED_FIELDS):
            sampled = post.read_case(directory, row["gas"]["mass_kg"], time=source["time"], mesh_time=source["mesh_time"])
            temperature = sampled.temperature / state["T0_K"]
            metadata["velocity"] = sampled.velocity
            metadata["velocity_sqrt_beta0"] = sampled.velocity * np.sqrt(state["beta0_s2_per_m2"])
        metadata["cartesian_centres_m"] = centres
        metadata["time"] = source["time"]
    else:
        raise InventoryError(f"Unsupported field source {source['kind']!r}")
    return ScalingCase(name=row["name"], kn=row["Kn"], points=points,
        density=density, volumes=volumes, speed_ratio=row["speed_ratio"],
        temperature=temperature, occupied=occupied, suitable=bool(row["suitable"]), metadata=metadata)
