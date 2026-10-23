r"""The interface track: ``interfaces.yaml``, its cases, and their manifest.

A parallel tree beside ``Cases/``, so the four-case direct-exit matrix and its
manifest are exactly what they were::

    cases/cai2012-kn-extension/
        study.yaml, Cases/, manifest.yaml     the direct exit (unchanged)
        interfaces.yaml                       this track's inputs
        Interface/<variant>/                  baseCase + baseCaseInterface/
        interface-manifest.yaml               this track's manifest

Three groups, each a different kind of evidence and labelled as such:

``production``
    the two low-Kn cases of the matrix (1e-3, 1e-4), source-flow injection at
    one threshold. Named apart from ``Cases/`` (``..._sphere``) so a direct
    and an interface run can never overwrite each other.
``validation``
    truncated runs at Kn = 0.01 -- Cai's, so excluded from ``Cases/`` -- one per
    Bird threshold, compared with the existing full-domain benchmark. Allowed
    here and only here; they never become production cases.
``smoke``
    tiny development cases. Provisional angles and radii are allowed only
    here, and nothing they produce is evidence.

Generation is :mod:`.generation`'s: stage, compare (``case.yaml`` and
``interface.yaml`` both inputs), keep identical cases, refuse conflicts,
never touch results.
"""

from __future__ import annotations

import copy
import datetime
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from plumetools.cai2012 import study as study_module
from plumetools.cai2012ext import (
    EXTENSION_KNUDSEN,
    GENERATOR,
    generation,
    interface,
    interface_case,
    policy,
)

#: The production Kn that may use the interface: the matrix's dense pair.
PRODUCTION_KNUDSEN = (1.0e-3, 1.0e-4)

#: Per-variant ``case.yaml`` overrides. Never ``exit``, ``gas`` or ``nozzle``:
#: D, S0, T0 and the gas are the experiment.
SCIENTIFIC_OVERRIDES = {"sections": ("geometry", "mesh"),
                        "keys": ("dsmc.courant_target", "dsmc.delta_t_s",
                                 "dsmc.n_subdomains")}
SMOKE_OVERRIDES = {"sections": ("geometry", "mesh", "post"),
                   "keys": ("dsmc.courant_target", "dsmc.delta_t_s",
                            "dsmc.n_subdomains", "dsmc.transient_basis",
                            "dsmc.transient_domain_transits",
                            "dsmc.sampling_domain_transits",
                            "dsmc.write_interval_s",
                            "resolution.target_particles_per_cell")}

TOP_KEYS = ("overlay", "tree", "manifest_name", "resources", "production",
            "validation", "smoke", "meta")


@dataclass(frozen=True)
class Variant:
    """One interface case. Duck-types the parent's ``StudyCase`` for generation."""

    name: str
    group: str
    knudsen: float
    knudsen_case: object
    particle_level: object
    overrides: dict
    interface: dict
    overlap: dict | None = None
    note: str = ""

    @property
    def multiplier(self) -> float:
        return float(self.particle_level.multiplier)


@dataclass
class InterfaceStudy:
    overlay: str
    tree: str
    manifest_name: str
    resources: dict
    variants: list = field(default_factory=list)
    validation_thresholds: tuple = ()
    meta: dict = field(default_factory=dict)


def _token(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _check_overrides(overrides, allowed, where) -> list[str]:
    problems = []
    if not isinstance(overrides, dict):
        return [f"{where}: overrides must be a mapping"]
    for section, values in overrides.items():
        if section in allowed["sections"]:
            continue
        keys = values if isinstance(values, dict) else {None: values}
        for key in keys:
            dotted = f"{section}.{key}" if key is not None else section
            if dotted not in allowed["keys"]:
                problems.append(
                    f"{where}: may not override {dotted!r}; allowed: sections "
                    f"{list(allowed['sections'])}, keys {list(allowed['keys'])}")
    try:
        study_module._reject_forbidden_overrides(overrides, where)
    except study_module.StudyError as exc:
        problems.append(str(exc))
    return problems


def _interface_block(raw, where, purpose) -> tuple[dict, list[str]]:
    block = dict(raw or {})
    block.setdefault("enabled", True)
    block.setdefault("purpose", purpose)
    problems = []
    if block["purpose"] != purpose:
        problems.append(f"{where}: interface.purpose must be {purpose!r}, "
                        f"got {block['purpose']!r}")
    try:
        interface.parse_interface(block, where)
    except interface.InterfaceConfigError as exc:
        problems.append(str(exc))
    return block, problems


def load(path: Path, study) -> InterfaceStudy:
    """Read and police ``interfaces.yaml`` against the direct-exit ``study``.

    Raises:
        policy.ExtensionError: listing every problem.
    """
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    problems = []
    unknown = sorted(set(raw) - set(TOP_KEYS))
    if unknown:
        problems.append(f"unknown top-level key(s) {unknown}; valid: {list(TOP_KEYS)}")
    level = next(p for p in study.enabled_particle_levels() if p.is_baseline)
    kn_cases = {c.name: c for c in study.enabled_cases()}
    resources = raw.get("resources") or {}
    try:
        interface_case._check_resources(resources, path)
    except interface.InterfaceConfigError as exc:
        problems.append(str(exc))

    variants = []
    for i, entry in enumerate(raw.get("production") or []):
        where = f"production[{i}]"
        name = entry.get("name")
        kn_case = kn_cases.get(entry.get("kn_case"))
        if kn_case is None:
            problems.append(f"{where}: kn_case {entry.get('kn_case')!r} is not an "
                            f"enabled case of study.yaml")
            continue
        if not any(policy.same_knudsen(kn_case.knudsen, k) for k in PRODUCTION_KNUDSEN):
            problems.append(f"{where}: Kn = {kn_case.knudsen:g}; the interface is "
                            f"enabled for production only at {list(PRODUCTION_KNUDSEN)}")
        block, more = _interface_block(entry.get("interface"), where, "production")
        problems += more
        problems += _check_overrides(entry.get("overrides") or {}, SCIENTIFIC_OVERRIDES, where)
        variants.append(Variant(name=str(name), group="production",
                                knudsen=float(kn_case.knudsen), knudsen_case=kn_case,
                                particle_level=level,
                                overrides=copy.deepcopy(entry.get("overrides") or {}),
                                interface=block, note=str(entry.get("note", ""))))

    thresholds = ()
    validation = raw.get("validation")
    if validation:
        where = "validation"
        knudsen = float(validation.get("knudsen", float("nan")))
        if not any(policy.same_knudsen(knudsen, k)
                   for k in interface_case.VALIDATION_KNUDSEN):
            problems.append(f"{where}: knudsen must be one of "
                            f"{list(interface_case.VALIDATION_KNUDSEN)} (the "
                            f"full-domain benchmark); got {knudsen:g}")
        thresholds = tuple(float(t) for t in validation.get("bird_thresholds") or [])
        if not thresholds or len(set(thresholds)) != len(thresholds):
            problems.append(f"{where}: bird_thresholds must be distinct numbers")
        if "bird_threshold" in (validation.get("interface") or {}):
            problems.append(f"{where}: set the thresholds in bird_thresholds, not "
                            f"interface.bird_threshold")
        problems += _check_overrides(validation.get("overrides") or {},
                                     SCIENTIFIC_OVERRIDES, where)
        overlap = validation.get("overlap")
        if overlap is None:
            problems.append(f"{where}: an overlap section is required")
        prefix = str(validation.get("name_prefix", f"val_Kn{_token(knudsen)}"))
        reference_case = study_module.KnudsenCase(
            name=f"Kn{_token(knudsen)}", knudsen=knudsen,
            note="interface validation against the full-domain benchmark")
        for threshold in thresholds:
            block, more = _interface_block(
                {**(validation.get("interface") or {}), "bird_threshold": threshold},
                f"{where} (B = {threshold:g})", "interface-validation")
            problems += more
            variants.append(Variant(
                name=f"{prefix}_B{_token(threshold)}", group="validation",
                knudsen=knudsen, knudsen_case=reference_case, particle_level=level,
                overrides=copy.deepcopy(validation.get("overrides") or {}),
                interface=block, overlap=copy.deepcopy(overlap),
                note=f"overlap validation, B = {threshold:g}"))

    for i, entry in enumerate(raw.get("smoke") or []):
        where = f"smoke[{i}]"
        knudsen = float(entry.get("knudsen", float("nan")))
        if not any(policy.same_knudsen(knudsen, k) for k in EXTENSION_KNUDSEN):
            problems.append(f"{where}: smoke cases run at an extension Kn "
                            f"{list(EXTENSION_KNUDSEN)}")
        block, more = _interface_block(entry.get("interface"), where,
                                       "development-smoke")
        problems += more
        problems += _check_overrides(entry.get("overrides") or {}, SMOKE_OVERRIDES, where)
        variants.append(Variant(
            name=str(entry.get("name")), group="smoke", knudsen=knudsen,
            knudsen_case=study_module.KnudsenCase(name=f"Kn{_token(knudsen)}",
                                                  knudsen=knudsen),
            particle_level=level, overrides=copy.deepcopy(entry.get("overrides") or {}),
            interface=block, note="DEVELOPMENT SMOKE -- not evidence"))

    names = [v.name for v in variants]
    direct = {c.name for c in study.expand()}
    for name in names:
        if names.count(name) > 1:
            problems.append(f"variant name {name!r} is used twice")
        if name in direct or name in ("None", ""):
            problems.append(f"variant name {name!r} is missing or collides with a "
                            f"direct-exit case; interface cases are named apart")
    if problems:
        raise policy.ExtensionError("interfaces.yaml:\n  - " + "\n  - ".join(
            dict.fromkeys(problems)))
    return InterfaceStudy(
        overlay=str(raw.get("overlay", "baseCaseInterface")),
        tree=str(raw.get("tree", interface_case.TREE)),
        manifest_name=str(raw.get("manifest_name", interface_case.MANIFEST)),
        resources=dict(resources), variants=variants,
        validation_thresholds=thresholds, meta=dict(raw.get("meta") or {}))


# --------------------------------------------------------------------------- #
# staging
# --------------------------------------------------------------------------- #

def write_interface_yaml(case_dir: Path, variant: Variant, resources: dict) -> Path:
    document = {"interface": variant.interface}
    if variant.overlap is not None:
        document["overlap"] = variant.overlap
    if resources:
        document["resources"] = resources
    path = Path(case_dir) / interface.SETTINGS_FILE
    path.write_text(
        f"# GENERATED by {GENERATOR.replace('generate_cases', 'generate_interfaces')}"
        " from ../../interfaces.yaml -- do not edit.\n"
        f"# {variant.name}: {variant.group}, Kn = {variant.knudsen:g}, interface "
        f"mode ({interface.MODEL_NAME}).\n"
        "# Edit ../../interfaces.yaml and re-run ../../generate_interfaces.py.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


def stage(root: Path, layout: InterfaceStudy, study, staging_root: Path,
          variant: Variant) -> Path:
    """``baseCase`` + the overlay, the parent's ``case.yaml``, ``interface.yaml``."""
    destination = Path(staging_root) / variant.name
    shutil.copytree(Path(root) / study.base_case, destination)
    overlay = Path(root) / layout.overlay
    for item in sorted(overlay.iterdir()):
        if item.is_file():
            shutil.copy2(item, destination / item.name)
    study_module.prune_generated_dictionaries(destination)
    for stale in ("topoSetDict", "createPatchDict"):
        (destination / "system" / stale).unlink(missing_ok=True)
    generation.write_case_yaml(destination, variant, study)
    path = destination / "case.yaml"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    data["meta"].update({"interface_group": variant.group,
                         "interface_model": interface.MODEL_NAME,
                         "case_name": variant.name})
    text = path.read_text(encoding="utf-8")
    header = text[:text.index("\n\n") + 2]
    path.write_text(header + yaml.safe_dump(data, sort_keys=False,
                                            default_flow_style=False),
                    encoding="utf-8", newline="\n")
    write_interface_yaml(destination, variant, layout.resources)
    return destination


def classify(root: Path, layout: InterfaceStudy, study, staging_root: Path) -> list:
    states = []
    for variant in layout.variants:
        destination = Path(root) / layout.tree / variant.name
        staged = stage(root, layout, study, staging_root, variant)
        if not destination.exists():
            states.append(generation.CaseState(variant, destination, staged,
                                               generation.STATUS_NEW))
            continue
        inputs, auxiliary = generation.compare(staged, destination)
        status = (generation.STATUS_CONFLICT if inputs else
                  generation.STATUS_STALE if auxiliary else generation.STATUS_IDENTICAL)
        states.append(generation.CaseState(
            variant, destination, staged, status,
            has_results=generation.has_solver_output(destination),
            inputs=inputs, auxiliary=auxiliary))
    return states


# --------------------------------------------------------------------------- #
# the manifest
# --------------------------------------------------------------------------- #

def validation_summary(entries: list, layout: InterfaceStudy) -> dict | None:
    rows = [e for e in entries if e["group"] == "validation"]
    if not rows:
        return None
    floor = [e.get("max_core_bird_at_r_min") for e in rows
             if e.get("max_core_bird_at_r_min") is not None]
    axis = [e.get("axis_bird_at_r_min") for e in rows
            if e.get("axis_bird_at_r_min") is not None]
    feasible = [e["name"] for e in rows if e["selection_status"] == interface.STATUS_ADMISSIBLE]
    return {
        "knudsen": rows[0]["Kn"],
        "thresholds": list(layout.validation_thresholds),
        "statuses": {e["name"]: e["selection_status"] for e in rows},
        "admissible": feasible,
        "B_at_smallest_admissible_radius": {
            "axis": min(axis) if axis else None,
            "max_over_core": min(floor) if floor else None,
            "note": "the least threshold that could admit any sphere at this "
                    "Kn, core angle and source-validity radius; independent "
                    "of the DSMC results",
        },
        "overlap_study_possible": bool(feasible),
    }


def write_manifest(path: Path, layout: InterfaceStudy, entries: list, *,
                   code: dict) -> Path:
    counts = {}
    for entry in entries:
        counts[entry["status"]] = counts.get(entry["status"], 0) + 1
    document = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "study": "Cai 2012 Knudsen extension -- INTERFACE TRACK",
        "model": interface.MODEL_NAME,
        "assumptions": list(interface.ASSUMPTIONS),
        "tree": layout.tree,
        "n_cases": len(entries),
        "status_counts": counts,
        "validation": validation_summary(entries, layout),
        "code": code,
        "meta": layout.meta,
        "cases": entries,
    }
    text = ("# GENERATED by generate_interfaces.py -- do not edit.\n"
            + yaml.safe_dump(interface_case._plain(document), sort_keys=False,
                             default_flow_style=False))
    path = Path(path)
    scratch = path.with_name(path.name + ".tmp")
    scratch.write_text(text, encoding="utf-8", newline="\n")
    os.replace(scratch, path)
    return path
