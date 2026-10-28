r"""Generate the extension's cases idempotently, and write its manifest.

The parent generator, ``cases/cai2012/generate_cases.py``, deletes and re-copies
every selected case that has no time directory yet, and refuses the whole run if
one has results. For a study whose dense cases take a long time to mesh that is
the wrong default: re-running generation should be free.

So every case is first **staged** -- built in a scratch directory exactly as it
would be generated, with the parent family's own
:func:`plumetools.cai2012.study.apply_case_parameters` -- and then compared with
what is on disk:

==============  ===========================================================
``new``         nothing there yet: the staged case is moved into place
``identical``   left alone, byte for byte -- mesh, dictionaries and results
                included. This is the normal outcome of a re-run.
``stale``       the physics is identical but a runner script or the
                ``meta`` block differs from the template. Left alone and
                reported; ``--refresh`` copies the new files in.
``conflict``    an input in ``case.yaml`` differs. **Nothing is touched**;
                ``--regenerate`` replaces the case only if it holds no solver
                output, and ``--overwrite`` replaces it even if it does --
                destroying those results, which is what the flag says.
==============  ===========================================================

A conflict is a configuration conflict, and the response to one is a report,
not a silent rewrite under somebody's results.
"""

from __future__ import annotations

import copy
import datetime
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from plumetools.cai2012 import PAPER, checks, inflow, mesh
from plumetools.cai2012 import study as study_module
from plumetools.cai2012.config import load_case_config
from plumetools.cai2012.dictionaries import AVERAGED_FIELDS
from plumetools.cai2012.geometry import from_config
from plumetools.cai2012ext import (
    EXCLUDED_KNUDSEN,
    EXTENSION_KNUDSEN,
    FAMILY_DIRECTORY,
    GENERATOR,
    PARTICLE_MULTIPLIERS,
    REFERENCE_KNUDSEN,
    feasibility,
    policy,
    provenance,
)
from plumetools.foamio.times import time_directories

#: Files whose content is a physical input, compared as parsed YAML.
#: ``interface.yaml`` exists only in interface-mode cases; a direct-exit case
#: has none, so for those this is the comparison it always was.
INPUT_FILES = ("case.yaml", "interface.yaml")

#: Keys of those files that are provenance, not inputs: a difference is
#: ``stale``, never a conflict. The overlap evidence a production interface case
#: points at is filled in after validation, on a case that may already be meshed.
AUXILIARY_KEYS = ("meta.", "interface.validation_evidence", "interface.label")

STATUS_NEW = "new"
STATUS_IDENTICAL = "identical"
STATUS_STALE = "stale"
STATUS_CONFLICT = "conflict"


# --------------------------------------------------------------------------- #
# solver output
# --------------------------------------------------------------------------- #

def has_solver_output(case_dir: Path) -> bool:
    """Whether a case holds anything a solver run produced.

    Wider than :func:`plumetools.cai2012.study.case_is_complete`, which looks
    only at the reconstructed case: a parallel run that was interrupted before
    ``reconstructPar`` has its results under ``processor*/`` only, and those are
    results too.
    """
    case_dir = Path(case_dir)
    if time_directories(case_dir, missing_ok=True):
        return True
    for processor in case_dir.glob("processor[0-9]*"):
        if time_directories(processor, missing_ok=True):
            return True
    return (case_dir / "log.dsmcFoam").is_file()


# --------------------------------------------------------------------------- #
# staging
# --------------------------------------------------------------------------- #

def write_case_yaml(case_dir: Path, case, study) -> dict:
    """Write a case's ``case.yaml`` -- the parent's way, then this family's meta.

    :func:`plumetools.cai2012.study.apply_case_parameters` sets the Knudsen
    number and the particle level structurally. It also names
    ``cases/cai2012/generate_cases.py`` as the generator, which is not true
    here, so the provenance is rewritten -- again structurally -- and the file
    re-emitted with this family's header. No physical key is touched after the
    parent function has written it.
    """
    applied = study_module.apply_case_parameters(case_dir, case, study)

    path = Path(case_dir) / "case.yaml"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    baseline = float(data["resolution"]["target_particles_per_cell"])
    multiplier = float(data["dsmc"]["numerical_particle_multiplier"])
    data["meta"].update({
        "generated_by": GENERATOR,
        "family": FAMILY_DIRECTORY,
        "study": "Cai 2012 Knudsen extension",
        "baseline_exit_particles_per_cell": baseline,
        "exit_particles_per_cell_target": baseline * multiplier,
        "knudsen_suite": policy.suite_name(study),
        "source_type": "direct-exit",
    })
    cfg = load_case_config(case_dir)
    # The case has not been re-emitted yet; the suite comes from the study.
    from dataclasses import replace
    cfg = replace(cfg, meta={**cfg.meta, "knudsen_suite": policy.suite_name(study)})
    data["meta"]["known_source_limitations"] = policy.source_limitations(cfg)
    header = (
        f"# GENERATED by {GENERATOR} -- do not edit.\n"
        "#\n"
        f"# {case.name}: Kn = {case.knudsen:g}, {multiplier:g}x numerical "
        f"particles ({baseline * multiplier:g} parcels\n"
        f"# in the exit cell). Cai extension profile: {policy.suite_name(study)}.\n"
        "#\n"
        f"# nEquivalentParticles = (n0 V_exit / {baseline:g}) / {multiplier:g}, "
        "derived at mesh time.\n"
        "#\n"
        "# Edit ../../baseCase/case.yaml (or ../../study.yaml) and re-run\n"
        "# ../../generate_cases.py. It never overwrites a case with results.\n"
        "\n"
    )
    path.write_text(
        header + yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    if data["meta"]["known_source_limitations"]:
        (Path(case_dir) / "KNOWN-LIMITATIONS").write_text(
            "PREPARED CASE: source limitation retained; no solver run is implied.\n\n"
            + "\n".join(data["meta"]["known_source_limitations"]) + "\n", encoding="utf-8")
    return applied


def stage(template: Path, staging_root: Path, case, study) -> Path:
    """Build one case in ``staging_root`` exactly as it would be generated."""
    destination = Path(staging_root) / case.name
    shutil.copytree(template, destination)
    study_module.prune_generated_dictionaries(destination)
    write_case_yaml(destination, case, study)
    return destination


def _files(root: Path) -> list[Path]:
    return sorted(p.relative_to(root) for p in Path(root).rglob("*")
                  if p.is_file())


def _flatten(data, prefix: str = "") -> dict:
    flat = {}
    if isinstance(data, dict):
        for key, value in data.items():
            flat.update(_flatten(value, f"{prefix}{key}."))
    else:
        flat[prefix[:-1]] = data
    return flat


def compare(staged: Path, existing: Path) -> tuple[list[str], list[str]]:
    """What differs between a staged case and the one on disk.

    Returns:
        ``(inputs, auxiliary)``. ``inputs`` are ``case.yaml`` keys outside
        ``meta`` whose values differ -- a physics conflict. ``auxiliary`` are
        ``meta`` keys and template files (the runner scripts) that differ or are
        missing. Files the case has and the template does not -- the mesh, the
        dictionaries, the results -- are the case's own and are not compared.
    """
    inputs, auxiliary = [], []
    for relative in _files(staged):
        want = Path(staged) / relative
        have = Path(existing) / relative
        name = str(relative).replace("\\", "/")
        if not have.is_file():
            (inputs if relative.name in INPUT_FILES else auxiliary).append(
                f"{name}: missing")
            continue
        if relative.name in INPUT_FILES:
            mine = _flatten(yaml.safe_load(want.read_text(encoding="utf-8")) or {})
            theirs = _flatten(yaml.safe_load(have.read_text(encoding="utf-8")) or {})
            for key in sorted(set(mine) | set(theirs)):
                if mine.get(key) == theirs.get(key):
                    continue
                entry = f"{key}: {theirs.get(key)!r} -> {mine.get(key)!r}"
                if relative.name != "case.yaml":
                    entry = f"{relative.name}: {entry}"
                (auxiliary if key.startswith(AUXILIARY_KEYS) else inputs).append(entry)
        elif want.read_bytes() != have.read_bytes():
            auxiliary.append(f"{name}: differs from the template")
    return inputs, auxiliary


@dataclass
class CaseState:
    """One case of the matrix, staged and compared with what is on disk."""

    case: object
    destination: Path
    staged: Path
    status: str
    has_results: bool = False
    inputs: list = field(default_factory=list)
    auxiliary: list = field(default_factory=list)
    action: str = ""

    @property
    def blocking(self) -> bool:
        return self.action == "refuse"


def classify(root: Path, study, staging_root: Path, cases: list) -> list[CaseState]:
    """Stage every case and decide what it is relative to the disk."""
    template = Path(root) / study.base_case
    states = []
    for case in cases:
        destination = Path(root) / study.case_path(case)
        staged = stage(template, staging_root, case, study)
        if not destination.exists():
            states.append(CaseState(case, destination, staged, STATUS_NEW))
            continue
        inputs, auxiliary = compare(staged, destination)
        if inputs:
            status = STATUS_CONFLICT
        elif auxiliary:
            status = STATUS_STALE
        else:
            status = STATUS_IDENTICAL
        states.append(CaseState(case, destination, staged, status,
                                has_results=has_solver_output(destination),
                                inputs=inputs, auxiliary=auxiliary))
    return states


def decide(states: list[CaseState], *, refresh: bool = False,
           regenerate: bool = False, overwrite: bool = False) -> list[CaseState]:
    """Set each state's ``action`` from its status and the flags given."""
    for state in states:
        if state.status == STATUS_NEW:
            state.action = "create"
        elif state.status == STATUS_IDENTICAL:
            state.action = "keep"
        elif state.status == STATUS_STALE:
            state.action = "refresh" if refresh else "keep-stale"
        elif state.has_results:
            state.action = "replace" if overwrite else "refuse"
        else:
            state.action = "replace" if (regenerate or overwrite) else "refuse"
    return states


def apply(states: list[CaseState]) -> None:
    """Carry out the decided actions. Call only when nothing is refused."""
    if any(state.blocking for state in states):
        raise policy.ExtensionError(
            "refusing to apply a generation plan with unresolved conflicts")
    for state in states:
        if state.action == "create":
            state.destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(state.staged), str(state.destination))
        elif state.action == "replace":
            shutil.rmtree(state.destination)
            shutil.move(str(state.staged), str(state.destination))
        elif state.action == "refresh":
            for relative in _files(state.staged):
                target = state.destination / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(state.staged / relative, target)


def stray_directories(root: Path, study, cases: list) -> list[dict]:
    """Directories in the cases tree that are not one of the study's cases.

    Reported, never deleted. One that holds an excluded Kn is called out by
    name, because it is the thing this extension must never run -- and since
    every runner reads the manifest rather than the directory, it never will.
    """
    tree = Path(root) / study.cases_dir
    if not tree.is_dir():
        return []
    known = {case.name for case in cases}
    strays = []
    for entry in sorted(tree.iterdir()):
        if not entry.is_dir() or entry.name in known:
            continue
        knudsen = None
        try:
            raw = yaml.safe_load((entry / "case.yaml").read_text(encoding="utf-8"))
            knudsen = float(raw["exit"]["knudsen"])
        except (OSError, TypeError, KeyError, ValueError, yaml.YAMLError):
            pass
        strays.append({
            "path": str(entry.relative_to(root)).replace("\\", "/"),
            "Kn": knudsen,
            "excluded": knudsen is not None and not policy.is_suite_knudsen(knudsen, study),
        })
    return strays


# --------------------------------------------------------------------------- #
# derivation and the manifest
# --------------------------------------------------------------------------- #

def derive(case_dir: Path):
    """Everything the case will run at, from its own ``case.yaml``.

    The same sequence as ``cases/cai2012/generate_cases.py``'s ``derive``.
    """
    cfg = load_case_config(case_dir)
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    particles = mesh.estimate_particles(
        plan, geom, exit_state, run.n_equivalent_particles)
    report = checks.run_checks(cfg, geom, exit_state, plan, run,
                               nozzle=nozzle, particles=particles)
    return cfg, geom, exit_state, plan, run, nozzle, particles, report


def build_entry(root: Path, study, study_path: Path, case, case_dir: Path, *,
                cost_model=None, cost_origin=None) -> tuple[dict, list]:
    """One manifest row: the parent family's, plus the extension's records.

    Returns:
        ``(entry, report_lines)``.
    """
    cfg, geom, exit_state, plan, run, nozzle, particles, report = derive(case_dir)
    entry = study_module.manifest_entry(
        case, study, cfg, geom, exit_state, plan, run,
        nozzle=nozzle, particles=particles)
    assessment = feasibility.assess(
        cfg, geom, exit_state, plan, run, nozzle, particles, report,
        cost_model=cost_model, cost_origin=cost_origin)
    entry.update({
        "knudsen_suite": policy.suite_name(study),
        "source_type": "direct-exit",
        "known_source_limitations": policy.source_limitations(cfg),
        "production_status": assessment["status"],
        "exit": provenance.exit_conditions(cfg, geom, exit_state),
        "gas": provenance.gas_parameters(cfg),
        "weight": provenance.particle_weight(cfg, run),
        "startup": provenance.startup(cfg, run, exit_state),
        "sampling": provenance.averaging(cfg, run),
        "sources": provenance.sources(root, case_dir, Path(root) / study.base_case,
                                      study_path),
        "dsmc_checks": report.as_dict(),
        "feasibility": assessment,
    })
    return entry, feasibility.describe(assessment)


#: Values that must be the same in every case, as manifest-row paths.
HELD_FIXED_ACROSS_STUDY = (
    "exit.D_m", "exit.S0", "exit.T0_K", "exit.characteristic_length",
    "exit.mean_free_path_convention",
    "gas.species_name", "gas.mass_kg", "gas.diameter_m", "gas.omega",
    "gas.t_ref_K", "gas.internal_degrees_of_freedom",
    "gas.binary_collision_model", "gas.collisions_enabled",
    "startup.average_start_s", "startup.reference_knudsen",
    "startup.requested_reference_collision_times",
    "sampling.requested_sampling_domain_transits",
    "sampling.requested_write_interval_s", "sampling.run_time_multiplier",
)

#: Values that must be the same across the particle levels of one Kn. With the
#: baseline only, each Kn is one case and this is trivially met; it is kept so
#: a second level could never slip in on a different mesh unnoticed.
HELD_FIXED_WITHIN_KNUDSEN = (
    "Kn", "n0_per_m3", "mesh_id", "n_cells", "deltaT_s", "average_start_s",
    "end_time_s", "write_interval_steps", "baseline_n_equivalent_particles",
    "weight.baseline_exit_particles_per_cell",
)


def _lookup(entry: dict, dotted: str):
    value = entry
    for part in dotted.split("."):
        if not isinstance(value, dict) or part not in value:
            return None
        value = value[part]
    return value


def held_fixed(entries: list) -> tuple[dict, dict, list]:
    """What is shared across the study and within each Kn, and what is not.

    Returns:
        ``(study_wide, per_knudsen, problems)``. A value that is not shared is
        recorded as ``"NOT SHARED -- [...]"`` and listed in ``problems``.
    """
    problems = []

    def shared(rows, keys, where):
        out = {}
        for key in keys:
            values = []
            for row in rows:
                value = _lookup(row, key)
                if value not in values:
                    values.append(value)
            if len(values) == 1:
                out[key] = values[0]
            else:
                out[key] = f"NOT SHARED -- {values}"
                problems.append(f"{where}: {key} differs: {values}")
        return out

    study_wide = shared(entries, HELD_FIXED_ACROSS_STUDY, "study")
    groups = {}
    for entry in entries:
        groups.setdefault(entry["cai_case"], []).append(entry)
    per_knudsen = {}
    for name, rows in groups.items():
        per_knudsen[name] = shared(rows, HELD_FIXED_WITHIN_KNUDSEN, name)
        per_knudsen[name]["particle_levels"] = [
            row["particle_level"] for row in rows]
        per_knudsen[name]["production_status"] = sorted(
            {row["production_status"] for row in rows})

        # Within a Kn the weights must be the baseline divided by the level.
        for row in rows:
            want = float(row["baseline_n_equivalent_particles"]) / float(
                row["numerical_particle_multiplier"])
            got = float(row["n_equivalent_particles"])
            if abs(got - want) > 1.0e-12 * abs(want):
                problems.append(
                    f"{row['name']}: nEquivalentParticles {got:.9e} is not the "
                    f"baseline / {row['numerical_particle_multiplier']:g} = "
                    f"{want:.9e}")
    return study_wide, per_knudsen, problems


def write_manifest(path: Path, study, entries: list, *, study_path: Path,
                   inheritance: dict | None, code: dict, strays: list) -> Path:
    """Write the extension's manifest.

    Keeps every key the parent family's manifest has -- so
    :func:`plumetools.cai2012.verify.verify_tree` and the parent's runner
    scripts read it unchanged -- and adds the extension's records. Written to a
    temporary file and renamed, so an interrupted write never leaves a
    half-manifest for a runner to schedule from.
    """
    path = Path(path)
    study_wide, per_knudsen, problems = held_fixed(entries)
    counts = {}
    for entry in entries:
        counts[entry["production_status"]] = counts.get(
            entry["production_status"], 0) + 1
    document = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "reference": PAPER["reference"],
        "study": f"Cai 2012 Knudsen extension ({policy.suite_name(study)})",
        "generator": GENERATOR,
        "family": FAMILY_DIRECTORY,
        "study_file": str(study_path.name),
        "base_case": study.base_case,
        "cases_dir": study.cases_dir,
        "n_generated": len(entries),
        "n_physical_cases": len(study.enabled_cases()),
        "n_particle_levels": len(study.enabled_particle_levels()),
        "knudsen_suite": policy.suite_name(study),
        "extension_knudsen": list(policy.suite_knudsen(study)),
        "excluded_knudsen": ([] if policy.suite_name(study) == "full_decades"
                              else list(EXCLUDED_KNUDSEN)),
        "particle_multipliers": list(PARTICLE_MULTIPLIERS),
        "reference_knudsen": REFERENCE_KNUDSEN,
        "run_time_multiplier": float(study.run_time_multiplier),
        "output_frequency_multiplier": float(study.output_frequency_multiplier),
        "production_status_counts": counts,
        "held_fixed": study_wide,
        "per_knudsen": per_knudsen,
        "held_fixed_problems": problems,
        "inheritance": inheritance,
        "code": code,
        "stray_directories": strays,
        "particle_levels": [
            {"name": p.name, "multiplier": float(p.multiplier),
             "enabled": bool(p.enabled), "note": p.note}
            for p in study.all_particle_levels()
        ],
        "averaged_fields": {
            name: {"mean": True, "prime2Mean": bool(prime2mean)}
            for name, prime2mean in AVERAGED_FIELDS.items()
        },
        "meta": copy.deepcopy(dict(study.meta)),
        "cases": entries,
        "disabled_cases": [
            {"name": c.name, "Kn": c.knudsen, "note": c.note}
            for c in study.all_cases() if not c.enabled
        ],
        "disabled_particle_levels": [
            {"name": p.name, "multiplier": float(p.multiplier), "note": p.note}
            for p in study.all_particle_levels() if not p.enabled
        ],
    }
    text = (f"# GENERATED by {GENERATOR} -- do not edit.\n"
            "# The auditable map from case name to inputs, provenance and "
            "production status.\n\n"
            + yaml.safe_dump(document, sort_keys=False, default_flow_style=False))
    path.parent.mkdir(parents=True, exist_ok=True)
    scratch = path.with_name(path.name + ".tmp")
    scratch.write_text(text, encoding="utf-8", newline="\n")
    os.replace(scratch, path)
    return path
