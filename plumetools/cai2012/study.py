r"""``study.yaml`` parsing, the Knudsen matrix, and the generation manifest.

A study is a base case plus a list of Knudsen cases. Each generated case is a
complete, self-contained OpenFOAM case with its own ``case.yaml``, and that
``case.yaml`` carries the exact Knudsen number the solver will run at.

How cases are generated
-----------------------
1. ``baseCase`` is **copied** into ``Cases/<name>`` -- :func:`clone_cases`;
2. each case-local ``case.yaml`` is rewritten **structurally**: loaded as YAML,
   edited as a data structure, re-emitted -- :func:`apply_case_parameters`;
3. every derived value is recorded in ``manifest.yaml``.

Step 2 is never string substitution. ``util/caseFoam/genCases.py`` does it the
other way, through ``'#!stringManipulation'``::

    'system/controlDict': {'#!stringManipulation':
                            {'deltaT          1.0E-05': '%s' % deltaT}}

which silently does nothing if the file is reformatted, if the value already
differs, or if the base case has moved on -- and none of those is visible in the
output.

Why this does not use CaseFoam
------------------------------
:mod:`plumetools.markelov1999.study` requires CaseFoam, and says why: that study
is a *hierarchy* (gap over pressure) and CaseFoam builds hierarchies. This study
is a **flat list**. There is no tree for CaseFoam to construct, the physical
values are applied structurally afterwards in both families, and requiring a
dependency that would do a ``copytree`` and nothing else would mean
``./generate_cases.py`` fails on a machine where the core install works. The
departure is deliberate and is recorded in ``cases/cai2012/README.md``.

Traceability
------------
A case labelled ``Kn0p1`` cannot run another case's density, because no density
is stored anywhere. ``case.yaml`` carries ``exit.knudsen`` and nothing else;
``n0`` comes out of :mod:`plumetools.cai2012.gas` at mesh time. ``study.yaml``
may not override ``exit.knudsen`` through the per-case ``overrides`` block --
that is rejected -- so the Knudsen number has exactly one source.
"""

from __future__ import annotations

import copy
import datetime
import shutil
from dataclasses import dataclass
from pathlib import Path

import yaml

#: Dictionaries a case generates from its own ``case.yaml``.
#:
#: A copy is faithful, so one of these left in ``baseCase`` by someone running
#: ``./Allmesh`` there would be inherited by every generated case. ``./Allmesh``
#: would overwrite it, but ``./Allrun`` only checks that
#: ``constant/dsmcProperties`` *exists* -- so an inherited one would let a case
#: run against the template's physics instead of its own.
GENERATED_DICTIONARIES = (
    "system/blockMeshDict",
    "system/topoSetDict",
    "system/createPatchDict",
    "system/meshQualityDict",
    "system/controlDict",
    "system/dsmcInitialiseDict",
    "system/decomposeParDict",
    "system/fvSchemes",
    "system/fvSolution",
    "constant/dsmcProperties",
)

#: Directories a generated case must not inherit from the template.
#:
#: ``constant/polyMesh`` because each case meshes itself; ``0/`` because it holds
#: fields derived from this case's density; ``results/`` and the logs because a
#: fresh case has no results.
GENERATED_DIRECTORIES = ("constant/polyMesh", "0", "results", "postProcessing")

#: ``case.yaml`` paths a study may **not** override.
#:
#: ``exit.knudsen`` is what the study *is*: allowing an override would let a case
#: called ``Kn0p1`` run at another Knudsen number. ``model`` selects the loader.
FORBIDDEN_OVERRIDES = ("exit.knudsen", "model")


class StudyError(ValueError):
    """``study.yaml`` is missing, malformed, or inconsistent."""


@dataclass(frozen=True)
class KnudsenCase:
    """One Knudsen-number case.

    Attributes:
        name: directory name, e.g. ``Kn0p01``.
        knudsen: the Knudsen number. **[PAPER]**
        enabled: whether :meth:`StudyConfig.enabled_cases` includes it. A
            disabled case stays in ``study.yaml`` as the record of the intended
            matrix -- deleting the entry would lose that.
        overrides: ``case.yaml`` values this case changes, as nested sections.
            Recorded verbatim in the manifest.
        note: free text, carried into the manifest.
    """

    name: str
    knudsen: float
    enabled: bool = True
    overrides: dict = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.overrides is None:
            object.__setattr__(self, "overrides", {})


@dataclass(frozen=True)
class StudyConfig:
    """A parsed ``study.yaml``."""

    base_case: str = "baseCase"
    cases_dir: str = "Cases"
    kn_cases: tuple = ()
    manifest_name: str = "manifest.yaml"
    meta: dict = None

    def __post_init__(self) -> None:
        if self.meta is None:
            object.__setattr__(self, "meta", {})

    def enabled_cases(self) -> list:
        """The enabled cases, most rarefied first.

        Sorted by **decreasing** Knudsen number, not by file order: reordering
        ``study.yaml`` must not change the generated tree. Decreasing Kn is also
        the development order -- the collisionless benchmark first, then the
        progressively more collisional cases, which is the sequence §14 of the
        case specification lays out.
        """
        return sorted((c for c in self.kn_cases if c.enabled),
                      key=lambda c: (-c.knudsen, c.name))

    def all_cases(self) -> list:
        """Every case, enabled or not, in the same deterministic order."""
        return sorted(self.kn_cases, key=lambda c: (-c.knudsen, c.name))

    def case_path(self, case: KnudsenCase) -> Path:
        """Path of one case relative to ``study.yaml``: ``Cases/Kn100``."""
        return Path(self.cases_dir) / case.name


def load_study(path: Path) -> StudyConfig:
    """Read and validate ``study.yaml``.

    Args:
        path: the ``study.yaml`` file, or the directory containing it.

    Returns:
        A validated :class:`StudyConfig`.

    Raises:
        StudyError: for a missing file, an unknown key, a duplicate case name or
            Knudsen number, a non-positive Knudsen number, or a forbidden
            override.
    """
    path = Path(path)
    if path.is_dir():
        path = path / "study.yaml"
    if not path.is_file():
        raise StudyError(f"no study.yaml at {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise StudyError(f"{path}: expected a mapping at the top level")

    known = {"base_case", "cases_dir", "kn_cases", "manifest_name", "meta"}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise StudyError(
            f"{path}: unknown top-level key(s) {unknown}; valid: {sorted(known)}")

    entries = raw.get("kn_cases") or []
    if not isinstance(entries, list) or not entries:
        raise StudyError(
            f"{path}: kn_cases must be a non-empty list of "
            f"{{name, Kn, enabled}} mappings")

    cases = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise StudyError(
                f"{path}: kn_cases[{i}] is a {type(entry).__name__}, "
                f"expected a mapping")
        case_keys = {"name", "Kn", "enabled", "overrides", "note"}
        extra = sorted(set(entry) - case_keys)
        if extra:
            raise StudyError(
                f"{path}: kn_cases[{i}] has unknown key(s) {extra}; "
                f"valid: {sorted(case_keys)}")
        for required in ("name", "Kn"):
            if required not in entry:
                raise StudyError(f"{path}: kn_cases[{i}] is missing {required!r}")
        knudsen = float(entry["Kn"])
        if not knudsen > 0.0:
            raise StudyError(
                f"{path}: kn_cases[{i}] has Kn {knudsen}; a Knudsen number is "
                f"positive, and n0 = 1/(sqrt(2) pi d^2 Kn D) would be negative "
                f"or infinite")

        overrides = entry.get("overrides") or {}
        if not isinstance(overrides, dict):
            raise StudyError(
                f"{path}: kn_cases[{i}].overrides must be a mapping of "
                f"case.yaml sections, got {type(overrides).__name__}")
        _reject_forbidden_overrides(overrides, f"{path}: kn_cases[{i}]")

        cases.append(KnudsenCase(
            name=str(entry["name"]),
            knudsen=knudsen,
            enabled=bool(entry.get("enabled", True)),
            overrides=overrides,
            note=str(entry.get("note", "")),
        ))

    names = [c.name for c in cases]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise StudyError(
            f"{path}: duplicate case name(s) {duplicates}; each generates a "
            f"directory, so the second would overwrite the first")

    knudsens = [c.knudsen for c in cases]
    duplicate_kn = sorted({k for k in knudsens if knudsens.count(k) > 1})
    if duplicate_kn:
        raise StudyError(
            f"{path}: duplicate Knudsen number(s) {duplicate_kn}; two cases at "
            f"the same Kn would differ in name only")

    return StudyConfig(
        base_case=str(raw.get("base_case", "baseCase")),
        cases_dir=str(raw.get("cases_dir", "Cases")),
        kn_cases=tuple(cases),
        manifest_name=str(raw.get("manifest_name", "manifest.yaml")),
        meta=raw.get("meta") or {},
    )


def _reject_forbidden_overrides(overrides: dict, where: str) -> None:
    """Raise if an override would change something the study owns."""
    for dotted in FORBIDDEN_OVERRIDES:
        section, _, key = dotted.partition(".")
        if not key:
            if section in overrides:
                raise StudyError(
                    f"{where}: overrides may not set {section!r}. It selects "
                    f"which schema loads the file.")
            continue
        if isinstance(overrides.get(section), dict) and key in overrides[section]:
            raise StudyError(
                f"{where}: overrides may not set {dotted!r}. The Knudsen number "
                f"is what a case IS -- overriding it would let a case named for "
                f"one Kn run at another. Change the 'Kn' entry instead.")


def filter_by_knudsen(cases: list, wanted) -> list:
    """Keep only the cases whose Knudsen number appears in ``wanted``.

    Raises:
        StudyError: if a requested Kn matches no case, so a typo on the command
            line stops the run instead of quietly generating nothing.
    """
    if not wanted:
        return list(cases)
    targets = {float(k) for k in wanted}
    kept = [c for c in cases if c.knudsen in targets]
    missing = sorted(targets - {c.knudsen for c in cases})
    if missing:
        raise StudyError(
            f"no case at Kn {missing}; the study defines "
            f"{sorted(c.knudsen for c in cases)}")
    return kept


# --------------------------------------------------------------------------- #
# generation
# --------------------------------------------------------------------------- #

def clone_cases(root: Path, study: StudyConfig, cases: list) -> list[Path]:
    """Copy ``baseCase`` into ``Cases/<name>`` for each case.

    Args:
        root: the study directory -- the one containing ``baseCase``.
        study: the study.
        cases: the cases to generate.

    Returns:
        The generated case directories, in the order given.

    Raises:
        StudyError: if the template is missing or has no ``case.yaml``.
    """
    root = Path(root)
    template = root / study.base_case
    if not template.is_dir():
        raise StudyError(f"base case {template} does not exist")
    if not (template / "case.yaml").is_file():
        raise StudyError(f"base case {template} has no case.yaml")

    generated = []
    for case in cases:
        destination = root / study.case_path(case)
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(template, destination)
        prune_generated_dictionaries(destination)
        generated.append(destination)
    return generated


def prune_generated_dictionaries(case_dir: Path) -> list[Path]:
    """Remove anything a case must generate from its own ``case.yaml``.

    Template hygiene, not part of the copy. See :data:`GENERATED_DICTIONARIES`
    for why an inherited dictionary is worse than a missing one.
    """
    removed = []
    for relative in GENERATED_DICTIONARIES:
        path = Path(case_dir) / relative
        if path.is_file():
            path.unlink()
            removed.append(path)
    for relative in GENERATED_DIRECTORIES:
        path = Path(case_dir) / relative
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(path)
    for log in Path(case_dir).glob("log.*"):
        log.unlink()
        removed.append(log)
    return removed


def case_is_complete(case_dir: Path) -> bool:
    """Whether a case looks like it has been run.

    True when the case has a time directory other than ``0``, which only the
    solver creates. Used to refuse to overwrite results: regenerating rewrites
    ``case.yaml``, and doing that under finished output would leave results
    whose inputs no longer describe them.
    """
    case_dir = Path(case_dir)
    if not case_dir.is_dir():
        return False
    for entry in case_dir.iterdir():
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if value > 0.0:
            return True
    return False


def _merge(target: dict, updates: dict, where: str) -> None:
    """Merge ``updates`` into ``target`` one section deep.

    One level, not arbitrary depth: ``case.yaml`` is two levels
    (section -> key) and a deeper merge would let an override create a section
    the schema does not have, which the loader would then reject with a message
    about the generated file rather than about ``study.yaml``.
    """
    for section, values in updates.items():
        if section not in target:
            raise StudyError(
                f"{where}: overrides name section {section!r}, which is not in "
                f"the base case.yaml. Sections: {sorted(target)}")
        if not isinstance(values, dict):
            target[section] = values
            continue
        if not isinstance(target[section], dict):
            raise StudyError(
                f"{where}: overrides give {section!r} a mapping, but the base "
                f"case.yaml has a {type(target[section]).__name__} there")
        for key, value in values.items():
            if key not in target[section]:
                raise StudyError(
                    f"{where}: overrides set {section}.{key}, which is not in "
                    f"the base case.yaml. A key the schema does not have would "
                    f"be rejected later, naming the generated file instead of "
                    f"study.yaml.")
            target[section][key] = value


def apply_case_parameters(case_dir: Path, case: KnudsenCase,
                          study: StudyConfig) -> dict:
    """Rewrite a case-local ``case.yaml`` with this case's Knudsen number.

    Args:
        case_dir: the generated case directory.
        case: the Knudsen case.
        study: the study, for the metadata.

    Returns:
        The parameters written, for the manifest.

    Raises:
        StudyError: if the file is missing or is not this family's template --
            rather than adding the missing sections, because a ``case.yaml``
            without an ``exit`` section is not a Cai case and silently growing
            one would hide that.

    **Structured, not textual.** The file is parsed to a dict, specific keys are
    set, and it is re-emitted. The cost is that comments do not survive the round
    trip, so a generated ``case.yaml`` has none -- the right trade for a
    *generated* file. ``baseCase/case.yaml`` keeps its comments and is the one
    anybody edits.
    """
    path = Path(case_dir) / "case.yaml"
    if not path.is_file():
        raise StudyError(f"no case.yaml at {path}")

    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    for section in ("exit", "meta", "nozzle"):
        if section not in data:
            raise StudyError(
                f"{path}: no {section!r} section. This does not look like the "
                f"cai2012 template; generation would silently add keys the rest "
                f"of the file does not expect.")

    data["exit"]["knudsen"] = float(case.knudsen)
    applied = {"knudsen": float(case.knudsen)}

    if case.overrides:
        _merge(data, copy.deepcopy(case.overrides), f"{path}")
        applied["overrides"] = copy.deepcopy(case.overrides)

    data["meta"] = dict(data.get("meta") or {})
    data["meta"].update({
        "case_name": case.name,
        "knudsen": float(case.knudsen),
        "generated_by": "cases/cai2012/generate_cases.py",
    })
    if case.note:
        data["meta"]["note"] = case.note

    header = (
        "# GENERATED by cases/cai2012/generate_cases.py -- do not edit.\n"
        "#\n"
        f"# {case.name}: Kn = {case.knudsen:g}.\n"
        "#\n"
        "# Everything derived from it -- lambda0, n0, the cell size, the particle\n"
        "# weight, the time step -- is computed at mesh time, not stored here.\n"
        "#\n"
        "# Edit ../../baseCase/case.yaml and regenerate. The comments in that file\n"
        "# explain every value; they do not survive the YAML round trip, and\n"
        "# manifest.yaml records what was applied here.\n"
        "\n"
    )
    path.write_text(
        header + yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return applied


def write_manifest(path: Path, study: StudyConfig, entries: list) -> Path:
    """Write the generation manifest.

    Args:
        path: destination file.
        study: the study.
        entries: one dict per generated case.

    Returns:
        The path written.

    The auditable map from case name to physical inputs: given a results
    directory, it says exactly what was generated, from what, at which Knudsen
    number, with which density, cell size, particle weight and time step. It also
    records the **skipped** cases, so a disabled case is visible as a decision
    rather than as an absence.
    """
    path = Path(path)
    document = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "reference": "Cai and Wang 2012",
        "generator": "cases/cai2012/generate_cases.py",
        "base_case": study.base_case,
        "cases_dir": study.cases_dir,
        "n_generated": len(entries),
        "meta": dict(study.meta),
        "cases": entries,
        "disabled_cases": [
            {"name": c.name, "Kn": c.knudsen, "note": c.note}
            for c in study.all_cases() if not c.enabled
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by cases/cai2012/generate_cases.py -- do not edit.\n"
        "# The auditable map from case name to physical inputs.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False,
                         default_style=None),
        encoding="utf-8", newline="\n")
    return path


def manifest_entry(case: KnudsenCase, study: StudyConfig, cfg, geom,
                   exit_state, plan, run, nozzle=None, particles=None) -> dict:
    """Build one manifest row from a fully derived case.

    Every quantity §7 of the case specification asks for -- ``Kn``, ``D``,
    ``lambda0``, ``T0``, ``U0``, ``n0``, the particle weight, the minimum cell
    size and ``deltaT`` -- plus what it took to get them.
    """
    entry = {
        "name": case.name,
        "path": str(study.case_path(case)).replace("\\", "/"),
        "Kn": float(case.knudsen),
        "D_m": float(geom.diameter_m),
        "lambda0_m": float(exit_state.mean_free_path_m),
        "T0_K": float(exit_state.T0_K),
        "S0": float(exit_state.speed_ratio),
        "U0_m_per_s": float(exit_state.velocity_m_per_s),
        "n0_per_m3": float(exit_state.number_density_per_m3),
        "n_equivalent_particles": float(run.n_equivalent_particles),
        "min_cell_size_m": float(plan.min_cell_size_m),
        "max_cell_size_m": float(plan.max_cell_size_m),
        "core_cell_size_m": float(plan.core_cell_size_m),
        "cell_over_mean_free_path": float(plan.cell_over_mfp),
        "deltaT_s": float(run.delta_t_s),
        "end_time_s": float(run.end_time_s),
        "average_start_s": float(run.average_start_s),
        "n_cells": int(plan.n_cells),
        "coarsened_for_budget": bool(plan.coarsened),
        "mean_free_path_convention": exit_state.convention,
    }
    if nozzle is not None:
        entry["nozzle_faces"] = int(nozzle.n_faces)
        entry["nozzle_area_m2"] = float(nozzle.area_m2)
        entry["nozzle_area_error"] = float(nozzle.area_error)
    if particles is not None:
        entry["estimated_particles"] = float(f"{particles['total_particles']:.6g}")
    if case.overrides:
        entry["overrides"] = copy.deepcopy(case.overrides)
    if case.note:
        entry["note"] = case.note
    return entry
