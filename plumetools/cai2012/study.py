r"""``study.yaml`` parsing, the Knudsen matrix, and the generation manifest.

A study is a base case, a list of Knudsen cases, and a list of numerical-particle
levels. Each generated case is a complete, self-contained OpenFOAM case with its
own ``case.yaml``, and that ``case.yaml`` carries the exact Knudsen number the
solver will run at and the exact parcel population it will run with.

The two axes
------------
.. code-block:: text

    kn_cases          the PHYSICAL matrix -- Kn = 100, 0.1, 0.01   [PAPER]
    particle_levels   the NUMERICAL axis  -- 1x, 2x, 5x, 10x parcels

The product is generated: ``Cases/Kn100_np1x``, ``Cases/Kn100_np2x``, and so on.
The physical definition of a case is untouched by its level -- same mesh, same
density, same time step, same transient -- because a numerical-particle
convergence study is only meaningful if the thing being resolved holds still.
The single knob the level turns is ``nEquivalentParticles``:

.. code-block:: text

    nEquivalentParticles = baseline / numerical_particle_multiplier

so ``np1x`` runs at exactly the weight the study has always run at and is the
control. That baseline is never typed in here: it is whatever
:func:`plumetools.cai2012.inflow.derive_run_settings` derives for the case, so
the sweep cannot drift away from the study it extends.

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

from plumetools.cai2012.dictionaries import AVERAGED_FIELDS

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
#:
#: The three multipliers are the same argument on the other axis. A per-case
#: override of ``numerical_particle_multiplier`` would let a directory named
#: ``_np2x`` run at another level; the run-time and output-frequency multipliers
#: are study-wide by construction, because a sweep whose members were averaged
#: over different windows and written at different rates is not comparable
#: member to member. All three have exactly one source: ``study.yaml``.
FORBIDDEN_OVERRIDES = (
    "exit.knudsen",
    "model",
    "dsmc.numerical_particle_multiplier",
    "dsmc.run_time_multiplier",
    "dsmc.output_frequency_multiplier",
)


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
class ParticleLevel:
    """One numerical-particle level: a parcel population, as a multiplier.

    Attributes:
        name: the suffix appended to a physical case name, e.g. ``np2x``.
        multiplier: how many times the baseline parcel population. The particle
            weight is *divided* by it -- one parcel stands for
            ``nEquivalentParticles`` molecules, so the count goes as the
            reciprocal of the weight.
        enabled: whether :meth:`StudyConfig.enabled_particle_levels` includes it.
            A disabled level stays in ``study.yaml`` as the record of the
            intended matrix.
        note: free text, carried into the manifest.
    """

    name: str
    multiplier: float
    enabled: bool = True
    note: str = ""

    @property
    def is_baseline(self) -> bool:
        """Whether this is the 1x level -- the control the others are read against."""
        return self.multiplier == 1.0


@dataclass(frozen=True)
class StudyCase:
    """One generated case: a physical Knudsen case at one particle level.

    The pair is the identity of a directory, and both halves are recoverable
    from it -- from the name, from ``case.yaml``'s ``meta``, and from the
    manifest row. Everything a generator needs of a :class:`KnudsenCase` is
    forwarded, so the two are interchangeable to :func:`clone_cases` and
    :func:`StudyConfig.case_path`.
    """

    knudsen_case: KnudsenCase
    particle_level: ParticleLevel

    @property
    def name(self) -> str:
        """``Kn100_np2x`` -- the physical case, then the level."""
        return f"{self.knudsen_case.name}_{self.particle_level.name}"

    @property
    def knudsen(self) -> float:
        return self.knudsen_case.knudsen

    @property
    def enabled(self) -> bool:
        return self.knudsen_case.enabled and self.particle_level.enabled

    @property
    def overrides(self) -> dict:
        return self.knudsen_case.overrides

    @property
    def note(self) -> str:
        return self.knudsen_case.note

    @property
    def multiplier(self) -> float:
        return self.particle_level.multiplier


@dataclass(frozen=True)
class StudyConfig:
    """A parsed ``study.yaml``."""

    base_case: str = "baseCase"
    cases_dir: str = "Cases"
    kn_cases: tuple = ()
    particle_levels: tuple = ()
    run_time_multiplier: float = 1.0
    output_frequency_multiplier: float = 1.0
    manifest_name: str = "manifest.yaml"
    meta: dict = None

    def __post_init__(self) -> None:
        if self.meta is None:
            object.__setattr__(self, "meta", {})
        if not self.particle_levels:
            # A study with no particle axis is the study as it was before the
            # axis existed: one level, 1x, which divides the weight by one.
            object.__setattr__(self, "particle_levels",
                               (ParticleLevel(name="np1x", multiplier=1.0),))

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

    def enabled_particle_levels(self) -> list:
        """The enabled levels, cheapest first.

        Ascending multiplier, so the 1x control -- the level that reproduces the
        study as it was -- is generated and run first. If the sweep has to be
        cut short, what survives is the part that is comparable with the
        published results.
        """
        return sorted((p for p in self.particle_levels if p.enabled),
                      key=lambda p: (p.multiplier, p.name))

    def all_particle_levels(self) -> list:
        """Every level, enabled or not, in the same deterministic order."""
        return sorted(self.particle_levels, key=lambda p: (p.multiplier, p.name))

    def expand(self, cases=None, levels=None) -> list:
        """The product of the physical cases and the particle levels.

        Args:
            cases: physical cases to expand; default the enabled ones.
            levels: levels to expand over; default the enabled ones.

        Returns:
            :class:`StudyCase` objects, physical case major and level minor --
            so a case's four variants are adjacent, which is the order they are
            compared in.
        """
        cases = self.enabled_cases() if cases is None else list(cases)
        levels = self.enabled_particle_levels() if levels is None else list(levels)
        return [StudyCase(knudsen_case=c, particle_level=p)
                for c in cases for p in levels]

    def expanded_cases(self) -> list:
        """Every enabled generated case: ``len(kn_cases) * len(particle_levels)``."""
        return self.expand()

    def case_path(self, case) -> Path:
        """Path of one case relative to ``study.yaml``: ``Cases/Kn100_np1x``."""
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

    known = {"base_case", "cases_dir", "kn_cases", "particle_levels",
             "run_time_multiplier", "output_frequency_multiplier",
             "manifest_name", "meta"}
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

    levels = _load_particle_levels(raw.get("particle_levels"), path)

    return StudyConfig(
        base_case=str(raw.get("base_case", "baseCase")),
        cases_dir=str(raw.get("cases_dir", "Cases")),
        kn_cases=tuple(cases),
        particle_levels=tuple(levels),
        run_time_multiplier=_positive_multiplier(
            raw, "run_time_multiplier", path),
        output_frequency_multiplier=_positive_multiplier(
            raw, "output_frequency_multiplier", path),
        manifest_name=str(raw.get("manifest_name", "manifest.yaml")),
        meta=raw.get("meta") or {},
    )


def _positive_multiplier(raw: dict, key: str, path) -> float:
    """Read a top-level multiplier, defaulting to 1.0 -- i.e. no change."""
    if key not in raw or raw[key] is None:
        return 1.0
    try:
        value = float(raw[key])
    except (TypeError, ValueError):
        raise StudyError(f"{path}: {key} is {raw[key]!r}, which is not a "
                         f"number") from None
    if not value > 0.0:
        raise StudyError(
            f"{path}: {key} is {value}; it scales a duration or a write "
            f"frequency, so it must be positive. 1.0 leaves the baseline "
            f"schedule alone.")
    return value


def _load_particle_levels(entries, path) -> list:
    """Parse ``particle_levels``, the numerical-particle axis.

    An absent or empty list is not an error: it means the study has no such
    axis, and :class:`StudyConfig` supplies the single 1x level that reproduces
    the study as it was before the axis existed.

    Raises:
        StudyError: for a malformed entry, an unknown key, a non-positive
            multiplier, or a duplicate name or multiplier.
    """
    if entries is None:
        return []
    if not isinstance(entries, list):
        raise StudyError(
            f"{path}: particle_levels must be a list of "
            f"{{name, multiplier, enabled}} mappings, got "
            f"{type(entries).__name__}")

    known = {"name", "multiplier", "enabled", "note"}
    levels = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise StudyError(
                f"{path}: particle_levels[{i}] is a {type(entry).__name__}, "
                f"expected a mapping")
        extra = sorted(set(entry) - known)
        if extra:
            raise StudyError(
                f"{path}: particle_levels[{i}] has unknown key(s) {extra}; "
                f"valid: {sorted(known)}")
        for required in ("name", "multiplier"):
            if required not in entry:
                raise StudyError(
                    f"{path}: particle_levels[{i}] is missing {required!r}")
        multiplier = float(entry["multiplier"])
        if not multiplier > 0.0:
            raise StudyError(
                f"{path}: particle_levels[{i}] has multiplier {multiplier}; it "
                f"DIVIDES nEquivalentParticles, so zero or negative would give "
                f"an infinite or negative particle weight")
        levels.append(ParticleLevel(
            name=str(entry["name"]),
            multiplier=multiplier,
            enabled=bool(entry.get("enabled", True)),
            note=str(entry.get("note", "")),
        ))

    names = [p.name for p in levels]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise StudyError(
            f"{path}: duplicate particle level name(s) {duplicates}; each "
            f"becomes a directory suffix, so the second would overwrite the "
            f"first")

    multipliers = [p.multiplier for p in levels]
    duplicate_multipliers = sorted(
        {m for m in multipliers if multipliers.count(m) > 1})
    if duplicate_multipliers:
        raise StudyError(
            f"{path}: duplicate particle multiplier(s) {duplicate_multipliers}; "
            f"two levels at the same parcel population would differ in name "
            f"only")

    if levels and not any(p.is_baseline for p in levels):
        raise StudyError(
            f"{path}: particle_levels has no 1x level. The 1x case IS the "
            f"existing study -- it runs at exactly the nEquivalentParticles the "
            f"physical case has always used -- and without it the sweep has no "
            f"control to be read against. Multipliers: {multipliers}")

    return levels


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


def filter_by_particles(levels: list, wanted) -> list:
    """Keep only the levels whose multiplier or name appears in ``wanted``.

    Accepts either form -- ``2`` or ``np2x`` -- because both are what the
    generated tree calls the level.

    Raises:
        StudyError: if a requested level matches none, so a typo on the command
            line stops the run instead of quietly generating nothing.
    """
    if not wanted:
        return list(levels)

    requested = [str(w) for w in wanted]
    kept, missing = [], []
    for token in requested:
        match = None
        for level in levels:
            if level.name == token:
                match = level
                break
            try:
                if level.multiplier == float(token):
                    match = level
                    break
            except ValueError:
                continue
        if match is None:
            missing.append(token)
        elif match not in kept:
            kept.append(match)

    if missing:
        raise StudyError(
            f"no particle level {missing}; the study defines "
            f"{[(p.name, p.multiplier) for p in levels]}")
    return [p for p in levels if p in kept]


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


def apply_case_parameters(case_dir: Path, case, study: StudyConfig) -> dict:
    """Rewrite a case-local ``case.yaml`` with this case's parameters.

    The Knudsen number, the numerical-particle multiplier, and the study-wide
    run-time and output-frequency multipliers. Nothing else: the density, the
    mesh, the time step and the transient are all derived downstream from what
    is already in the template, which is what makes the four variants of a case
    the same physical problem.

    Args:
        case_dir: the generated case directory.
        case: a :class:`StudyCase`, or a bare :class:`KnudsenCase` for a study
            with no particle axis.
        study: the study, for the multipliers and the metadata.

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

    # The numerical axis. Written here rather than left to the template so that
    # a generated case.yaml states its own level: the directory name and the
    # file agree, and neither has to be trusted over the other.
    level = getattr(case, "particle_level", None)
    multiplier = float(level.multiplier) if level is not None else 1.0
    if not isinstance(data.get("dsmc"), dict):
        raise StudyError(
            f"{path}: no 'dsmc' section. This does not look like the cai2012 "
            f"template; the particle multiplier has nowhere to go.")
    data["dsmc"]["numerical_particle_multiplier"] = multiplier
    data["dsmc"]["run_time_multiplier"] = float(study.run_time_multiplier)
    data["dsmc"]["output_frequency_multiplier"] = float(
        study.output_frequency_multiplier)
    applied.update({
        "numerical_particle_multiplier": multiplier,
        "run_time_multiplier": float(study.run_time_multiplier),
        "output_frequency_multiplier": float(study.output_frequency_multiplier),
    })
    if level is not None:
        applied["particle_level"] = level.name

    if case.overrides:
        _merge(data, copy.deepcopy(case.overrides), f"{path}")
        applied["overrides"] = copy.deepcopy(case.overrides)

    data["meta"] = dict(data.get("meta") or {})
    data["meta"].update({
        "case_name": case.name,
        "knudsen": float(case.knudsen),
        "generated_by": "cases/cai2012/generate_cases.py",
    })
    if level is not None:
        # Both halves of the identity, recoverable without parsing the name.
        data["meta"].update({
            "cai_case": case.knudsen_case.name,
            "particle_level": level.name,
            "numerical_particle_multiplier": multiplier,
        })
        if level.note:
            data["meta"]["particle_level_note"] = level.note
    if case.note:
        data["meta"]["note"] = case.note

    header = (
        "# GENERATED by cases/cai2012/generate_cases.py -- do not edit.\n"
        "#\n"
        f"# {case.name}: Kn = {case.knudsen:g}, {multiplier:g}x numerical particles.\n"
        "#\n"
        f"# nEquivalentParticles = baseline / {multiplier:g}. The BASELINE is derived\n"
        "# from this case's own density and cell volume and is not stored, so the\n"
        "# 1x case runs at exactly the weight this study has always used.\n"
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
        "n_physical_cases": len(study.enabled_cases()),
        "n_particle_levels": len(study.enabled_particle_levels()),
        "run_time_multiplier": float(study.run_time_multiplier),
        "output_frequency_multiplier": float(study.output_frequency_multiplier),
        "particle_levels": [
            {"name": p.name, "multiplier": float(p.multiplier),
             "enabled": bool(p.enabled), "note": p.note}
            for p in study.all_particle_levels()
        ],
        "averaged_fields": {
            name: {"mean": True, "prime2Mean": bool(prime2mean)}
            for name, prime2mean in AVERAGED_FIELDS.items()
        },
        "meta": dict(study.meta),
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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by cases/cai2012/generate_cases.py -- do not edit.\n"
        "# The auditable map from case name to physical inputs.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False,
                         default_style=None),
        encoding="utf-8", newline="\n")
    return path


def mesh_identifier(plan) -> str:
    """A short, stable name for *this* mesh.

    Two cases share a mesh when this string matches. It is what makes the "same
    mesh across the four particle levels" claim checkable from the manifest
    alone, without re-deriving anything or diffing ``blockMeshDict``.
    """
    return (f"{int(plan.n_cells)}cells"
            f"-core{plan.core_cell_size_m:.6e}m"
            f"-min{plan.min_cell_size_m:.6e}m")


def manifest_entry(case, study: StudyConfig, cfg, geom,
                   exit_state, plan, run, nozzle=None, particles=None) -> dict:
    """Build one manifest row from a fully derived case.

    Every quantity §7 of the case specification asks for -- ``Kn``, ``D``,
    ``lambda0``, ``T0``, ``U0``, ``n0``, the particle weight, the minimum cell
    size and ``deltaT`` -- plus what it took to get them, plus everything needed
    to say which physical case and which particle level a directory holds.
    """
    level = getattr(case, "particle_level", None)
    physical = getattr(case, "knudsen_case", case)
    entry = {
        "name": case.name,
        "path": str(study.case_path(case)).replace("\\", "/"),
        "cai_case": physical.name,
        "particle_level": level.name if level is not None else "np1x",
        "numerical_particle_multiplier": (
            float(level.multiplier) if level is not None else 1.0),
        "template": study.base_case,
        "mesh_id": mesh_identifier(plan),
        "Kn": float(case.knudsen),
        "D_m": float(geom.diameter_m),
        "lambda0_m": float(exit_state.mean_free_path_m),
        "T0_K": float(exit_state.T0_K),
        "S0": float(exit_state.speed_ratio),
        "U0_m_per_s": float(exit_state.velocity_m_per_s),
        "n0_per_m3": float(exit_state.number_density_per_m3),
        "n_equivalent_particles": float(run.n_equivalent_particles),
        "baseline_n_equivalent_particles": float(
            run.baseline_n_equivalent_particles),
        "exit_particles_per_cell": float(run.exit_particles_per_cell),
        "min_cell_size_m": float(plan.min_cell_size_m),
        "max_cell_size_m": float(plan.max_cell_size_m),
        "core_cell_size_m": float(plan.core_cell_size_m),
        "cell_over_mean_free_path": float(plan.cell_over_mfp),
        "deltaT_s": float(run.delta_t_s),
        "deltaT_source": run.delta_t_source,
        "start_time_s": 0.0,
        "end_time_s": float(run.end_time_s),
        "baseline_end_time_s": float(run.baseline_end_time_s),
        "run_time_multiplier": float(run.run_time_multiplier),
        "n_steps": int(run.n_steps),
        "write_control": "timeStep",
        "write_interval_steps": int(run.write_interval_steps),
        "write_interval_s": float(run.write_interval_s),
        "baseline_write_interval_steps": int(run.baseline_write_interval_steps),
        "output_frequency_multiplier": float(run.output_frequency_multiplier),
        "achieved_output_frequency_multiplier": float(
            run.achieved_output_frequency_multiplier),
        "n_writes": int(run.n_writes),
        "average_start_s": float(run.average_start_s),
        "sampling_time_s": float(run.sampling_time_s),
        "n_sampled_writes": int(run.n_sampled_writes),
        "averaging": {
            "function_object": "fieldAverage1",
            "type": "fieldAverage",
            "write_control": "writeTime",
            "time_start_s": float(run.average_start_s),
            "transient_basis": run.transient_basis,
            "fields": {
                name: {"mean": True, "prime2Mean": bool(prime2mean)}
                for name, prime2mean in AVERAGED_FIELDS.items()
            },
        },
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
    if level is not None and level.note:
        entry["particle_level_note"] = level.note
    return entry
