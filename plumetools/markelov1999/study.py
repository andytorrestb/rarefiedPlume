r"""``study.yaml`` parsing, the case matrix, and the generation manifest.

A study is a base case plus a list of pressure cases. Each generated case is a
complete, self-contained OpenFOAM case with its own ``case.yaml``.

How cases are generated
-----------------------
1. **CaseFoam** clones ``baseCase`` into the hierarchy -- :func:`clone_cases`.
2. A **structured Python step** rewrites each case-local ``case.yaml`` -- loaded
   as YAML, edited as a data structure, re-emitted.
3. **plumetools** generates the OpenFOAM dictionaries from that ``case.yaml``.

CaseFoam (https://github.com/DLR-RY/caseFOAM) is a **dependency**, not an
optional accelerant::

    pip install -e ".[cases]"

There is deliberately no built-in substitute. Reimplementing a case generator the
repository already depends on would mean maintaining two, and the one that is not
exercised is the one that drifts. ``util/caseFoam/`` shows this project has used
CaseFoam since before the refactor.

Why the parameters are applied outside CaseFoam
-----------------------------------------------
Step 2 is not string substitution. ``util/caseFoam/genCases.py`` does it the
other way, through CaseFoam's own ``caseData``::

    'system/controlDict': {'#!stringManipulation':
                            {'deltaT          1.0E-05': '%s' % deltaT}}

which silently does nothing if the file is reformatted, if the value already
differs, or if the base case has moved on -- and none of those is visible in the
output.

``caseData`` has two forms and neither fits ``case.yaml``. The dictionary-aware
form goes through PyFoam's ``ParsedParameterFile``, which reads OpenFOAM
dictionaries, not YAML. The other is ``'#!stringManipulation'``, the substitution
above. So CaseFoam does the cloning and the hierarchy -- what it is good at -- and
the physical values are applied structurally by :func:`apply_case_parameters`.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from plumetools.markelov1999.constants import psi_to_pa

#: Dictionaries a case generates from its own ``case.yaml``.
#:
#: CaseFoam copies the template faithfully, so one of these left in ``baseCase``
#: by someone running ``./Allmesh`` there would be inherited by every generated
#: case. See :func:`prune_generated_dictionaries`.
GENERATED_DICTIONARIES = (
    "system/blockMeshDict",
    "system/snappyHexMeshDict",
    "system/meshQualityDict",
    "system/controlDict",
    "system/dsmcInitialiseDict",
    "system/decomposeParDict",
    "constant/dsmcProperties",
)

#: Directories a generated case must not inherit from the template.
#:
#: ``constant/polyMesh`` because each case meshes itself -- inheriting one would
#: silently give every case the template's geometry. ``0/`` because it holds
#: results, and a case that starts with someone else's results is not a fresh
#: case.
GENERATED_DIRECTORIES = ("constant/polyMesh", "0")


class StudyError(ValueError):
    """``study.yaml`` is missing, malformed, or inconsistent."""


@dataclass(frozen=True)
class PressureCase:
    """One reservoir-pressure case.

    Attributes:
        name: directory name, e.g. ``p005psi``.
        pressure_psi: reservoir pressure [psi].
        enabled: whether :func:`enabled_cases` includes it. A disabled case stays
            in ``study.yaml`` as a record of the intended matrix -- deleting the
            entry would lose that.
        note: free text, carried into the manifest.
    """

    name: str
    pressure_psi: float
    enabled: bool = True
    note: str = ""

    @property
    def pressure_pa(self) -> float:
        """Reservoir pressure [Pa], from the exact psi definition."""
        return psi_to_pa(self.pressure_psi)


@dataclass(frozen=True)
class StudyConfig:
    """A parsed ``study.yaml``.

    Attributes:
        base_case: path to the template case, relative to ``study.yaml``.
        cases_dir: where the hierarchy is written, relative to ``study.yaml``.
        gap_in: the cylinder-to-plate separation this study covers [in]. One
            value: the 12 in cases are explicitly out of scope.
        gap_dir: the directory level naming the gap, e.g. ``gap06in``.
        pressure_cases: every case in the matrix, enabled or not.
        manifest_name: filename of the generation manifest.
    """

    base_case: str = "baseCase"
    cases_dir: str = "Cases"
    gap_in: float = 6.0
    gap_dir: str = "gap06in"
    pressure_cases: tuple = ()
    manifest_name: str = "manifest.yaml"

    def enabled_cases(self) -> list[PressureCase]:
        """The enabled cases, in deterministic order.

        Sorted by pressure, not by the order they appear in the file: reordering
        ``study.yaml`` must not change the generated tree, or two people editing
        the same study produce different-looking output from the same matrix.
        """
        return sorted((c for c in self.pressure_cases if c.enabled),
                      key=lambda c: (c.pressure_psi, c.name))

    def all_cases(self) -> list[PressureCase]:
        """Every case, enabled or not, in the same deterministic order."""
        return sorted(self.pressure_cases, key=lambda c: (c.pressure_psi, c.name))

    def case_path(self, case: PressureCase) -> Path:
        """Path of one case relative to ``study.yaml``: ``Cases/gap06in/p025psi``."""
        return Path(self.cases_dir) / self.gap_dir / case.name


def load_study(path: Path) -> StudyConfig:
    """Read and validate ``study.yaml``.

    Args:
        path: the ``study.yaml`` file, or the directory containing it.

    Returns:
        A validated :class:`StudyConfig`.

    Raises:
        StudyError: for a missing file, an unknown key, a duplicate case name or
            pressure, or a non-positive pressure.
    """
    path = Path(path)
    if path.is_dir():
        path = path / "study.yaml"
    if not path.is_file():
        raise StudyError(f"no study.yaml at {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise StudyError(f"{path}: expected a mapping at the top level")

    known = {"base_case", "cases_dir", "gap_in", "gap_dir", "pressure_cases",
             "manifest_name", "meta"}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise StudyError(
            f"{path}: unknown top-level key(s) {unknown}; valid: {sorted(known)}")

    entries = raw.get("pressure_cases") or []
    if not isinstance(entries, list) or not entries:
        raise StudyError(
            f"{path}: pressure_cases must be a non-empty list of "
            f"{{name, pressure_psi, enabled}} mappings")

    cases = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise StudyError(
                f"{path}: pressure_cases[{i}] is a {type(entry).__name__}, "
                f"expected a mapping")
        case_keys = {"name", "pressure_psi", "enabled", "note"}
        extra = sorted(set(entry) - case_keys)
        if extra:
            raise StudyError(
                f"{path}: pressure_cases[{i}] has unknown key(s) {extra}; "
                f"valid: {sorted(case_keys)}")
        for required in ("name", "pressure_psi"):
            if required not in entry:
                raise StudyError(
                    f"{path}: pressure_cases[{i}] is missing {required!r}")
        pressure = float(entry["pressure_psi"])
        if pressure <= 0.0:
            raise StudyError(
                f"{path}: pressure_cases[{i}] has pressure_psi {pressure}; "
                f"a reservoir pressure must be positive")
        cases.append(PressureCase(
            name=str(entry["name"]),
            pressure_psi=pressure,
            enabled=bool(entry.get("enabled", True)),
            note=str(entry.get("note", "")),
        ))

    names = [c.name for c in cases]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise StudyError(
            f"{path}: duplicate case name(s) {duplicates}; each generates a "
            f"directory, so the second would overwrite the first")

    pressures = [c.pressure_psi for c in cases]
    duplicate_p = sorted({p for p in pressures if pressures.count(p) > 1})
    if duplicate_p:
        raise StudyError(
            f"{path}: duplicate pressure(s) {duplicate_p} psi; two cases at the "
            f"same pressure would differ in name only")

    return StudyConfig(
        base_case=str(raw.get("base_case", "baseCase")),
        cases_dir=str(raw.get("cases_dir", "Cases")),
        gap_in=float(raw.get("gap_in", 6.0)),
        gap_dir=str(raw.get("gap_dir", "gap06in")),
        pressure_cases=tuple(cases),
        manifest_name=str(raw.get("manifest_name", "manifest.yaml")),
    )


def filter_by_pressure(cases: list, wanted) -> list:
    """Keep only the cases whose pressure appears in ``wanted``.

    Args:
        cases: the candidate cases.
        wanted: pressures in psi. ``None`` or empty keeps everything.

    Raises:
        StudyError: if a requested pressure matches no case, so a typo on the
            command line stops the run instead of quietly generating nothing.
    """
    if not wanted:
        return list(cases)
    targets = {float(p) for p in wanted}
    kept = [c for c in cases if c.pressure_psi in targets]
    missing = sorted(targets - {c.pressure_psi for c in cases})
    if missing:
        raise StudyError(
            f"no case at {missing} psi; the study defines "
            f"{sorted(c.pressure_psi for c in cases)}")
    return kept


# --------------------------------------------------------------------------- #
# generation
# --------------------------------------------------------------------------- #

def require_casefoam():
    """Import and return ``casefoam``, or raise with the install command.

    A hard dependency, not a capability to probe. There is no built-in
    substitute: reimplementing a case generator this repository already depends
    on would mean maintaining two of them, and the one that is not exercised
    would be the one that drifts.
    """
    try:
        import casefoam
    except ImportError as exc:
        raise StudyError(
            "casefoam is not importable, and case generation requires it.\n"
            '    pip install -e ".[cases]"\n'
            "  See https://github.com/DLR-RY/caseFOAM"
        ) from exc
    return casefoam


def clone_cases(root: Path, study: StudyConfig, cases: list) -> list[Path]:
    """Build the case hierarchy with CaseFoam.

    Args:
        root: the study directory -- the one containing ``baseCase``.
        study: the study.
        cases: the pressure cases to generate, in order.

    Returns:
        The generated case directories, in the order given.

    Raises:
        StudyError: if CaseFoam is absent, the template is missing, or the
            expected hierarchy did not appear.

    ``casefoam.mkCases`` is called the way it is designed to be called::

        mkCases(<the baseCase directory>, [[gap], [case, ...]],
                caseData, hierarchy="tree", writeDir=<Cases>)

    It copies the template to ``writeDir``, moves that content down into
    ``writeDir/baseCase``, and creates ``writeDir/<gap>/<case>`` from it -- so
    the resulting ``Cases/baseCase`` alongside ``Cases/gap06in/p005psi`` is
    CaseFoam's own layout, the same one ``cases/caseFoamEx`` has. Its
    ``rmCases``, ``Allrun`` and ``Allclean`` land in ``Cases/`` too. None of that
    is cleaned up: it belongs to CaseFoam.

    ``caseData`` is empty because the physical values are applied afterwards by
    :func:`apply_case_parameters`. CaseFoam's own mechanism for a non-OpenFOAM
    file like ``case.yaml`` is ``'#!stringManipulation'``, which is the
    whitespace-sensitive substitution this design exists to avoid; its
    dictionary-aware path goes through PyFoam's parser, which does not read YAML.
    So CaseFoam does the cloning and the hierarchy, and the parameters are
    applied structurally.
    """
    casefoam = require_casefoam()

    root = Path(root)
    template = root / study.base_case
    if not template.is_dir():
        raise StudyError(f"base case {template} does not exist")
    if not (template / "case.yaml").is_file():
        raise StudyError(f"base case {template} has no case.yaml")

    structure = [[study.gap_dir], [c.name for c in cases]]
    data = {name: {} for level in structure for name in level}

    write_dir = root / study.cases_dir
    if write_dir.exists():
        # mkCases swallows FileExistsError from its copytree and would then
        # restructure whatever is already there.
        shutil.rmtree(write_dir)

    previous_cwd = Path.cwd()
    try:
        os.chdir(root)
        casefoam.mkCases(str(template), structure, data,
                         hierarchy="tree", writeDir=str(write_dir))
    finally:
        os.chdir(previous_cwd)

    generated = []
    for case in cases:
        destination = root / study.case_path(case)
        if not destination.is_dir():
            raise StudyError(
                f"casefoam did not produce {study.case_path(case)}. Expected the "
                f"'tree' hierarchy {study.cases_dir}/{study.gap_dir}/<case>.")
        if not (destination / "case.yaml").is_file():
            raise StudyError(
                f"{study.case_path(case)} has no case.yaml after cloning")
        prune_generated_dictionaries(destination)
        generated.append(destination)

    return generated


def prune_generated_dictionaries(case_dir: Path) -> list[Path]:
    """Remove dictionaries a case must generate from its own ``case.yaml``.

    Args:
        case_dir: a freshly cloned case.

    Returns:
        The paths removed.

    Template hygiene, not part of the cloning. CaseFoam copies the template
    faithfully -- as it should -- so a dictionary left in ``baseCase`` by someone
    running ``./Allmesh`` there would be inherited by every case. ``./Allmesh``
    would overwrite it, but ``./Allrun`` checks only that
    ``constant/dsmcProperties`` exists, so an inherited one would let a case run
    against the template's physics instead of its own.

    They are gitignored, so this only ever fires on a working copy where someone
    has meshed the template in place.
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
    return removed


def case_is_complete(case_dir: Path) -> bool:
    """Whether a case looks like it has been run.

    True when the case has a time directory other than ``0``, which only the
    solver creates. Used to refuse to overwrite results: regenerating a case
    rewrites its ``case.yaml``, and doing that under finished output would leave
    results whose inputs no longer describe them.
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


def apply_case_parameters(case_dir: Path, case: PressureCase,
                          study: StudyConfig,
                          n_equivalent_particles: float | None = None) -> dict:
    """Rewrite a case-local ``case.yaml`` with this case's parameters.

    Args:
        case_dir: the generated case directory.
        case: the pressure case.
        study: the study, for the gap metadata.
        n_equivalent_particles: the particle weight to record. ``None`` leaves
            the existing value, which means "derive at mesh time".

    Returns:
        The parameters written, for the manifest.

    Raises:
        StudyError: if the file is missing or does not have the sections this
            edits -- rather than adding them, because a ``case.yaml`` without a
            ``stagnation`` section is not this case family's template and
            silently growing one would hide that.

    **Structured, not textual.** The file is parsed to a dict, specific keys are
    set, and it is re-emitted. A key that is not there is an error; a key that is
    there is replaced whatever its formatting.

    The cost is that comments do not survive round-tripping through PyYAML, so a
    generated ``case.yaml`` has none. That is the right trade for a *generated*
    file -- ``baseCase/case.yaml`` keeps its comments and is the one anybody
    edits -- and the manifest records every value anyway.
    """
    path = Path(case_dir) / "case.yaml"
    if not path.is_file():
        raise StudyError(f"no case.yaml at {path}")

    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    for section in ("stagnation", "dsmc", "meta"):
        if section not in data:
            raise StudyError(
                f"{path}: no {section!r} section. This does not look like the "
                f"markelov1999 template; generation would silently add keys the "
                f"rest of the file does not expect.")

    applied = {
        "p0_pa": case.pressure_pa,
        "pressure_psi": case.pressure_psi,
    }
    data["stagnation"]["p0_pa"] = case.pressure_pa

    if n_equivalent_particles is not None:
        data["dsmc"]["n_equivalent_particles"] = float(n_equivalent_particles)
        applied["n_equivalent_particles"] = float(n_equivalent_particles)

    data["meta"] = dict(data["meta"])
    data["meta"].update({
        "case_name": case.name,
        "pressure_psi": case.pressure_psi,
        "pressure_pa": case.pressure_pa,
        "gap_in": study.gap_in,
        "generated_by": "cases/markelov1999/generate_cases.py",
    })
    if case.note:
        data["meta"]["note"] = case.note

    header = (
        "# GENERATED by cases/markelov1999/generate_cases.py -- do not edit.\n"
        "#\n"
        f"# {case.name}: reservoir pressure {case.pressure_psi:g} psi "
        f"= {case.pressure_pa:.6f} Pa.\n"
        "#\n"
        "# Edit ../../../baseCase/case.yaml and regenerate. The comments in that\n"
        "# file explain every value; they do not survive the YAML round trip, and\n"
        "# manifest.yaml records what was applied here.\n"
        "\n"
    )
    path.write_text(
        header + yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return applied


def casefoam_version() -> str:
    """The installed CaseFoam version, for the manifest.

    Recorded because the generated tree is CaseFoam's output: if its layout ever
    changes, the manifest says which version produced what is on disk.
    """
    try:
        from importlib.metadata import version
        return version("casefoam")
    except Exception:  # pragma: no cover - metadata is present in any install
        return "unknown"


def write_manifest(path: Path, study: StudyConfig, entries: list) -> Path:
    """Write the generation manifest.

    Args:
        path: destination file.
        study: the study.
        entries: one dict per generated case.

    Returns:
        The path written.

    The manifest is the auditable map from case name to physical inputs that
    requirement 2 asks for: given a results directory, it says exactly what was
    generated, from what, with which pressure and which particle weight. It also
    records the cases that were *skipped*, so a disabled case is visible as a
    decision rather than as an absence.
    """
    import datetime

    path = Path(path)
    document = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "reference": "AIAA 99-3455",
        "generator": f"casefoam {casefoam_version()}",
        "base_case": study.base_case,
        "cases_dir": study.cases_dir,
        "gap_in": study.gap_in,
        "gap_dir": study.gap_dir,
        "n_generated": len(entries),
        "cases": entries,
        "disabled_cases": [
            {"name": c.name, "pressure_psi": c.pressure_psi, "note": c.note}
            for c in study.all_cases() if not c.enabled
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by cases/markelov1999/generate_cases.py -- do not edit.\n"
        "# The auditable map from case name to physical inputs.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path
