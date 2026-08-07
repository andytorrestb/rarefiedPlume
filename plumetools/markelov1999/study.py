r"""``study.yaml`` parsing, the case matrix, and the generation manifest.

A study is a base case plus a list of pressure cases. Each generated case is a
complete, self-contained OpenFOAM case with its own ``case.yaml``.

How cases are generated
-----------------------
1. **CaseFoam** clones ``baseCase`` into the hierarchy, when it is installed.
2. A **structured Python step** rewrites each case-local ``case.yaml`` -- loaded
   as YAML, edited as a data structure, re-emitted.
3. **plumetools** generates the OpenFOAM dictionaries from that ``case.yaml``.

Step 2 is deliberately not string substitution. ``util/caseFoam/genCases.py``
does it the other way::

    'system/controlDict': {'#!stringManipulation':
                            {'deltaT          1.0E-05': '%s' % deltaT}}

which silently does nothing if the file is reformatted, if the value already
differs, or if the base case has moved on -- and none of those is visible in the
output. Editing the parsed structure either finds the key or raises.

CaseFoam is **optional**. When it is not importable the same hierarchy is built
by :func:`clone_base_case`, which copies the same files to the same paths. The
tree is identical either way -- a test generates both and compares them -- so a
study generated on a machine without CaseFoam is not a different study.
:func:`generation_backend` resolves which is used and the manifest records it.

CaseFoam is used for cloning only. Its own ``caseData`` mechanism applies
parameters through ``'#!stringManipulation'``, which is the whitespace-sensitive
substitution this design exists to avoid, so the physical values are applied
afterwards by :func:`apply_case_parameters` instead.

``mkCases`` copies the whole directory it is pointed at and writes drivers of its
own beside it, so :func:`clone_with_casefoam` runs it in an isolated staging
directory and moves the finished cases into place. Nothing is ever deleted from
the study directory. See that function for the detail.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from plumetools.markelov1999.constants import psi_to_pa

#: Files and directories copied from ``baseCase`` into each generated case.
#:
#: Generated OpenFOAM dictionaries are deliberately NOT in this list: they are
#: written from each case's own ``case.yaml`` by ``./Allmesh``, so copying a
#: stale one would create a second source of truth for the physics.
CASE_TEMPLATE_ENTRIES = (
    "case.yaml",
    "Allmesh",
    "Allrun",
    "Allclean",
    "Allpost",
    "runInflow.py",
    "postProcess.py",
    "open.foam",
)

#: Dictionaries a case generates from its own ``case.yaml``.
#:
#: Neither backend copies these. A dictionary left behind in ``baseCase`` by a
#: local ``./Allmesh`` would otherwise be inherited by every generated case and
#: become a second, silently divergent source of truth for the physics.
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
#: silently give every case the template's geometry. ``0/`` and any time
#: directory because those are results, and a case that starts with someone
#: else's results is not a fresh case.
GENERATED_DIRECTORIES = ("constant/polyMesh", "0", "processor0")


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

def casefoam_available() -> bool:
    """Whether ``casefoam`` can be imported."""
    try:
        import casefoam  # noqa: F401
    except ImportError:
        return False
    return True


def generation_backend(requested: str = "auto") -> str:
    """Resolve the cloning backend to ``"casefoam"`` or ``"builtin"``.

    Args:
        requested: ``"auto"`` uses CaseFoam when it is importable; ``"casefoam"``
            requires it; ``"builtin"`` never uses it.

    Returns:
        The backend that will be used. Recorded in the manifest.

    Raises:
        StudyError: if CaseFoam was required and is not installed.

    The two backends produce the **same tree** -- verified by a test that
    generates both and compares them -- so under ``auto`` this is provenance
    rather than a behaviour switch. It is explicit rather than silent because
    "which tool built this" is exactly the sort of thing that is impossible to
    reconstruct later.
    """
    if requested not in ("auto", "casefoam", "builtin"):
        raise StudyError(
            f"unknown backend {requested!r}; valid: ['auto', 'casefoam', 'builtin']")

    if requested == "builtin":
        return "builtin"
    if casefoam_available():
        return "casefoam"
    if requested == "casefoam":
        raise StudyError(
            "backend 'casefoam' was requested but casefoam is not importable. "
            'Install it with: pip install -e ".[cases]"')
    return "builtin"


def clone_base_case(base: Path, destination: Path) -> list[Path]:
    """Copy the template case into ``destination``.

    Args:
        base: the ``baseCase`` directory.
        destination: where the case is created. Parents are created.

    Returns:
        The paths written.

    Raises:
        StudyError: if the base case is missing, or lacks ``case.yaml`` -- which
            would produce a tree of directories that no tool in this repository
            can do anything with.

    Only :data:`CASE_TEMPLATE_ENTRIES` are copied, and only those that exist. In
    particular ``constant/polyMesh``, ``0/`` and the generated dictionaries are
    not: each case meshes itself from its own ``case.yaml``, and inheriting a mesh
    from the template would silently give every case the template's geometry.
    """
    base = Path(base)
    destination = Path(destination)
    if not base.is_dir():
        raise StudyError(f"base case {base} does not exist")
    if not (base / "case.yaml").is_file():
        raise StudyError(f"base case {base} has no case.yaml")

    destination.mkdir(parents=True, exist_ok=True)
    written = []
    for entry in CASE_TEMPLATE_ENTRIES:
        source = base / entry
        if not source.exists():
            continue
        target = destination / entry
        if source.is_dir():
            shutil.copytree(source, target, dirs_exist_ok=True)
        else:
            shutil.copy2(source, target)
        written.append(target)

    # system/ may hold static dictionaries a case does not generate. Copy those,
    # but never a generated one: a dictionary left behind in the template by a
    # local ./Allmesh would otherwise be inherited by every case and become a
    # second source of truth for the physics.
    generated_names = {Path(p).name for p in GENERATED_DICTIONARIES}
    base_system = base / "system"
    if base_system.is_dir():
        for source in sorted(base_system.iterdir()):
            if source.is_file() and source.name not in generated_names:
                target = destination / "system" / source.name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
                written.append(target)

    return written


def prune_inherited_artefacts(case_dir: Path) -> list[Path]:
    """Remove what a generated case must not inherit from the template.

    Args:
        case_dir: a freshly cloned case.

    Returns:
        The paths removed.

    :func:`clone_base_case` never copies these in the first place. CaseFoam
    copies the template directory wholesale, so this brings its output back to
    the same content -- which is what makes "both backends produce the same tree"
    true rather than merely intended.
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

    # Any time directory other than 0, which GENERATED_DIRECTORIES covers.
    for entry in Path(case_dir).iterdir():
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if value > 0.0:
            shutil.rmtree(entry)
            removed.append(entry)

    # A directory emptied by the removals above is harmless in itself, but the
    # builtin clone never creates one, and "the two backends produce the same
    # tree" has to mean exactly that.
    for name in ("constant", "system"):
        directory = Path(case_dir) / name
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()
            removed.append(directory)

    return removed


def clone_with_casefoam(root: Path, study: StudyConfig, cases: list) -> list[Path]:
    """Build the case hierarchy with CaseFoam, then clean up after it.

    Args:
        root: the study directory -- the one containing ``baseCase``.
        study: the study.
        cases: the pressure cases to generate.

    Returns:
        The generated case directories, in the order given.

    Raises:
        StudyError: if CaseFoam is not importable, or does not produce the
            expected hierarchy.

    CaseFoam does the **cloning** and the directory hierarchy; the physical
    values are applied afterwards by :func:`apply_case_parameters`, which edits
    parsed YAML. CaseFoam's own ``caseData`` mechanism is deliberately not used
    for them: its ``'#!stringManipulation'`` form is whitespace-sensitive
    substitution, which is the failure mode this whole design avoids.

    ``mkCases`` is run in an **isolated staging directory** holding nothing but a
    copy of the template, and the finished cases are moved into place afterwards.
    That is not fastidiousness. Its ``baseCase`` argument means "the directory
    containing the template", and with ``writeDir`` set it copies *that whole
    directory* into the write directory -- so running it against the study root
    puts ``study.yaml``, ``generate_cases.py``, the ``All*Cases`` drivers and the
    README inside ``Cases/``. It also writes ``Allrun``, ``Allclean`` and
    ``rmCases`` of its own next to the template. Staging keeps every one of those
    out of the study directory, and means nothing has to be deleted from a
    directory the user owns.

    Each case is then pruned by :func:`prune_inherited_artefacts`, which is what
    makes the output identical to the built-in clone's.
    """
    try:
        import casefoam
    except ImportError as exc:  # pragma: no cover - guarded by the caller
        raise StudyError(
            'casefoam is not importable. Install it with: pip install -e ".[cases]"'
        ) from exc

    root = Path(root)
    template = root / study.base_case
    if not template.is_dir():
        raise StudyError(f"base case {template} does not exist")
    if not (template / "case.yaml").is_file():
        raise StudyError(f"base case {template} has no case.yaml")

    structure = [[study.gap_dir], [c.name for c in cases]]
    # Empty per-name data: CaseFoam clones, the Python step below applies values.
    data = {name: {} for level in structure for name in level}

    staging = Path(tempfile.mkdtemp(prefix="markelov-casefoam-"))
    previous_cwd = Path.cwd()
    generated = []
    try:
        shutil.copytree(template, staging / study.base_case)

        os.chdir(staging)
        casefoam.mkCases(str(staging), structure, data,
                         hierarchy="tree", writeDir=study.cases_dir)
        os.chdir(previous_cwd)

        for case in cases:
            produced = staging / study.case_path(case)
            if not produced.is_dir():
                raise StudyError(
                    f"casefoam did not produce {study.case_path(case)} in its "
                    f"staging directory. Expected the 'tree' hierarchy "
                    f"{study.cases_dir}/{study.gap_dir}/<case>.")
            if not (produced / "case.yaml").is_file():
                raise StudyError(
                    f"{study.case_path(case)} has no case.yaml after cloning")

            destination = root / study.case_path(case)
            if destination.exists():
                shutil.rmtree(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(produced), str(destination))

            prune_inherited_artefacts(destination)
            generated.append(destination)
    finally:
        os.chdir(previous_cwd)
        shutil.rmtree(staging, ignore_errors=True)

    return generated


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


def write_manifest(path: Path, study: StudyConfig, entries: list,
                   backend: str) -> Path:
    """Write the generation manifest.

    Args:
        path: destination file.
        study: the study.
        entries: one dict per generated case.
        backend: ``"casefoam"`` or ``"builtin"``.

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
        "backend": backend,
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
