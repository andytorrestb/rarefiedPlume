r"""The statistical-health sweep: a 3 x 3 matrix over parcels and averaging.

``cases/cai2012`` measures how well DSMC reproduces Cai's collisionless
solution. It never measures how much **statistical noise** its own answer
carries. :func:`plumetools.cai2012.checks._check_occupancy` prints an *estimated*
exit-cell occupancy and says outright that only a post-run audit of the sampled
``dsmcRhoN`` field can say what it actually was. This family is that audit.

Two knobs set statistical quality, and this sweeps both at fixed Kn = 100 -- the
most rarefied case in that study, and the one whose imagery is visibly speckled:

.. code-block:: text

    resolution.target_particles_per_cell   how many parcels a cell holds
    dsmc.sampling_domain_transits          how long fieldAverage accumulates

Everything else is held fixed -- the same mesh, the same Knudsen number, the same
time step, the same transient -- so every case is comparable image to image.
That is not a convention here, it is enforced: this module's ``study.yaml`` has
**no override mechanism at all**. The two axes are the entire per-case input.
Anything else worth changing is changed in ``baseCase/case.yaml``, where it
necessarily applies to all nine and the matrix stays comparable by construction.

Why this one uses CaseFoam and ``cases/cai2012`` does not
--------------------------------------------------------
:mod:`plumetools.cai2012.study` says why that family does not: it is a **flat
list** of three Knudsen numbers, and requiring a dependency to do a
``copytree`` would make ``generate_cases.py`` fail on a machine where the core
install works.

This family is a **hierarchy** -- weight over sampling duration -- and building
hierarchies is what CaseFoam is for. It earns the dependency exactly the way
``cases/markelov1999``'s ``<gap>/<pressure>`` does, and
:func:`plumetools.markelov1999.study.clone_cases` is the call this follows.
That function cannot be reused directly: its level 1 is a single ``gap_dir``,
and it ``rmtree``\ s the write directory on entry, so three calls to build three
weight rows would each destroy the last. Its ``require_casefoam`` is reused,
and so is ``mkCases(template, [[...], [...]], data, hierarchy="tree")`` itself.

Costing before running
----------------------
The matrix is expensive and the expense is wildly uneven: run time goes as
parcels x total transits, so the top weight row is most of the bill. The cost
model here is fitted from cases that have actually run
(:func:`fit_cost_model`) and is printed by ``./generate_cases.py --dry-run``
before anything is committed to. Its two terms are separated deliberately --
per-parcel work scales with the weight axis and per-cell work does not, so a
single measurement extrapolates badly down the axis.
"""

from __future__ import annotations

import copy
import datetime
import math
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from plumetools.cai2012.study import (
    StudyError,
    case_is_complete,
    prune_generated_dictionaries,
)

#: The Knudsen number this family is defined at.
#:
#: Not a per-case input and not overridable: the sweep exists to explain the
#: speckle in ``cases/cai2012/Cases/Kn100``'s imagery, and an image taken at
#: another Knudsen number cannot be compared with it. ``baseCase/case.yaml``
#: carries it, and :func:`load_study` checks the template still says so.
FAMILY_KNUDSEN = 100.0

#: ``case.yaml`` values the sweep writes, per case. Exactly two, one per axis.
SWEPT_KEYS = (
    "resolution.target_particles_per_cell",
    "dsmc.sampling_domain_transits",
)


@dataclass(frozen=True)
class WeightLevel:
    """One rung of the particle-weight axis -- a row of the matrix.

    Attributes:
        name: directory name, e.g. ``ppc020``.
        particles_per_cell: ``resolution.target_particles_per_cell``. The
            particle weight is derived from it so the **exit cell** -- the
            densest in the domain -- reaches this occupancy; everywhere else
            follows from the density ratio. Standard ``dsmcFoam`` has one
            global weight (``DSMCCloud::nParticle_`` is a single scalar), so
            the target can be met in exactly one place.
        enabled: whether :meth:`HealthStudy.enabled_cases` includes the row.
        note: free text, carried into the manifest.
    """

    name: str
    particles_per_cell: float
    enabled: bool = True
    note: str = ""


@dataclass(frozen=True)
class SamplingLevel:
    """One rung of the averaging-duration axis -- a column of the matrix.

    Attributes:
        name: directory name, e.g. ``s1p5``.
        domain_transits: ``dsmc.sampling_domain_transits``; how long
            ``fieldAverage`` accumulates, in transits of the domain at ``U0``.
        enabled: whether :meth:`HealthStudy.enabled_cases` includes the column.
        note: free text, carried into the manifest.
    """

    name: str
    domain_transits: float
    enabled: bool = True
    note: str = ""


@dataclass(frozen=True)
class HealthCase:
    """One cell of the matrix: a weight row crossed with a sampling column."""

    weight: WeightLevel
    sampling: SamplingLevel

    @property
    def name(self) -> str:
        """``ppc020/s1p5`` -- the path, which is also how a case is named."""
        return f"{self.weight.name}/{self.sampling.name}"

    @property
    def particles_per_cell(self) -> float:
        return self.weight.particles_per_cell

    @property
    def domain_transits(self) -> float:
        return self.sampling.domain_transits


@dataclass(frozen=True)
class HealthStudy:
    """A parsed ``study.yaml`` for the health family.

    Attributes:
        base_case: the template, relative to ``study.yaml``.
        cases_dir: where the hierarchy is written.
        weights: every weight row, enabled or not.
        samplings: every sampling column, enabled or not.
        manifest_name: filename of the generation manifest.
        meta: free-form metadata, copied into the manifest.
    """

    base_case: str = "baseCase"
    cases_dir: str = "Cases"
    weights: tuple = ()
    samplings: tuple = ()
    manifest_name: str = "manifest.yaml"
    meta: dict = field(default_factory=dict)

    def enabled_weights(self) -> list:
        """Enabled rows, cheapest first.

        Sorted by occupancy rather than by file order, so reordering
        ``study.yaml`` cannot change the generated tree -- and so the matrix is
        generated in increasing cost, which is the order anyone wants to run it
        in when the top row is most of the bill.
        """
        return sorted((w for w in self.weights if w.enabled),
                      key=lambda w: (w.particles_per_cell, w.name))

    def enabled_samplings(self) -> list:
        """Enabled columns, shortest first, for the same reasons."""
        return sorted((s for s in self.samplings if s.enabled),
                      key=lambda s: (s.domain_transits, s.name))

    def enabled_cases(self) -> list[HealthCase]:
        """The matrix, cheapest first: whole weight rows in increasing cost."""
        return [HealthCase(weight=w, sampling=s)
                for w in self.enabled_weights()
                for s in self.enabled_samplings()]

    def all_cases(self) -> list[HealthCase]:
        """Every cell, enabled or not, in the same deterministic order."""
        weights = sorted(self.weights, key=lambda w: (w.particles_per_cell, w.name))
        samplings = sorted(self.samplings, key=lambda s: (s.domain_transits, s.name))
        return [HealthCase(weight=w, sampling=s)
                for w in weights for s in samplings]

    def case_path(self, case: HealthCase) -> Path:
        """``Cases/ppc020/s1p5`` -- relative to ``study.yaml``."""
        return Path(self.cases_dir) / case.weight.name / case.sampling.name

    @property
    def baseline(self) -> HealthCase | None:
        """The cell whose settings match ``cases/cai2012/Cases/Kn100``.

        20 parcels per cell and 1.5 transits of sampling: the same run, at a
        finer write interval. Its centreline errors must reproduce that case's
        published numbers, which makes it the control for the whole sweep --
        and it is where the pinned colour ranges come from, since nine cases
        auto-scaled to nine different ranges cannot be compared at all.
        """
        for case in self.enabled_cases():
            if (math.isclose(case.particles_per_cell, 20.0)
                    and math.isclose(case.domain_transits, 1.5)):
                return case
        return None


def load_study(path: Path) -> HealthStudy:
    """Read and validate the health family's ``study.yaml``.

    Args:
        path: the ``study.yaml`` file, or the directory containing it.

    Returns:
        A validated :class:`HealthStudy`.

    Raises:
        StudyError: for a missing file, an unknown key, a duplicate or
            non-positive level, an empty axis, or an ``overrides`` block --
            which this family deliberately does not have.
    """
    path = Path(path)
    if path.is_dir():
        path = path / "study.yaml"
    if not path.is_file():
        raise StudyError(f"no study.yaml at {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise StudyError(f"{path}: expected a mapping at the top level")

    known = {"base_case", "cases_dir", "weights", "samplings", "manifest_name",
             "meta"}
    unknown = sorted(set(raw) - known)
    if unknown:
        if "overrides" in unknown or "kn_cases" in unknown:
            raise StudyError(
                f"{path}: this study has no {sorted(set(unknown) & {'overrides', 'kn_cases'})} "
                f"mechanism. The two axes are the whole per-case input, because "
                f"everything else has to be IDENTICAL across the matrix for the "
                f"images to be comparable. Change baseCase/case.yaml instead -- "
                f"that applies to all of them at once, which is the point.")
        raise StudyError(
            f"{path}: unknown top-level key(s) {unknown}; valid: {sorted(known)}")

    weights = tuple(_load_weights(raw.get("weights"), path))
    samplings = tuple(_load_samplings(raw.get("samplings"), path))

    study = HealthStudy(
        base_case=str(raw.get("base_case", "baseCase")),
        cases_dir=str(raw.get("cases_dir", "Cases")),
        weights=weights,
        samplings=samplings,
        manifest_name=str(raw.get("manifest_name", "manifest.yaml")),
        meta=raw.get("meta") or {},
    )
    if not study.enabled_weights() or not study.enabled_samplings():
        raise StudyError(
            f"{path}: every level of an axis is disabled, so the matrix is "
            f"empty. A disabled level stays in the file as the record of the "
            f"intended matrix; disabling a whole axis is not a sweep.")
    return study


def _load_weights(entries, path) -> list[WeightLevel]:
    """Parse and validate the particle-weight axis."""
    if not isinstance(entries, list) or not entries:
        raise StudyError(
            f"{path}: weights must be a non-empty list of "
            f"{{name, target_particles_per_cell}} mappings")
    keys = {"name", "target_particles_per_cell", "enabled", "note"}
    levels = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise StudyError(f"{path}: weights[{i}] is a {type(entry).__name__}, "
                             f"expected a mapping")
        extra = sorted(set(entry) - keys)
        if extra:
            raise StudyError(f"{path}: weights[{i}] has unknown key(s) {extra}; "
                             f"valid: {sorted(keys)}")
        for required in ("name", "target_particles_per_cell"):
            if required not in entry:
                raise StudyError(f"{path}: weights[{i}] is missing {required!r}")
        value = float(entry["target_particles_per_cell"])
        if not value > 0.0:
            raise StudyError(
                f"{path}: weights[{i}] asks for {value} parcels per cell. The "
                f"particle weight is n0 * V_exit / this, so a non-positive "
                f"target is a non-positive or infinite weight.")
        levels.append(WeightLevel(
            name=str(entry["name"]),
            particles_per_cell=value,
            enabled=bool(entry.get("enabled", True)),
            note=str(entry.get("note", "")),
        ))
    _reject_duplicates(levels, "weights", lambda w: w.particles_per_cell,
                       "parcels per cell", path)
    return levels


def _load_samplings(entries, path) -> list[SamplingLevel]:
    """Parse and validate the averaging-duration axis."""
    if not isinstance(entries, list) or not entries:
        raise StudyError(
            f"{path}: samplings must be a non-empty list of "
            f"{{name, sampling_domain_transits}} mappings")
    keys = {"name", "sampling_domain_transits", "enabled", "note"}
    levels = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise StudyError(f"{path}: samplings[{i}] is a "
                             f"{type(entry).__name__}, expected a mapping")
        extra = sorted(set(entry) - keys)
        if extra:
            raise StudyError(f"{path}: samplings[{i}] has unknown key(s) "
                             f"{extra}; valid: {sorted(keys)}")
        for required in ("name", "sampling_domain_transits"):
            if required not in entry:
                raise StudyError(f"{path}: samplings[{i}] is missing {required!r}")
        value = float(entry["sampling_domain_transits"])
        if not value > 0.0:
            raise StudyError(
                f"{path}: samplings[{i}] averages for {value} domain transits. "
                f"A run that ends before it starts sampling produces no *Mean "
                f"field at all, and there would be nothing to measure.")
        levels.append(SamplingLevel(
            name=str(entry["name"]),
            domain_transits=value,
            enabled=bool(entry.get("enabled", True)),
            note=str(entry.get("note", "")),
        ))
    _reject_duplicates(levels, "samplings", lambda s: s.domain_transits,
                       "sampling duration", path)
    return levels


def _reject_duplicates(levels, axis, value_of, what, path) -> None:
    """No two levels of an axis may share a name or a value."""
    names = [level.name for level in levels]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise StudyError(
            f"{path}: duplicate {axis} name(s) {duplicates}; each names a "
            f"directory, so the second would overwrite the first")
    values = [value_of(level) for level in levels]
    repeated = sorted({v for v in values if values.count(v) > 1})
    if repeated:
        raise StudyError(
            f"{path}: duplicate {what} {repeated} in {axis}; two levels at the "
            f"same value would differ in name only, and the matrix would have a "
            f"column that measures nothing")


# --------------------------------------------------------------------------- #
# generation
# --------------------------------------------------------------------------- #

def clone_cases(root: Path, study: HealthStudy, cases: list) -> list[Path]:
    """Build the ``<weight>/<sampling>`` hierarchy with CaseFoam.

    Args:
        root: the study directory -- the one containing ``baseCase``.
        study: the study.
        cases: the cells to generate. Their axes decide the tree, so a
            partial selection still produces a well-formed hierarchy.

    Returns:
        The generated case directories, in the order given.

    Raises:
        StudyError: if CaseFoam is absent, the template is missing, or the
            expected hierarchy did not appear.

    ``casefoam.mkCases`` is called the way it is designed to be called, and the
    way :func:`plumetools.markelov1999.study.clone_cases` calls it::

        mkCases(<baseCase>, [[ppc005, ppc020, ppc040], [s0p5, s1p5, s4p5]],
                caseData, hierarchy="tree", writeDir=<Cases>)

    The difference from that function, and the reason this one exists, is that
    its first level is a **list**: markelov's is the single ``gap_dir``, so
    three calls would be needed here and each ``rmtree``\ s the write directory
    on entry, destroying the row before it.

    ``caseData`` is empty for the same reason it is there: CaseFoam's mechanism
    for a non-OpenFOAM file like ``case.yaml`` is ``'#!stringManipulation'``,
    the whitespace-sensitive substitution both families exist to avoid, and its
    dictionary-aware path goes through PyFoam's parser, which does not read
    YAML. So CaseFoam does the cloning and the hierarchy, and
    :func:`apply_case_parameters` applies the values structurally.

    CaseFoam also leaves its own ``rmCases``, ``Allrun`` and ``Allclean`` in
    ``Cases/``, along with a ``Cases/baseCase``. None of that is cleaned up: it
    belongs to CaseFoam, and ``Cases/`` is a gitignored build product.
    """
    from plumetools.markelov1999.study import require_casefoam

    casefoam = require_casefoam()

    root = Path(root)
    template = root / study.base_case
    if not template.is_dir():
        raise StudyError(f"base case {template} does not exist")
    if not (template / "case.yaml").is_file():
        raise StudyError(f"base case {template} has no case.yaml")

    structure = [
        sorted({case.weight.name for case in cases}),
        sorted({case.sampling.name for case in cases}),
    ]
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
                f"casefoam did not produce {study.case_path(case)}. Expected "
                f"the 'tree' hierarchy {study.cases_dir}/<weight>/<sampling>.")
        if not (destination / "case.yaml").is_file():
            raise StudyError(
                f"{study.case_path(case)} has no case.yaml after cloning")
        prune_generated_dictionaries(destination)
        generated.append(destination)

    return generated


def apply_case_parameters(case_dir: Path, case: HealthCase,
                          study: HealthStudy) -> dict:
    """Rewrite a case-local ``case.yaml`` with this cell's two swept values.

    Args:
        case_dir: the generated case directory.
        case: the matrix cell.
        study: the study, for the metadata.

    Returns:
        The parameters written, for the manifest.

    Raises:
        StudyError: if the file is missing, is not a cai2012 template, or is
            not at :data:`FAMILY_KNUDSEN`.

    **Structured, not textual**, as in both other families: the file is parsed
    to a dict, two keys are set, and it is re-emitted. Comments do not survive
    the round trip, which is the right trade for a generated file --
    ``baseCase/case.yaml`` keeps them and is the one anybody edits.

    Exactly two keys are written. Not the Knudsen number, not the mesh, not the
    time step, not the transient: those are what "everything else held fixed"
    means, and writing them per case is how a matrix stops being comparable.
    """
    path = Path(case_dir) / "case.yaml"
    if not path.is_file():
        raise StudyError(f"no case.yaml at {path}")

    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    for section in ("exit", "dsmc", "resolution", "meta"):
        if section not in data:
            raise StudyError(
                f"{path}: no {section!r} section. This does not look like the "
                f"cai2012 template; generation would silently add keys the rest "
                f"of the file does not expect.")

    knudsen = float(data["exit"].get("knudsen", float("nan")))
    if not math.isclose(knudsen, FAMILY_KNUDSEN, rel_tol=1e-12):
        raise StudyError(
            f"{path}: exit.knudsen is {knudsen:g}, not {FAMILY_KNUDSEN:g}. This "
            f"family sweeps statistical quality at ONE Knudsen number -- the "
            f"most rarefied case in cases/cai2012, whose imagery this exists to "
            f"explain. At another Kn the density, the parcel count and the "
            f"picture all change, and nothing in the matrix would be comparable "
            f"with that case.")

    data["resolution"]["target_particles_per_cell"] = float(
        case.particles_per_cell)
    data["dsmc"]["sampling_domain_transits"] = float(case.domain_transits)
    applied = {
        "target_particles_per_cell": float(case.particles_per_cell),
        "sampling_domain_transits": float(case.domain_transits),
    }

    data["meta"] = dict(data.get("meta") or {})
    data["meta"].update({
        "case_name": case.name,
        "weight_level": case.weight.name,
        "sampling_level": case.sampling.name,
        "target_particles_per_cell": float(case.particles_per_cell),
        "sampling_domain_transits": float(case.domain_transits),
        "generated_by": "cases/cai2012-health/generate_cases.py",
    })
    notes = [n for n in (case.weight.note, case.sampling.note) if n]
    if notes:
        data["meta"]["note"] = " ".join(notes)

    header = (
        "# GENERATED by cases/cai2012-health/generate_cases.py -- do not edit.\n"
        "#\n"
        f"# {case.name}: {case.particles_per_cell:g} parcels per exit cell, "
        f"{case.domain_transits:g} domain transits of averaging.\n"
        "#\n"
        f"# Kn = {FAMILY_KNUDSEN:g}, the mesh, the time step and the transient "
        "are the SAME in\n"
        "# every case of this matrix. Only the two values above vary, which is "
        "what makes\n"
        "# the nine cases comparable image to image.\n"
        "#\n"
        "# Edit ../../../baseCase/case.yaml and regenerate. The comments in "
        "that file\n"
        "# explain every value; they do not survive the YAML round trip, and\n"
        "# manifest.yaml records what was applied here.\n"
        "\n"
    )
    path.write_text(
        header + yaml.safe_dump(data, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return applied


# --------------------------------------------------------------------------- #
# cost
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class CostModel:
    """Wall-clock seconds for a run, as two separable terms.

    .. code-block:: text

        seconds = (per_parcel_s * parcels + per_cell_s * cells) * steps

    Attributes:
        per_parcel_s: seconds per parcel per step -- move, index, collide.
        per_cell_s: seconds per cell per step -- the loop over cells that runs
            whether or not there is anything in them.
        source: where the coefficients came from, for the report.
        n_measurements: how many runs were fitted.

    **Why two terms.** The weight axis multiplies the parcel work by sixteen and
    leaves the cell work alone. A single-term model calibrated on the baseline
    would predict the ppc005 row at a quarter of the baseline's cost, when in
    practice the cell loop puts a floor under it -- so the cheap row looks
    cheaper than it is, and the extrapolation is wrong in the direction that
    matters least. Fitted the other way round, from the cheap row upward, a
    one-term model *understates* the expensive row, which is the direction that
    matters most.
    """

    per_parcel_s: float
    per_cell_s: float
    source: str = "prior"
    n_measurements: int = 0

    def seconds(self, *, parcels: float, cells: float, steps: float) -> float:
        """Predicted wall clock for one run [s]."""
        return (self.per_parcel_s * float(parcels)
                + self.per_cell_s * float(cells)) * float(steps)

    def as_dict(self) -> dict:
        return {
            "per_parcel_s": self.per_parcel_s,
            "per_cell_s": self.per_cell_s,
            "source": self.source,
            "n_measurements": self.n_measurements,
        }


#: Cost model before anything has been measured.
#:
#: From the completed ``cases/cai2012/Cases/Kn100`` run: 1.31e6 cells, ~1.2e6
#: parcels, 7212 steps, about an hour on 12 ranks (``cases/cai2012/README.md``).
#: One measurement cannot separate the two terms, so the split is assumed --
#: 80% of the work on the parcels -- and it is labelled ``prior`` everywhere it
#: is printed. Two runs replace it with a fit; see :func:`fit_cost_model`.
PRIOR_COST_MODEL = CostModel(
    per_parcel_s=0.8 * 3600.0 / (1.2e6 * 7212.0),
    per_cell_s=0.2 * 3600.0 / (1.31e6 * 7212.0),
    source="prior (cases/cai2012 Kn100, ~1 h on 12 ranks)",
    n_measurements=0,
)


def fit_cost_model(measurements: list) -> CostModel:
    """Least-squares fit of :class:`CostModel` to runs that have happened.

    Args:
        measurements: dicts with ``seconds``, ``parcels``, ``cells`` and
            ``steps``. Anything else in them is ignored, so a run record can
            carry its case name and settings too.

    Returns:
        The fitted model, or :data:`PRIOR_COST_MODEL` if there is nothing
        usable to fit.

    With one measurement the two coefficients are not separable, so the prior's
    *ratio* between them is kept and only the scale is fitted. That is stated in
    ``source`` rather than presented as a two-parameter fit of one number.
    """
    usable = [m for m in measurements
              if float(m.get("seconds", 0.0)) > 0.0
              and float(m.get("steps", 0.0)) > 0.0]
    if not usable:
        return PRIOR_COST_MODEL

    if len(usable) == 1:
        one = usable[0]
        predicted = PRIOR_COST_MODEL.seconds(
            parcels=one["parcels"], cells=one["cells"], steps=one["steps"])
        scale = float(one["seconds"]) / predicted if predicted > 0.0 else 1.0
        return CostModel(
            per_parcel_s=PRIOR_COST_MODEL.per_parcel_s * scale,
            per_cell_s=PRIOR_COST_MODEL.per_cell_s * scale,
            source="scaled from 1 measured run (the two terms are not "
                   "separable from one point; the prior's ratio is kept)",
            n_measurements=1,
        )

    # Ordinary least squares on seconds = a*(parcels*steps) + b*(cells*steps).
    # numpy rather than a hand-rolled normal equation: the two columns differ by
    # orders of magnitude and lstsq is the one that stays conditioned.
    import numpy as np

    design = np.array([[float(m["parcels"]) * float(m["steps"]),
                        float(m["cells"]) * float(m["steps"])]
                       for m in usable], dtype=np.float64)
    observed = np.array([float(m["seconds"]) for m in usable], dtype=np.float64)
    solution, *_ = np.linalg.lstsq(design, observed, rcond=None)

    per_parcel, per_cell = (float(solution[0]), float(solution[1]))
    if per_parcel <= 0.0 or per_cell <= 0.0:
        # A negative coefficient is a fit, not a cost. It happens when the
        # measurements are collinear -- two runs at the same weight, say -- and
        # a model that says a bigger mesh runs faster would mislead worse than
        # the prior does.
        scale = float(np.mean([
            m["seconds"] / PRIOR_COST_MODEL.seconds(
                parcels=m["parcels"], cells=m["cells"], steps=m["steps"])
            for m in usable]))
        return CostModel(
            per_parcel_s=PRIOR_COST_MODEL.per_parcel_s * scale,
            per_cell_s=PRIOR_COST_MODEL.per_cell_s * scale,
            source=f"scaled from {len(usable)} measured run(s); the two-term "
                   f"fit was degenerate (measurements do not separate parcel "
                   f"work from cell work)",
            n_measurements=len(usable),
        )

    return CostModel(
        per_parcel_s=per_parcel,
        per_cell_s=per_cell,
        source=f"fitted to {len(usable)} measured run(s)",
        n_measurements=len(usable),
    )


def load_cost_model(path: Path) -> CostModel:
    """Read the recorded run measurements and fit a model to them.

    Args:
        path: ``results/cost-model.yaml``, written by ``AllrunCases``.

    Returns:
        The fitted model, or the prior if the file is absent or unreadable. A
        missing measurement file is the normal state before anything has run,
        not an error.
    """
    path = Path(path)
    if not path.is_file():
        return PRIOR_COST_MODEL
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError:
        return PRIOR_COST_MODEL
    return fit_cost_model(list(document.get("measurements") or []))


def record_measurement(path: Path, entry: dict) -> Path:
    """Append one measured run to ``results/cost-model.yaml``.

    Args:
        path: the file, created if absent.
        entry: at least ``case``, ``seconds``, ``parcels``, ``cells``,
            ``steps``.

    Returns:
        The path written.

    Re-running a case replaces its entry rather than adding a second one: a
    resumed run's wall clock is not the cost of the whole case, and averaging
    the two would quietly halve the estimate.
    """
    path = Path(path)
    document = {}
    if path.is_file():
        try:
            document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError:
            document = {}
    measurements = [m for m in (document.get("measurements") or [])
                    if m.get("case") != entry.get("case")]
    measurements.append(dict(entry))
    measurements.sort(key=lambda m: str(m.get("case", "")))

    model = fit_cost_model(measurements)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by cases/cai2012-health/AllrunCases -- do not edit.\n"
        "# Measured wall clock per case, and the cost model fitted to it.\n"
        "# Used by ./generate_cases.py --dry-run to cost the matrix.\n\n"
        + yaml.safe_dump({"measurements": measurements,
                          "model": model.as_dict()},
                         sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


@dataclass(frozen=True)
class CaseCost:
    """What one cell of the matrix is predicted to cost."""

    case: str
    particles_per_cell: float
    domain_transits: float
    parcels: float
    cells: int
    steps: int
    seconds: float
    writes: int
    gigabytes: float

    @property
    def hours(self) -> float:
        return self.seconds / 3600.0


def estimate_case(case: HealthCase, *, parcels: float, cells: int, steps: int,
                  writes: int, bytes_per_write: float,
                  model: CostModel) -> CaseCost:
    """Predicted wall clock and disk for one cell.

    Args:
        case: the matrix cell.
        parcels: steady-state parcel count.
        cells: mesh size.
        steps: total time steps, transient included.
        writes: time directories the run will produce.
        bytes_per_write: reconstructed size of one of them.
        model: the cost model.

    Returns:
        A :class:`CaseCost`.
    """
    return CaseCost(
        case=case.name,
        particles_per_cell=case.particles_per_cell,
        domain_transits=case.domain_transits,
        parcels=float(parcels),
        cells=int(cells),
        steps=int(steps),
        seconds=model.seconds(parcels=parcels, cells=cells, steps=steps),
        writes=int(writes),
        gigabytes=float(writes) * float(bytes_per_write) / 1024.0 ** 3,
    )


def cost_report(costs: list, model: CostModel) -> list[str]:
    """The matrix's cost, as printed by ``--dry-run``.

    Args:
        costs: one :class:`CaseCost` per cell, in generation order.
        model: the model that produced them, for the provenance line.

    Returns:
        Report lines.

    Printed **before** anything is generated, because the top weight row is
    most of the bill and finding that out from a progress bar is finding it out
    too late.
    """
    lines = [
        "Predicted cost",
        f"  model      {model.source}",
        f"             {model.per_parcel_s:.3e} s/parcel/step, "
        f"{model.per_cell_s:.3e} s/cell/step",
        "",
        f"  {'case':<16} {'ppc':>6} {'transits':>9} {'parcels':>10} "
        f"{'steps':>8} {'hours':>7} {'writes':>7} {'GB':>7}",
    ]
    for cost in costs:
        lines.append(
            f"  {cost.case:<16} {cost.particles_per_cell:>6g} "
            f"{cost.domain_transits:>9g} {cost.parcels:>10.2e} "
            f"{cost.steps:>8,} {cost.hours:>7.2f} {cost.writes:>7} "
            f"{cost.gigabytes:>7.1f}")

    total_hours = sum(c.hours for c in costs)
    total_gb = sum(c.gigabytes for c in costs)
    lines += [
        "",
        f"  TOTAL      {total_hours:.1f} h of solver, {total_gb:.0f} GB of "
        f"time directories",
    ]
    if model.n_measurements == 0:
        lines.append(
            "  These are PRIOR estimates -- nothing in this family has run yet. "
            "Run the\n"
            "  cheapest case first; AllrunCases records its wall clock and this "
            "table is\n"
            "  refitted from it.")
    return lines


# --------------------------------------------------------------------------- #
# the manifest
# --------------------------------------------------------------------------- #

def manifest_entry(case: HealthCase, study: HealthStudy, cfg, geom, exit_state,
                   plan, run, *, nozzle=None, particles=None,
                   cost: CaseCost | None = None) -> dict:
    """One manifest row from a fully derived cell.

    Everything ``cases/cai2012``'s manifest records, plus what this family
    varies and what it is predicted to cost. The mesh and time-step figures are
    included **precisely because they should be identical in every row** --
    a manifest where they are not is the evidence that the matrix stopped being
    comparable.
    """
    entry = {
        "name": case.name,
        "path": str(study.case_path(case)).replace("\\", "/"),
        "weight_level": case.weight.name,
        "sampling_level": case.sampling.name,
        "target_particles_per_cell": float(case.particles_per_cell),
        "sampling_domain_transits": float(case.domain_transits),
        "Kn": float(exit_state.knudsen),
        "n0_per_m3": float(exit_state.number_density_per_m3),
        "n_equivalent_particles": float(run.n_equivalent_particles),
        "exit_particles_per_cell_estimated": float(run.exit_particles_per_cell),
        "min_cell_size_m": float(plan.min_cell_size_m),
        "core_cell_size_m": float(plan.core_cell_size_m),
        "n_cells": int(plan.n_cells),
        "deltaT_s": float(run.delta_t_s),
        "transient_domain_transits": float(cfg.dsmc.transient_domain_transits),
        "average_start_s": float(run.average_start_s),
        "end_time_s": float(run.end_time_s),
        "write_interval_s": float(run.write_interval_s),
        "write_interval_steps": int(run.write_interval_steps),
        "n_steps": int(run.n_steps),
        "n_writes": int(run.n_writes),
        "n_sampled_writes": int(sampled_writes(run)),
    }
    if nozzle is not None:
        entry["nozzle_faces"] = int(nozzle.n_faces)
        entry["nozzle_area_error"] = float(nozzle.area_error)
    if particles is not None:
        entry["estimated_particles"] = float(f"{particles['total_particles']:.6g}")
    if cost is not None:
        entry["predicted_hours"] = float(f"{cost.hours:.4g}")
        entry["predicted_gigabytes"] = float(f"{cost.gigabytes:.4g}")
    notes = [n for n in (case.weight.note, case.sampling.note) if n]
    if notes:
        entry["note"] = " ".join(notes)
    return entry


def sampled_writes(run) -> int:
    """Time directories written **at or after** ``average_start_s``.

    The frames that carry a ``*Mean`` field, and therefore the only ones this
    family can measure or draw. The rest are written during the transient, when
    ``fieldAverage`` has not started and there is no average to look at.
    """
    if run.write_interval_s <= 0.0:
        return 0
    first = math.ceil(run.average_start_s / run.write_interval_s - 1.0e-9)
    return max(0, int(run.n_writes) - first + 1)


def write_manifest(path: Path, study: HealthStudy, entries: list,
                   model: CostModel | None = None) -> Path:
    """Write the generation manifest.

    The auditable map from case name to inputs, as in both other families, plus
    the fixed values every case shares. A reader who wants to know whether the
    matrix really was run on one mesh at one time step can find out here without
    opening nine ``case.yaml`` files.
    """
    path = Path(path)
    entries = list(entries)

    shared = {}
    if entries:
        for key in ("Kn", "n_cells", "core_cell_size_m", "deltaT_s",
                    "transient_domain_transits", "average_start_s",
                    "write_interval_s", "write_interval_steps"):
            values = {entry.get(key) for entry in entries}
            if len(values) == 1:
                shared[key] = entries[0].get(key)
            else:
                # Loud, in the generated artefact, rather than a silent column
                # of differing numbers nobody reads.
                shared[key] = (f"NOT SHARED -- {sorted(v for v in values if v is not None)}. "
                               f"The matrix is no longer comparable.")

    document = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "reference": "Cai and Wang 2012",
        "generator": "cases/cai2012-health/generate_cases.py",
        "study": "statistical health sweep: parcels per cell x averaging time",
        "base_case": study.base_case,
        "cases_dir": study.cases_dir,
        "n_generated": len(entries),
        "held_fixed": shared,
        "meta": dict(study.meta),
        "cost_model": model.as_dict() if model is not None else None,
        "cases": entries,
        "disabled": [
            {"name": case.name,
             "target_particles_per_cell": case.particles_per_cell,
             "sampling_domain_transits": case.domain_transits}
            for case in study.all_cases()
            if not (case.weight.enabled and case.sampling.enabled)
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by cases/cai2012-health/generate_cases.py -- do not edit.\n"
        "# The auditable map from case name to inputs, and what is held fixed.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False,
                         default_style=None),
        encoding="utf-8", newline="\n")
    return path


__all__ = [
    "FAMILY_KNUDSEN",
    "SWEPT_KEYS",
    "PRIOR_COST_MODEL",
    "CaseCost",
    "CostModel",
    "HealthCase",
    "HealthStudy",
    "SamplingLevel",
    "StudyError",
    "WeightLevel",
    "apply_case_parameters",
    "case_is_complete",
    "clone_cases",
    "cost_report",
    "estimate_case",
    "fit_cost_model",
    "load_cost_model",
    "load_study",
    "manifest_entry",
    "record_measurement",
    "sampled_writes",
    "write_manifest",
]
