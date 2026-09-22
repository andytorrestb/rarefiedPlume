r"""Check a generated Cai 2012 tree against the study it claims to be.

:mod:`plumetools.cai2012.checks` asks whether *one* case is a valid DSMC
calculation -- Courant number, cell size, occupancy, nozzle area. This module
asks a different question: whether the **matrix** on disk is the matrix
``study.yaml`` describes, and in particular whether the four numerical-particle
variants of a physical case differ in the particle weight **and in nothing
else**.

That distinction is the whole experiment. A convergence study in which one
variant also picked up a different mesh, a different time step or a different
density does not measure statistical resolution; it measures an accident, and it
measures it after the HPC hours have been spent. So this runs before them.

What is compared, and against what
----------------------------------
The **generated OpenFOAM dictionaries**, not the Python that produced them:

.. code-block:: text

    constant/dsmcProperties   nEquivalentParticles, species, VHS, collisions
    system/controlDict        deltaT, endTime, writeControl, writeInterval,
                              fieldAverage timeStart and fields
    system/blockMeshDict      the mesh, compared byte for byte
    case.yaml                 the physical inputs, compared key for key

Re-deriving the values in Python and comparing them with themselves would pass
whatever the case directory happened to contain. The dictionaries are what
``dsmcFoam`` reads, so they are what is read here.

The 1x case is the reference for its own group. It is not compared with a
number stored anywhere: whatever weight ``Kn0p1_np1x`` carries is *by definition*
the baseline for ``Kn0p1``, and the other three are required to be that value
divided by their multiplier. This is why the sweep cannot drift away from the
study it extends -- there is nothing to drift from.

Two tolerances, and the difference matters
------------------------------------------
Comparing two cases with each other is comparing two values written by the same
formatter, so they agree exactly or they are genuinely different:
:data:`RELATIVE_TOLERANCE` is there only to keep the arithmetic of dividing by
the multiplier honest. ``nEquivalentParticles`` is written with ``%.6e`` --
seven significant figures -- so a ratio of two written values agrees with the
exact ratio to about 1e-7, and any real mis-scaling is out by a factor of two.

Comparing a case with ``manifest.yaml`` is different: ``controlDict`` writes
``deltaT`` and ``endTime`` with ``%g``, six significant figures, while the
manifest carries the full double. ``1.3728817139518423e-06`` against
``1.37288e-06`` is a 1.2e-6 disagreement that means nothing, so that comparison
uses :data:`WRITTEN_TOLERANCE` -- wide enough for the formatting, still far too
tight to admit a value that was actually computed differently.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from plumetools.cai2012.dictionaries import AVERAGED_FIELDS

#: Relative agreement required between two values written by the same formatter.
#:
#: Set by the ``%.6e`` used for ``nEquivalentParticles``: seven significant
#: figures, so ratios are good to ~1e-7. Anything looser would still catch a
#: factor of two.
RELATIVE_TOLERANCE = 1.0e-6

#: Relative agreement required between a written value and the full-precision one.
#:
#: ``controlDict`` writes ``deltaT`` and ``endTime`` with ``%g`` -- six
#: significant figures, so up to ~5e-6 of rounding against the double in
#: ``manifest.yaml``. Wide enough for that and nothing else.
WRITTEN_TOLERANCE = 1.0e-5

#: ``case.yaml`` keys the particle axis is ALLOWED to differ on.
#:
#: Everything else in the file must match across a group, and that is the check.
#: Listing what may differ, rather than what may not, is what makes it a real
#: test: a key added to ``case.yaml`` next year is compared by default, and has
#: to be named here to be exempted.
PARTICLE_AXIS_KEYS = ("dsmc.numerical_particle_multiplier",)

#: ``case.yaml`` keys that are provenance rather than input.
IGNORED_KEYS = ("meta",)

#: Fields whose variance the study needs, and why.
REQUIRED_PRIME2MEAN = ("rhoN", "dsmcRhoN")


class VerificationError(ValueError):
    """A generated tree does not match the study that describes it."""


# --------------------------------------------------------------------------- #
# reading what OpenFOAM will read
# --------------------------------------------------------------------------- #

def foam_entry(path: Path, key: str) -> str | None:
    """The value of a top-level ``key value;`` entry, as written text.

    Deliberately the same shallow scan ``baseCase/Allrun`` does with ``sed``.
    A real dictionary parser is not needed for scalar entries and would have to
    be kept correct against OpenFOAM's syntax; this reads exactly the lines the
    generator writes, and returns ``None`` for anything it does not find rather
    than guessing.

    Args:
        path: the dictionary file.
        key: the entry name.

    Returns:
        The text between the key and the ``;``, stripped, or ``None``.
    """
    path = Path(path)
    if not path.is_file():
        return None
    pattern = re.compile(r"^\s*" + re.escape(key) + r"\s+([^;]*);")
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("//"):
            continue
        match = pattern.match(line)
        if match:
            return match.group(1).strip()
    return None


def foam_float(path: Path, key: str) -> float | None:
    """:func:`foam_entry`, as a float. ``None`` if absent or unparseable."""
    text = foam_entry(path, key)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def field_average_block(path: Path) -> dict:
    """The ``fieldAverage1`` entry of a ``controlDict``, as a plain dict.

    Returns:
        ``{"timeStart": float | None, "fields": {name: {"mean": bool,
        "prime2Mean": bool}}}``. An absent function object gives an empty
        ``fields``, which every caller treats as a failure rather than as a
        default.
    """
    path = Path(path)
    result = {"timeStart": None, "fields": {}}
    if not path.is_file():
        return result

    text = path.read_text(encoding="utf-8")
    start = text.find("fieldAverage1")
    if start < 0:
        return result
    block = text[start:]

    match = re.search(r"^\s*timeStart\s+([^;]*);", block, re.MULTILINE)
    if match:
        try:
            result["timeStart"] = float(match.group(1).strip())
        except ValueError:
            pass

    # Scan the `fields ( ... )` list, not the whole function object. Scanning
    # from `fieldAverage1` would match ITS opening brace first and, because the
    # body match is lazy, swallow everything up to the first field's closing
    # brace -- silently losing the first field, which here is rhoN.
    fields = re.search(r"^\s*fields\s*$\s*\((.*)\)\s*;?", block,
                       re.MULTILINE | re.DOTALL)
    if not fields:
        return result

    # Each field is `name { mean on; prime2Mean off; base time; }`. Matching the
    # name to the braces that follow it is enough here because the generator
    # writes exactly that shape and nothing nests inside it.
    for entry in re.finditer(r"^\s*(\w+)\s*$\s*\{(.*?)\}", fields.group(1),
                             re.MULTILINE | re.DOTALL):
        name, body = entry.group(1), entry.group(2)
        if "mean" not in body:
            continue
        result["fields"][name] = {
            "mean": bool(re.search(r"\bmean\s+on\s*;", body)),
            "prime2Mean": bool(re.search(r"\bprime2Mean\s+on\s*;", body)),
        }
    return result


def _flatten(data, prefix: str = "") -> dict:
    """A nested mapping as ``{"dsmc.delta_t_s": value}``."""
    flat = {}
    if isinstance(data, dict):
        for key, value in data.items():
            flat.update(_flatten(value, f"{prefix}{key}."))
        return flat
    flat[prefix.rstrip(".")] = data
    return flat


def case_inputs(case_dir: Path) -> dict:
    """A generated ``case.yaml``, flattened, with provenance removed.

    What is left is every physical and numerical input the case will be built
    from -- which is exactly what has to be identical across a particle group
    apart from :data:`PARTICLE_AXIS_KEYS`.
    """
    path = Path(case_dir) / "case.yaml"
    if not path.is_file():
        raise VerificationError(f"no case.yaml at {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    for key in IGNORED_KEYS:
        data.pop(key, None)
    return _flatten(data)


# --------------------------------------------------------------------------- #
# one case, read off disk
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class CaseFacts:
    """What a generated case directory actually says it will run.

    Every field is read from the case, never re-derived, so a mismatch between
    this and the manifest is a real disagreement rather than two evaluations of
    the same expression.
    """

    name: str
    path: Path
    cai_case: str
    particle_level: str
    multiplier: float
    n_equivalent_particles: float | None
    delta_t_s: float | None
    start_time: str | None
    end_time_s: float | None
    write_control: str | None
    write_interval: float | None
    average_start_s: float | None
    averaged_fields: dict
    block_mesh_digest: str | None
    inputs: dict

    @property
    def has_dictionaries(self) -> bool:
        """Whether ``./Allmesh`` has written this case's dictionaries yet."""
        return (self.n_equivalent_particles is not None
                and self.delta_t_s is not None)


def read_case(case_dir: Path, *, name: str | None = None) -> CaseFacts:
    """Read one generated case directory.

    Args:
        case_dir: the case.
        name: its name; defaults to the directory name.

    Returns:
        A :class:`CaseFacts`.

    Raises:
        VerificationError: if the directory has no ``case.yaml``, i.e. is not a
            generated case at all.
    """
    case_dir = Path(case_dir)
    inputs = case_inputs(case_dir)
    meta = yaml.safe_load(
        (case_dir / "case.yaml").read_text(encoding="utf-8")).get("meta") or {}

    dsmc_properties = case_dir / "constant" / "dsmcProperties"
    control_dict = case_dir / "system" / "controlDict"
    block_mesh = case_dir / "system" / "blockMeshDict"

    digest = None
    if block_mesh.is_file():
        import hashlib
        digest = hashlib.sha256(block_mesh.read_bytes()).hexdigest()[:16]

    averaging = field_average_block(control_dict)

    return CaseFacts(
        name=name or case_dir.name,
        path=case_dir,
        cai_case=str(meta.get("cai_case", "")),
        particle_level=str(meta.get("particle_level", "")),
        multiplier=float(inputs.get("dsmc.numerical_particle_multiplier", 1.0)),
        n_equivalent_particles=foam_float(dsmc_properties,
                                          "nEquivalentParticles"),
        delta_t_s=foam_float(control_dict, "deltaT"),
        start_time=foam_entry(control_dict, "startTime"),
        end_time_s=foam_float(control_dict, "endTime"),
        write_control=foam_entry(control_dict, "writeControl"),
        write_interval=foam_float(control_dict, "writeInterval"),
        average_start_s=averaging["timeStart"],
        averaged_fields=averaging["fields"],
        block_mesh_digest=digest,
        inputs=inputs,
    )


# --------------------------------------------------------------------------- #
# the report
# --------------------------------------------------------------------------- #

@dataclass
class VerificationReport:
    """Findings, in the order they were made."""

    failures: list = field(default_factory=list)
    passes: list = field(default_factory=list)

    def fail(self, where: str, message: str) -> None:
        self.failures.append(f"{where}: {message}")

    def ok(self, where: str, message: str) -> None:
        self.passes.append(f"{where}: {message}")

    @property
    def passed(self) -> bool:
        return not self.failures

    def lines(self, *, verbose: bool = False) -> list[str]:
        out = []
        if verbose:
            out += [f"  ok    {p}" for p in self.passes]
        out += [f"  FAIL  {f}" for f in self.failures]
        return out


def _close(a: float, b: float, tolerance: float = RELATIVE_TOLERANCE) -> bool:
    """Relative comparison that does not divide by zero."""
    if a is None or b is None:
        return False
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=0.0)


def verify_group(cases: list, report: VerificationReport, *,
                 tolerance: float = RELATIVE_TOLERANCE) -> None:
    """Check one physical case's particle variants against each other.

    ``cases`` are the :class:`CaseFacts` for one Cai case at every level. The 1x
    member is the reference; the rest must be it, divided.
    """
    if not cases:
        return
    group = cases[0].cai_case or cases[0].name
    baselines = [c for c in cases if c.multiplier == 1.0]
    if not baselines:
        report.fail(group, "no 1x variant, so there is no baseline to compare "
                           "the others against")
        return
    if len(baselines) > 1:
        report.fail(group, f"{len(baselines)} variants claim to be 1x: "
                           f"{[c.name for c in baselines]}")
        return
    reference = baselines[0]

    if reference.n_equivalent_particles is None:
        report.fail(reference.name,
                    "no nEquivalentParticles in constant/dsmcProperties. The "
                    "dictionaries are generated by ./Allmesh; run "
                    "'./AllmeshCases --dict-only' first.")
        return
    nep0 = reference.n_equivalent_particles

    # --- the axis itself -----------------------------------------------------
    for case in cases:
        expected = nep0 / case.multiplier
        actual = case.n_equivalent_particles
        if actual is None:
            report.fail(case.name, "no nEquivalentParticles in "
                                   "constant/dsmcProperties")
        elif not _close(actual, expected, tolerance):
            report.fail(
                case.name,
                f"nEquivalentParticles is {actual:.6e}, expected "
                f"{expected:.6e} = {nep0:.6e} / {case.multiplier:g}. The "
                f"parcel population would be "
                f"{nep0 / actual:.4g}x the baseline, not {case.multiplier:g}x.")
        else:
            report.ok(case.name,
                      f"nEquivalentParticles {actual:.6e} = {nep0:.6e} / "
                      f"{case.multiplier:g}  ({case.multiplier:g}x parcels)")

    # --- everything that must NOT move --------------------------------------
    for case in cases:
        if case is reference:
            continue

        if not _close(case.delta_t_s, reference.delta_t_s, tolerance):
            report.fail(
                case.name,
                f"deltaT is {case.delta_t_s!r} against the 1x case's "
                f"{reference.delta_t_s!r}. This sweep varies the particle "
                f"count ONLY; a different time step would confound statistical "
                f"resolution with discretisation.")
        else:
            report.ok(case.name, f"deltaT {case.delta_t_s:g}, as 1x")

        if case.block_mesh_digest != reference.block_mesh_digest:
            report.fail(
                case.name,
                f"system/blockMeshDict differs from the 1x case "
                f"({case.block_mesh_digest} against "
                f"{reference.block_mesh_digest}); the variants are not on the "
                f"same mesh")
        elif case.block_mesh_digest is not None:
            report.ok(case.name,
                      f"same mesh as 1x (blockMeshDict {case.block_mesh_digest})")

        # The physical inputs, key by key. Only the particle multiplier may
        # differ -- see PARTICLE_AXIS_KEYS for why this is stated as an
        # exemption list rather than a comparison list.
        differing = sorted(
            key for key in set(case.inputs) | set(reference.inputs)
            if case.inputs.get(key) != reference.inputs.get(key))
        unexpected = [k for k in differing if k not in PARTICLE_AXIS_KEYS]
        if unexpected:
            detail = ", ".join(
                f"{k}: {reference.inputs.get(k)!r} -> {case.inputs.get(k)!r}"
                for k in unexpected[:6])
            report.fail(
                case.name,
                f"case.yaml differs from the 1x case in {len(unexpected)} key(s) "
                f"beyond the particle multiplier: {detail}")
        else:
            report.ok(case.name,
                      "case.yaml identical to 1x apart from "
                      "dsmc.numerical_particle_multiplier (same mesh, gas, "
                      "geometry, boundary conditions, collision model, "
                      "density and deltaT)")

        if not _close(case.average_start_s, reference.average_start_s, tolerance):
            report.fail(
                case.name,
                f"fieldAverage timeStart is {case.average_start_s!r} against "
                f"the 1x case's {reference.average_start_s!r}; the variants "
                f"would be averaged over different transients")

        if not _close(case.end_time_s, reference.end_time_s, tolerance):
            report.fail(
                case.name,
                f"endTime is {case.end_time_s!r} against the 1x case's "
                f"{reference.end_time_s!r}; the variants would be sampled for "
                f"different lengths of time")

        if case.write_control != reference.write_control:
            report.fail(
                case.name,
                f"writeControl is {case.write_control!r} against the 1x case's "
                f"{reference.write_control!r}")


def verify_schedule(case: CaseFacts, expected: dict,
                    report: VerificationReport, *,
                    tolerance: float = WRITTEN_TOLERANCE) -> None:
    """Check one case's run length, write schedule and averaging.

    Args:
        case: the case, read off disk.
        expected: the manifest row for it.
        report: accumulates findings.
        tolerance: relative agreement required. Defaults to
            :data:`WRITTEN_TOLERANCE`, because this compares ``%g``-formatted
            dictionary entries with the full-precision manifest.

    The manifest is the record of what generation intended; the dictionaries are
    what the solver will do. Disagreement between them means a case was edited,
    or half-regenerated, after the manifest was written.
    """
    where = case.name

    for key, actual, label in (
            ("deltaT_s", case.delta_t_s, "deltaT"),
            ("end_time_s", case.end_time_s, "endTime"),
            ("average_start_s", case.average_start_s,
             "fieldAverage timeStart"),
            ("write_interval_steps", case.write_interval, "writeInterval"),
            ("n_equivalent_particles", case.n_equivalent_particles,
             "nEquivalentParticles")):
        want = expected.get(key)
        if want is None or actual is None:
            report.fail(where, f"{label} is missing from the case or from the "
                               f"manifest")
        elif not _close(float(actual), float(want), tolerance):
            report.fail(where, f"{label} is {actual!r} but manifest.yaml says "
                               f"{want!r}; the case and the manifest disagree")

    want_control = expected.get("write_control")
    if want_control and case.write_control != want_control:
        report.fail(where, f"writeControl is {case.write_control!r} but "
                           f"manifest.yaml says {want_control!r}")

    if case.start_time is not None and float(case.start_time) != 0.0:
        report.fail(where, f"startTime is {case.start_time!r}, expected 0")

    # --- the requested schedule ---------------------------------------------
    baseline_end = expected.get("baseline_end_time_s")
    run_multiplier = expected.get("run_time_multiplier")
    if baseline_end and run_multiplier and case.end_time_s:
        achieved = case.end_time_s / float(baseline_end)
        # The end time is snapped UP to a whole number of write intervals, so it
        # lands at or just above the request -- never below it, and never by
        # more than one interval.
        interval = (case.write_interval or 0) * (case.delta_t_s or 0)
        slack = 1.0 + (interval / float(baseline_end) if baseline_end else 0.0)
        if not float(run_multiplier) <= achieved <= float(run_multiplier) * slack:
            report.fail(
                where,
                f"endTime is {achieved:.4f}x the baseline {baseline_end:g} s, "
                f"but run_time_multiplier is {run_multiplier:g}")
        else:
            report.ok(where, f"endTime {case.end_time_s:g} s = "
                             f"{achieved:.4f}x the baseline (requested "
                             f"{float(run_multiplier):g}x)")

    baseline_interval = expected.get("baseline_write_interval_steps")
    output_multiplier = expected.get("output_frequency_multiplier")
    if baseline_interval and output_multiplier and case.write_interval:
        achieved = float(baseline_interval) / case.write_interval
        # writeControl timeStep takes an integer, so the request is met to
        # within that rounding. Half a step of the baseline interval is the
        # widest the rounding can be.
        allowed = float(output_multiplier) * (
            1.0 + 0.5 / max(1.0, case.write_interval))
        if not abs(achieved - float(output_multiplier)) <= abs(
                allowed - float(output_multiplier)) + 1.0e-9:
            report.fail(
                where,
                f"writeInterval {case.write_interval:g} steps is "
                f"{achieved:.4f}x more frequent than the baseline "
                f"{baseline_interval:g}, but output_frequency_multiplier is "
                f"{output_multiplier:g}")
        else:
            report.ok(where, f"writeInterval {case.write_interval:g} steps = "
                             f"{achieved:.4f}x more frequent than the baseline "
                             f"{baseline_interval:g} (requested "
                             f"{float(output_multiplier):g}x)")

    # --- averaging -----------------------------------------------------------
    if not case.averaged_fields:
        report.fail(where, "system/controlDict has no fieldAverage1 function "
                           "object; nothing would be averaged and every written "
                           "field would be one timestep of shot noise")
        return

    for name in AVERAGED_FIELDS:
        if name not in case.averaged_fields:
            report.fail(where, f"fieldAverage does not average {name!r}")
        elif not case.averaged_fields[name]["mean"]:
            report.fail(where, f"fieldAverage has mean off for {name!r}")

    for name in REQUIRED_PRIME2MEAN:
        entry = case.averaged_fields.get(name)
        if entry is None or not entry["prime2Mean"]:
            report.fail(
                where,
                f"fieldAverage has prime2Mean off for {name!r}, so no "
                f"{name}Prime2Mean is written and the statistical scatter this "
                f"study exists to measure -- CV = sqrt({name}Prime2Mean) / "
                f"{name}Mean -- cannot be evaluated afterwards")
        else:
            report.ok(where, f"{name}Mean and {name}Prime2Mean will be written")


def verify_tree(root: Path, manifest: dict, *,
                tolerance: float = RELATIVE_TOLERANCE) -> VerificationReport:
    """Check every case the manifest lists.

    Args:
        root: the study directory -- the one holding ``manifest.yaml``.
        manifest: the loaded manifest.
        tolerance: relative agreement required of values that should be exact.

    Returns:
        A :class:`VerificationReport`. Nothing is raised: the caller decides
        what a failure means, and a report naming every bad case is more useful
        than an exception naming the first.
    """
    root = Path(root)
    report = VerificationReport()

    entries = manifest.get("cases") or []
    if not entries:
        report.fail("manifest.yaml", "lists no cases")
        return report

    names = [e.get("name") for e in entries]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        report.fail("manifest.yaml",
                    f"duplicate case name(s) {duplicates}; two entries would "
                    f"share a directory and the second would overwrite the "
                    f"first")

    paths = [e.get("path") for e in entries]
    duplicate_paths = sorted({p for p in paths if paths.count(p) > 1})
    if duplicate_paths:
        report.fail("manifest.yaml", f"duplicate case path(s) {duplicate_paths}")

    groups: dict = {}
    for entry in entries:
        case_dir = root / entry["path"]
        if not case_dir.is_dir():
            report.fail(entry["name"], f"{case_dir} does not exist")
            continue
        try:
            facts = read_case(case_dir, name=entry["name"])
        except VerificationError as exc:
            report.fail(entry["name"], str(exc))
            continue

        if not facts.has_dictionaries:
            report.fail(
                entry["name"],
                "has no generated OpenFOAM dictionaries. They are written by "
                "./Allmesh from case.yaml; run './AllmeshCases --dict-only' "
                "(no OpenFOAM needed) or './AllmeshCases' first.")
            continue

        # The name and the file must agree about what this case is. If they do
        # not, everything downstream is comparing the wrong things.
        claimed = entry.get("cai_case")
        if claimed and facts.cai_case and facts.cai_case != claimed:
            report.fail(entry["name"],
                        f"case.yaml says it belongs to {facts.cai_case!r}, the "
                        f"manifest says {claimed!r}")
        claimed_multiplier = entry.get("numerical_particle_multiplier")
        if claimed_multiplier is not None and not _close(
                facts.multiplier, float(claimed_multiplier), tolerance):
            report.fail(entry["name"],
                        f"case.yaml has numerical_particle_multiplier "
                        f"{facts.multiplier:g}, the manifest says "
                        f"{float(claimed_multiplier):g}")

        verify_schedule(facts, entry, report, tolerance=WRITTEN_TOLERANCE)
        groups.setdefault(facts.cai_case or entry["name"], []).append(facts)

    for group in groups.values():
        verify_group(sorted(group, key=lambda c: c.multiplier), report,
                     tolerance=tolerance)

    return report


def load_manifest(path: Path) -> dict:
    """Read ``manifest.yaml``.

    Raises:
        VerificationError: if it is missing, so the caller says "generate the
            cases first" rather than reporting an empty study as a clean one.
    """
    path = Path(path)
    if path.is_dir():
        path = path / "manifest.yaml"
    if not path.is_file():
        raise VerificationError(
            f"no manifest.yaml at {path}. Generate the cases first:\n"
            f"    ./generate_cases.py")
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
