r"""Matrix profiles, compatible overrides and inherited numerical settings.

With no declaration, the legacy profile retains exactly Kn 1e-4, 1e-3, 1, 10
and refuses Cai's paper Kn even when disabled. The explicit
``meta.knudsen_suite: full_decades`` profile requires all seven decade cases
100, 10, 1, .1, .01, .001, .0001. Both use one baseline particle level.

The profile is carried from study to manifest and generated case metadata.
Generation, template checks, manifest runners, HPC guards and direct-exit
post-processing share it. Unknown profiles are refused. Matrix membership
does not waive numerical feasibility or resolve the known very-low-Kn source
surface limitation.

The template inherits the health study's ppc040/s4p5 settings, with the declared
overrides below. Per-Kn mesh, time-step and parcel-target changes are allowed;
gas, geometry, speed ratio, startup and averaging remain fixed.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import sys
from pathlib import Path

import yaml

from plumetools.cai2012.study import StudyError
from plumetools.cai2012ext import (
    EXCLUDED_KNUDSEN,
    EXTENSION_KNUDSEN,
    FULL_SUITE_KNUDSEN,
    PARTICLE_MULTIPLIERS,
    REFERENCE_KNUDSEN,
)


class ExtensionError(StudyError):
    """The extension's matrix, template or manifest breaks its policy."""


#: ``case.yaml`` sections, in the order :class:`~plumetools.cai2012.config.CaiCaseConfig` declares them.
SECTIONS = ("nozzle", "exit", "gas", "geometry", "mesh", "dsmc", "resolution",
            "checks", "output", "post")

#: The keys in which the extension template may differ from
#: ``cases/cai2012-health/baseCase/case.yaml``, and why.
#:
#: Every other key is inherited, and :func:`inheritance_report` fails on any
#: that is not. ``meta`` is provenance and is never compared.
DELIBERATE_OVERRIDES = {
    "exit.knudsen": (
        "written per case by generate_cases.py. The health family pins 100, "
        "which is one of the excluded Kn; the template carries 10.0, an "
        "extension value, so that even ./Allmesh run in baseCase/ cannot "
        "produce an excluded case."),
    "resolution.target_particles_per_cell": (
        "40 -- the health study's highest-resolution weight row, ppc040 "
        "(cases/cai2012-health/study.yaml). The health TEMPLATE says 20 "
        "because its generator writes the row value per case. This is the "
        "baseline, and the only level the extension runs: 40 parcels in the "
        "exit cell."),
    "dsmc.sampling_domain_transits": (
        "4.5 -- the health study's longest averaging column, s4p5 "
        "(cases/cai2012-health/study.yaml). The template says 1.5 for the same "
        "reason as above. Averaging begins AFTER the startup below."),
    "dsmc.transient_basis": (
        "'cai': Cai's 10^4 reference collision times, the reference being the "
        "Kn = 0.01 exit state (resolution.reference_knudsen). NOT the health "
        "family's 2-domain-transit startup, which that family keeps to match "
        "its historical comparison and documents as limiting the far field "
        "(docs/cai2012-case.md section 5.1)."),
    "dsmc.n_subdomains": (
        "192 -- one HPC node of 192 cores, on which each case runs as its own "
        "PBS job (runCase.sh). ./Allmesh writes it into "
        "system/decomposeParDict and ./Allrun reads it from there, so no job "
        "passes a rank count. The health template's 12 is its development "
        "machine's physical cores. A machine setting, not a physical one."),
}

#: Where the two swept health values were read from, for the record.
HEALTH_STUDY_ROWS = {
    "resolution.target_particles_per_cell": ("weights", "ppc040", 40.0),
    "dsmc.sampling_domain_transits": ("samplings", "s4p5", 4.5),
}

#: ``case.yaml`` keys a Knudsen case in ``study.yaml`` may override.
#:
#: An ALLOW-list, because the invariants the extension promises -- the same D,
#: S0, T0 and gas everywhere, the same startup policy, the same averaging, the
#: same estimators -- are exactly the keys a per-Kn override would break. The
#: mesh budget and the time step are allowed because resolving a dense case is
#: precisely what a per-Kn mesh override is for. The exit-cell parcel target is
#: allowed because a mesh override changes the exit cell: on a uniform mesh the
#: weight has to be set by the far field (20 parcels at n/n0 = 1e-2) rather than
#: the 40 exit-cell parcels the graded mesh's growing outer cells make enough.
ALLOWED_OVERRIDE_SECTIONS = ("mesh",)
ALLOWED_OVERRIDE_KEYS = ("dsmc.courant_target", "dsmc.delta_t_s",
                         "dsmc.n_subdomains",
                         "resolution.target_particles_per_cell")

#: Relative tolerance for "the same Knudsen number". Values are typed as
#: decimals and pass through YAML, so exact float equality would be a trap.
KNUDSEN_RTOL = 1.0e-9
LEGACY_SUITE = "legacy_extension"
FULL_SUITE = "full_decades"


# --------------------------------------------------------------------------- #
# Knudsen numbers
# --------------------------------------------------------------------------- #

def same_knudsen(a, b) -> bool:
    """Whether two Knudsen numbers are the same number."""
    return math.isclose(float(a), float(b), rel_tol=KNUDSEN_RTOL, abs_tol=0.0)


def is_excluded(knudsen) -> bool:
    """Whether Kn is on the legacy paper-Kn exclusion list."""
    return any(same_knudsen(knudsen, k) for k in EXCLUDED_KNUDSEN)


def is_extension_knudsen(knudsen) -> bool:
    """Whether Kn belongs to the legacy four-case extension."""
    return any(same_knudsen(knudsen, k) for k in EXTENSION_KNUDSEN)


def suite_name(source=None) -> str:
    """Resolve the explicit matrix profile from a study, manifest or case.

    A missing declaration preserves the historical four-case policy. A string
    may be supplied by internal callers, but arbitrary profiles never default
    to the permissive seven-decade matrix.
    """
    if source is None:
        value = LEGACY_SUITE
    elif isinstance(source, str):
        value = source
    elif isinstance(source, dict):
        value = (source.get("meta") or source).get("knudsen_suite", LEGACY_SUITE)
    else:
        value = (getattr(source, "meta", None) or {}).get("knudsen_suite", LEGACY_SUITE)
    if value not in (LEGACY_SUITE, FULL_SUITE):
        raise ExtensionError(f"Unknown knudsen_suite {value!r}; expected {LEGACY_SUITE!r} or {FULL_SUITE!r}")
    return value


def suite_knudsen(source=None) -> tuple:
    """The allowed Kn values for the selected explicit matrix profile."""
    return FULL_SUITE_KNUDSEN if suite_name(source) == FULL_SUITE else EXTENSION_KNUDSEN


def is_suite_knudsen(knudsen, source=None) -> bool:
    return any(same_knudsen(knudsen, k) for k in suite_knudsen(source))


def source_limitations(cfg) -> list[str]:
    """Known direct-exit source issue, independent of mesh/time-step checks.

    The very-low-Kn inflow-surface correction is deliberately deferred by the
    user. This warning records that choice; it does not validate the source.
    """
    if suite_name(cfg) == FULL_SUITE and float(cfg.exit.knudsen) <= 0.001:
        return ["known_low_kn_inflow_surface: direct nozzle-exit inflow is retained; "
                "the correct inflow surface for very low Kn_e is deferred to future work"]
    return []


def direct_exit_problem(cfg, where="case.yaml") -> str | None:
    """Common per-case gate for direct-exit run and post-processing paths."""
    try:
        profile = suite_name(cfg)
    except ExtensionError as exc:
        return str(exc)
    if is_excluded(cfg.exit.knudsen) and profile != FULL_SUITE:
        return _excluded_message(cfg.exit.knudsen, where)
    if profile == FULL_SUITE and not is_suite_knudsen(cfg.exit.knudsen, cfg):
        return f"{where}: Kn = {cfg.exit.knudsen:g} is outside the explicit full-decade suite"
    return None


def _excluded_message(knudsen, where: str) -> str:
    return (f"{where}: Kn = {float(knudsen):g} is one of Cai's own Knudsen "
            f"numbers {list(EXCLUDED_KNUDSEN)}. Those are the conditions his "
            f"paper already compares, and the legacy extension profile never generates, "
            f"meshes, runs or post-processes them -- whether or not results "
            f"for them exist elsewhere. They belong to cases/cai2012.")


# --------------------------------------------------------------------------- #
# the study
# --------------------------------------------------------------------------- #

def check_study(study) -> None:
    """Refuse a study that is not exactly its declared matrix profile.

    Args:
        study: a :class:`~plumetools.cai2012.study.StudyConfig`, loaded with
            the parent family's own loader.

    Raises:
        ExtensionError: listing **every** problem, not the first.
    """
    problems = []

    wanted = suite_knudsen(study)
    full_suite = suite_name(study) == FULL_SUITE

    for case in study.all_cases():
        if is_excluded(case.knudsen) and not full_suite:
            state = "enabled" if case.enabled else "disabled"
            problems.append(_excluded_message(
                case.knudsen, f"kn_cases {case.name!r} ({state})"))

    enabled = [c.knudsen for c in study.enabled_cases()]
    missing = [k for k in wanted
               if not any(same_knudsen(k, e) for e in enabled)]
    extra = [e for e in enabled if not is_suite_knudsen(e, study)
             and (full_suite or not is_excluded(e))]
    if len(enabled) != len(set(float(k) for k in enabled)):
        problems.append("kn_cases enables duplicate Kn values")
    if missing:
        problems.append(
            f"kn_cases does not enable Kn {missing}. The extension is exactly "
            f"{list(wanted)}; a matrix with one missing is not the "
            f"study the manifest would claim it is.")
    if extra:
        problems.append(
            f"kn_cases enables Kn {extra}, which is not part of the extension "
            f"{list(wanted)}")

    multipliers = [float(p.multiplier) for p in study.enabled_particle_levels()]
    if sorted(multipliers) != sorted(PARTICLE_MULTIPLIERS):
        problems.append(
            f"particle_levels enables multipliers {multipliers}; the extension "
            f"runs the baseline only, {list(PARTICLE_MULTIPLIERS)}. There is no "
            f"parcel-population sweep in this study.")

    # The parent study runs 3x longer and writes 4x more often through these.
    # Inherited by accident they would multiply the 4.5-transit averaging window
    # and the health write cadence, so they are pinned at one here.
    for key in ("run_time_multiplier", "output_frequency_multiplier"):
        value = float(getattr(study, key))
        if value != 1.0:
            problems.append(
                f"{key} is {value:g}. The extension averages for exactly "
                f"dsmc.sampling_domain_transits after the startup and writes at "
                f"the health cadence; a study-wide multiplier would silently "
                f"change both. It must be 1.0.")

    for case in study.all_cases():
        for section, values in (case.overrides or {}).items():
            if section in ALLOWED_OVERRIDE_SECTIONS:
                continue
            keys = values if isinstance(values, dict) else {None: values}
            for key in keys:
                dotted = f"{section}.{key}" if key is not None else section
                if dotted not in ALLOWED_OVERRIDE_KEYS:
                    problems.append(
                        f"kn_cases {case.name!r} overrides {dotted!r}. Only "
                        f"{list(ALLOWED_OVERRIDE_SECTIONS)} sections and "
                        f"{list(ALLOWED_OVERRIDE_KEYS)} may be set per Kn: "
                        f"everything else is held fixed across the whole "
                        f"extension (D, S0, T0, gas, startup, averaging, "
                        f"estimators).")

    if problems:
        raise ExtensionError(
            "study.yaml is not the extension matrix:\n  - "
            + "\n  - ".join(problems))


def expected_case_count(source=None) -> int:
    """Case count for the matrix profile (legacy four, explicit full seven)."""
    return len(suite_knudsen(source)) * len(PARTICLE_MULTIPLIERS)


# --------------------------------------------------------------------------- #
# the template
# --------------------------------------------------------------------------- #

def check_template(cfg, where: str = "baseCase/case.yaml", *, suite=None) -> None:
    """Refuse a template that would break the extension's fixed policies.

    Args:
        cfg: the loaded template, a
            :class:`~plumetools.cai2012.config.CaiCaseConfig`.
        where: for messages.

    Raises:
        ExtensionError: listing every problem.
    """
    problems = []
    profile = cfg if suite is None else suite
    profile_name = suite_name(profile)
    if is_excluded(cfg.exit.knudsen) and profile_name != FULL_SUITE:
        problems.append(_excluded_message(cfg.exit.knudsen, "exit.knudsen"))
    elif profile_name == FULL_SUITE and not is_suite_knudsen(cfg.exit.knudsen, profile):
        problems.append(f"exit.knudsen = {cfg.exit.knudsen:g} is outside the explicit full-decade suite")
    if not same_knudsen(cfg.resolution.reference_knudsen, REFERENCE_KNUDSEN):
        problems.append(
            f"resolution.reference_knudsen is "
            f"{cfg.resolution.reference_knudsen:g}. Cai's reference collision "
            f"time -- what the 10^4 t0 startup is counted in -- is the "
            f"Kn = {REFERENCE_KNUDSEN:g} exit state, for every case.")
    if cfg.dsmc.transient_basis != "cai":
        problems.append(
            f"dsmc.transient_basis is {cfg.dsmc.transient_basis!r}; the "
            f"extension starts sampling after Cai's 10^4 reference collision "
            f"times ('cai')")
    for key in ("average_start_s", "end_time_s", "n_equivalent_particles"):
        if getattr(cfg.dsmc, key) is not None:
            problems.append(
                f"dsmc.{key} is pinned at {getattr(cfg.dsmc, key)!r}. It must "
                f"be derived (null): a pinned value would not follow the case's "
                f"own density or the fixed Cai startup.")
    for key in ("numerical_particle_multiplier", "run_time_multiplier",
                "output_frequency_multiplier"):
        if float(getattr(cfg.dsmc, key)) != 1.0:
            problems.append(
                f"dsmc.{key} is {getattr(cfg.dsmc, key)!r} in the template. "
                f"The particle level is applied once, by generation, and the "
                f"run length and write cadence are not scaled at all.")
    if not cfg.dsmc.collisions_enabled:
        problems.append("dsmc.collisions_enabled is false; collisions stay on at "
                        "every Kn")
    if problems:
        raise ExtensionError(f"{where} breaks the extension policy:\n  - "
                             + "\n  - ".join(problems))


def flatten_config(cfg) -> dict:
    """``{"section.key": value}`` for every input key of a loaded config.

    From the loaded dataclasses rather than the raw YAML, so a key one file
    states and the other leaves to its default still compares like for like.
    """
    flat = {}
    for section in SECTIONS:
        values = getattr(cfg, section)
        for item in dataclasses.fields(values):
            value = getattr(values, item.name)
            if isinstance(value, tuple):
                value = list(value)
            flat[f"{section}.{item.name}"] = value
    return flat


def inheritance_report(extension_cfg, health_cfg) -> dict:
    """Compare the extension template with the health template, key by key.

    Args:
        extension_cfg: ``cases/cai2012-kn-extension/baseCase/case.yaml``, loaded.
        health_cfg: ``cases/cai2012-health/baseCase/case.yaml``, loaded.

    Returns:
        ``{"inherited": {key: value}, "overridden": {key: {...}},
        "undeclared": {key: {...}}}``. ``undeclared`` must be empty: it lists
        keys that differ without being in :data:`DELIBERATE_OVERRIDES`.
    """
    ours = flatten_config(extension_cfg)
    theirs = flatten_config(health_cfg)
    inherited, overridden, undeclared = {}, {}, {}
    for key in sorted(set(ours) | set(theirs)):
        mine, health = ours.get(key), theirs.get(key)
        if key in DELIBERATE_OVERRIDES:
            overridden[key] = {"extension": mine, "health": health,
                               "reason": DELIBERATE_OVERRIDES[key]}
        elif mine == health:
            inherited[key] = mine
        else:
            undeclared[key] = {"extension": mine, "health": health}
    return {"inherited": inherited, "overridden": overridden,
            "undeclared": undeclared}


def check_health_study_rows(health_study_path: Path) -> list[str]:
    """Problems with the claim that 40 and 4.5 are health's ppc040 / s4p5.

    Reads the health ``study.yaml`` directly. Returns an empty list when both
    rows exist with the values :data:`HEALTH_STUDY_ROWS` records.
    """
    raw = yaml.safe_load(Path(health_study_path).read_text(encoding="utf-8")) or {}
    problems = []
    value_keys = {"weights": "target_particles_per_cell",
                  "samplings": "sampling_domain_transits"}
    for key, (axis, name, value) in HEALTH_STUDY_ROWS.items():
        rows = {row.get("name"): row for row in raw.get(axis) or []}
        if name not in rows:
            problems.append(f"{health_study_path}: no {axis} row {name!r} (the "
                            f"source of {key})")
            continue
        found = float(rows[name][value_keys[axis]])
        if found != value:
            problems.append(f"{health_study_path}: {axis}[{name}] is {found:g}, "
                            f"but the extension inherits {key} = {value:g} from it")
    return problems


# --------------------------------------------------------------------------- #
# the manifest -- the execution paths
# --------------------------------------------------------------------------- #

def check_manifest(manifest: dict) -> list[str]:
    """Policy problems with a manifest that is about to be executed.

    Returns a list of messages; empty means every listed case is one of the
    extension's. An excluded Kn anywhere in the manifest is reported whatever else
    is right, because the manifest is what every runner script reads.
    """
    problems = []
    try:
        wanted = suite_knudsen(manifest)
        full_suite = suite_name(manifest) == FULL_SUITE
    except ExtensionError as exc:
        return [str(exc)]
    entries = manifest.get("cases") or []
    for entry in entries:
        kn = entry.get("Kn")
        name = entry.get("name", "?")
        if kn is None:
            problems.append(f"{name}: no Kn in the manifest row")
        elif is_excluded(kn) and not full_suite:
            problems.append(_excluded_message(kn, f"manifest case {name!r}"))
        elif not any(same_knudsen(kn, value) for value in wanted):
            problems.append(f"manifest case {name!r} is at Kn = {float(kn):g}, "
                            f"not one of {list(wanted)}")
    return problems


def production_status(entry: dict) -> str:
    """The feasibility verdict recorded in a manifest row."""
    return str((entry.get("feasibility") or {}).get("status", "unassessed"))


def schedulable(manifest: dict, *, allow_unsuitable: bool = False):
    """Split the manifest into production candidates and the rest.

    Returns:
        ``(runnable, refused)`` -- ``runnable`` the candidates, or every case
        with ``allow_unsuitable``; ``refused`` lists ``(entry, reasons)``.

    Raises:
        ExtensionError: if the manifest contains an excluded Kn at all. That is
            not a case to skip; it means the manifest is not this study's.
    """
    problems = check_manifest(manifest)
    if problems:
        raise ExtensionError("manifest.yaml breaks the extension policy:\n  - "
                             + "\n  - ".join(problems))
    runnable, refused = [], []
    for entry in manifest.get("cases") or []:
        status = production_status(entry)
        if status == "production-candidate" or allow_unsuitable:
            runnable.append(entry)
        else:
            reasons = (entry.get("feasibility") or {}).get("unresolved") or [
                f"feasibility status is {status!r}"]
            refused.append((entry, reasons))
    return runnable, refused


def main(argv=None) -> int:
    """List the manifest's case paths for a runner script, one per line.

    .. code-block:: text

        python -m plumetools.cai2012ext.policy manifest.yaml          candidates
        python -m plumetools.cai2012ext.policy manifest.yaml --all    every case

    ``AllmeshCases`` and ``AllpostCases`` use ``--all``. Exit status 2 if the
    manifest contains an excluded Kn (nothing is listed), otherwise 0. Without
    ``--all``, the cases that are not production candidates are reported on
    stderr.
    """
    parser = argparse.ArgumentParser(description=main.__doc__.splitlines()[0])
    parser.add_argument("manifest", nargs="?",
                        help="manifest path; otherwise resolve manifest_name from --study")
    parser.add_argument("--study", type=Path,
                        help="study YAML whose manifest_name identifies the active suite")
    parser.add_argument("--all", action="store_true",
                        help="list every case (meshing and post-processing); "
                             "the Kn policy is still enforced")
    args = parser.parse_args(argv)
    if args.manifest is None and args.study is None:
        parser.error("provide a manifest or --study PATH")
    manifest_path = Path(args.manifest) if args.manifest is not None else None
    try:
        study = None
        if args.study is not None:
            from plumetools.cai2012.study import load_study
            study_path = args.study.resolve()
            study = load_study(study_path)
            check_study(study)
            if manifest_path is None:
                manifest_path = study_path.parent / study.manifest_name
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
        if study is not None and suite_name(manifest) != suite_name(study):
            raise ExtensionError("manifest suite differs from the selected study; regenerate it")
        runnable, refused = schedulable(manifest, allow_unsuitable=args.all)
    except (StudyError, OSError, yaml.YAMLError) as exc:
        print(f"policy: {exc}", file=sys.stderr)
        return 2
    for entry in runnable:
        print(entry["path"])
    for entry, reasons in refused:
        print(f"policy: not listed, {entry['path']} -- unsuitable for "
              f"production:", file=sys.stderr)
        for reason in reasons:
            print(f"            {reason}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
