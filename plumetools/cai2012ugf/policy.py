"""What the cai2012-unigas family must share with the kn-extension.

The family changes the solver, never the cases (spec §3): its template
``case.yaml`` equals the kn-extension's in every key outside ``meta``, and its
``study.yaml`` has the same case matrix, overrides, particle level and
multipliers. These functions name every difference; an empty list is the pass.
"""

from __future__ import annotations

from pathlib import Path

from plumetools.cai2012ext.policy import flatten_config

#: The kn-extension family, relative to the repository root.
KN_EXTENSION = Path("cases") / "cai2012-kn-extension"
FAMILY = Path("cases") / "cai2012-unigas"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def kn_extension_root() -> Path:
    return repo_root() / KN_EXTENSION


def template_differences(ugf_cfg, ext_cfg) -> list[str]:
    """Dotted ``case.yaml`` keys whose values differ, ignoring ``meta.*``."""
    mine, theirs = flatten_config(ugf_cfg), flatten_config(ext_cfg)
    return sorted(key for key in set(mine) | set(theirs)
                  if not key.startswith("meta.") and key != "meta"
                  and mine.get(key) != theirs.get(key))


def _cases(study) -> dict:
    return {c.name: c for c in study.kn_cases}


def _levels(study) -> dict:
    return {p.name: p for p in study.particle_levels}


def study_differences(ugf_study, ext_study) -> list[str]:
    """How the family's study matrix departs from the kn-extension's."""
    problems = []
    mine, theirs = _cases(ugf_study), _cases(ext_study)
    for name in sorted(set(mine) ^ set(theirs)):
        where = "this family" if name in mine else "the kn-extension"
        problems.append(f"kn_cases.{name}: only in {where}")
    for name in sorted(set(mine) & set(theirs)):
        a, b = mine[name], theirs[name]
        for attr in ("knudsen", "enabled", "overrides"):
            if (getattr(a, attr) or None) != (getattr(b, attr) or None):
                problems.append(f"kn_cases.{name}.{attr}: {getattr(a, attr)!r} "
                                f"vs kn-extension {getattr(b, attr)!r}")
    mine, theirs = _levels(ugf_study), _levels(ext_study)
    for name in sorted(set(mine) ^ set(theirs)):
        problems.append(f"particle_levels.{name}: not in both studies")
    for name in sorted(set(mine) & set(theirs)):
        for attr in ("multiplier", "enabled"):
            if getattr(mine[name], attr) != getattr(theirs[name], attr):
                problems.append(f"particle_levels.{name}.{attr} differs")
    for attr in ("run_time_multiplier", "output_frequency_multiplier", "base_case"):
        if getattr(ugf_study, attr) != getattr(ext_study, attr):
            problems.append(f"{attr}: {getattr(ugf_study, attr)!r} vs kn-extension "
                            f"{getattr(ext_study, attr)!r}")
    if (ugf_study.meta or {}).get("knudsen_suite") != (ext_study.meta or {}).get("knudsen_suite"):
        problems.append("meta.knudsen_suite differs")
    return problems
