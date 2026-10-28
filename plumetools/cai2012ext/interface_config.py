"""``interface.yaml``: the interface track's configuration, read and validated.

Split out of :mod:`plumetools.cai2012ext.interface`, which keeps the model --
the equivalent source, Bird's parameter, the sphere selection, the inlet faces
-- and needs nothing from here but :class:`InterfaceConfig` itself. This module
is the only one that reads the sidecar file; :func:`load_interface` is the one
place the mode (direct exit or interface) is decided.

Every name here is also importable from :mod:`plumetools.cai2012ext.interface`,
as it always was.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

#: The sidecar file, beside ``case.yaml``. Absent means the direct exit.
SETTINGS_FILE = "interface.yaml"

# --------------------------------------------------------------------------- #
# configuration: interface.yaml
# --------------------------------------------------------------------------- #

#: What a run is for. Only ``production`` can ever be a scientific result;
#: the other two are labelled so wherever their output goes.
PURPOSES = ("production", "interface-validation", "development-smoke")

#: ``production`` refuses to run with an unresolved criterion; ``exploratory``
#: runs anyway and is labelled NON-PRODUCTION for good.
RUN_MODES = ("production", "exploratory")

BIRD_DEFINITIONS = ("general", "radial_approx")

#: Only a sphere in this increment. An angle-dependent surface is out of scope.
SURFACES = ("sphere",)

#: ``exclude``: zero density beyond the core angle (the only setting a
#: validation can pass with). ``inject_exploratory``: the continuum tail is
#: injected too, which is a different model and is labelled as one.
TAIL_TREATMENTS = ("exclude", "inject_exploratory")



class InterfaceConfigError(ValueError):
    """``interface.yaml`` is malformed. Never ignored, never defaulted around."""


@dataclass(frozen=True)
class InterfaceConfig:
    """A validated ``interface.yaml`` ``interface:`` section.

    Attributes:
        enabled: ``False`` keeps the direct-exit case exactly as it was.
        purpose: one of :data:`PURPOSES`.
        run_mode: one of :data:`RUN_MODES`.
        surface: ``sphere``.
        bird_definition: ``general`` (default) or ``radial_approx``.
        bird_threshold: ``B_target``. A study choice, not a universal
            continuum-breakdown criterion.
        core_half_angle_deg: ``theta_c``: the cone in which the continuum state
            is required valid and is injected. **Required** for production and
            validation; no default is inferred from results.
        core_angle_provisional: ``True`` marks a value chosen without the
            Kn = 0.01 profiles, allowed only for development smoke runs.
        source_validity_min_radius_over_D: the smallest radius at which the
            point-source mapping is trusted. No paper certifies one; required
            for production and validation.
        source_validity_provisional: as ``core_angle_provisional``.
        source_offset_over_D: ``(x_s - x_e)/D``. Only ``(0, 0, 0)`` (CP-1) is
            implemented; the field exists so the source and the exit origin are
            never the same variable.
        angular_exponent_offset: the paper's ``0.41``.
        tails: one of :data:`TAIL_TREATMENTS`.
        radius_scan_points: radii sampled before the root is bracketed.
        theta_scan_points: polar angles sampled across the core.
        domain_clearance_over_D: least gap between the sphere and any domain
            face it does not sit on.
        diagnostic_radius_over_D: the angular-comparison sphere (Allcompare);
            10 reproduces Cai's r/D = 10, anything else is recorded as a change.
        diagnostic_clearance_over_D: least gap between the interface and that
            sphere.
        meshed_flux_tolerance: relative mismatch allowed between the gross flux
            through the meshed patch and the core target before the coupling
            is called failed (faceting plus the gross/net difference).
        validation_evidence: path, relative to the case, of the
            ``overlap.yaml`` that validates this threshold; production only.
        label: free text carried into every report.
    """

    enabled: bool = False
    purpose: str = "production"
    run_mode: str = "production"
    surface: str = "sphere"
    bird_definition: str = "general"
    bird_threshold: float = 0.05
    core_half_angle_deg: float | None = None
    core_angle_provisional: bool = False
    source_validity_min_radius_over_D: float | None = None
    source_validity_provisional: bool = False
    source_offset_over_D: tuple = (0.0, 0.0, 0.0)
    angular_exponent_offset: float = 0.41
    tails: str = "exclude"
    radius_scan_points: int = 400
    theta_scan_points: int = 91
    domain_clearance_over_D: float = 0.5
    diagnostic_radius_over_D: float = 10.0
    diagnostic_clearance_over_D: float = 0.5
    meshed_flux_tolerance: float = 0.02
    validation_evidence: str | None = None
    label: str = ""

    @property
    def scientific(self) -> bool:
        """Whether the run is meant to produce evidence (not a smoke test)."""
        return self.purpose != "development-smoke"

    @property
    def core_half_angle_rad(self) -> float | None:
        if self.core_half_angle_deg is None:
            return None
        return math.radians(float(self.core_half_angle_deg))

    def missing_decisions(self) -> list[str]:
        """Choices a scientific run needs and this file has not made."""
        missing = []
        if self.core_half_angle_deg is None:
            missing.append(
                "core_half_angle_deg is null. No Kn = 0.01 profile is available "
                "locally to justify a plume-core angle, so a scientific run "
                "needs an explicit value (docs/cai2012-kn-interface.md, section 4).")
        elif self.scientific and self.core_angle_provisional:
            missing.append(
                "core_half_angle_deg is marked provisional; provisional values "
                "are for development-smoke runs only.")
        if self.source_validity_min_radius_over_D is None:
            missing.append(
                "source_validity_min_radius_over_D is null. No paper certifies a "
                "minimum radius for this mapping; a scientific run needs an "
                "explicit, justified value.")
        elif self.scientific and self.source_validity_provisional:
            missing.append(
                "source_validity_min_radius_over_D is marked provisional; "
                "provisional values are for development-smoke runs only.")
        return missing

    def as_dict(self) -> dict:
        return {
            "enabled": bool(self.enabled), "purpose": self.purpose,
            "run_mode": self.run_mode, "surface": self.surface,
            "bird_definition": self.bird_definition,
            "bird_threshold": float(self.bird_threshold),
            "core_half_angle_deg": None if self.core_half_angle_deg is None
            else float(self.core_half_angle_deg),
            "core_angle_provisional": bool(self.core_angle_provisional),
            "source_validity_min_radius_over_D":
                None if self.source_validity_min_radius_over_D is None
                else float(self.source_validity_min_radius_over_D),
            "source_validity_provisional": bool(self.source_validity_provisional),
            "source_offset_over_D": [float(v) for v in self.source_offset_over_D],
            "angular_exponent_offset": float(self.angular_exponent_offset),
            "tails": self.tails,
            "radius_scan_points": int(self.radius_scan_points),
            "theta_scan_points": int(self.theta_scan_points),
            "domain_clearance_over_D": float(self.domain_clearance_over_D),
            "diagnostic_radius_over_D": float(self.diagnostic_radius_over_D),
            "diagnostic_clearance_over_D": float(self.diagnostic_clearance_over_D),
            "meshed_flux_tolerance": float(self.meshed_flux_tolerance),
            "validation_evidence": self.validation_evidence,
            "label": self.label,
        }


#: Top-level sections ``interface.yaml`` may contain. ``overlap`` is read by
#: :mod:`plumetools.cai2012ext.overlap`, ``resources`` by
#: :mod:`plumetools.cai2012ext.interface_case`; each validates its own.
SECTIONS = ("interface", "overlap", "resources")


def _number(value, key: str, *, positive=False, minimum=None, maximum=None,
            allow_none=False):
    if value is None:
        if allow_none:
            return None
        raise InterfaceConfigError(f"interface.{key} must be a number, not null")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InterfaceConfigError(
            f"interface.{key} must be a number, got {value!r}")
    value = float(value)
    if not math.isfinite(value):
        raise InterfaceConfigError(f"interface.{key} must be finite, got {value}")
    if positive and not value > 0.0:
        raise InterfaceConfigError(f"interface.{key} must be positive, got {value}")
    if minimum is not None and value < minimum:
        raise InterfaceConfigError(
            f"interface.{key} must be >= {minimum}, got {value}")
    if maximum is not None and value > maximum:
        raise InterfaceConfigError(
            f"interface.{key} must be <= {maximum}, got {value}")
    return value


def _choice(value, key: str, choices) -> str:
    if value not in choices:
        raise InterfaceConfigError(
            f"interface.{key} is {value!r}; valid: {list(choices)}")
    return str(value)


def _flag(value, key: str) -> bool:
    if not isinstance(value, bool):
        raise InterfaceConfigError(
            f"interface.{key} must be true or false, got {value!r}")
    return value


def parse_interface(raw: dict, where: str = SETTINGS_FILE) -> InterfaceConfig:
    """Validate the ``interface:`` mapping. Unknown keys are errors."""
    if not isinstance(raw, dict):
        raise InterfaceConfigError(f"{where}: 'interface' must be a mapping")
    known = {f.name for f in InterfaceConfig.__dataclass_fields__.values()}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise InterfaceConfigError(
            f"{where}: unknown interface key(s) {unknown}; valid: {sorted(known)}")
    if "enabled" not in raw:
        raise InterfaceConfigError(
            f"{where}: interface.enabled must be stated (true or false); a "
            f"sidecar that does not say whether it is on is ambiguous")
    values = dict(raw)
    out = {"enabled": _flag(values["enabled"], "enabled")}
    out["purpose"] = _choice(values.get("purpose", "production"), "purpose", PURPOSES)
    out["run_mode"] = _choice(values.get("run_mode", "production"), "run_mode",
                              RUN_MODES)
    out["surface"] = _choice(values.get("surface", "sphere"), "surface", SURFACES)
    out["bird_definition"] = _choice(values.get("bird_definition", "general"),
                                     "bird_definition", BIRD_DEFINITIONS)
    out["bird_threshold"] = _number(values.get("bird_threshold", 0.05),
                                    "bird_threshold", positive=True)
    out["core_half_angle_deg"] = _number(
        values.get("core_half_angle_deg"), "core_half_angle_deg",
        positive=True, maximum=90.0, allow_none=True)
    out["core_angle_provisional"] = _flag(
        values.get("core_angle_provisional", False), "core_angle_provisional")
    out["source_validity_min_radius_over_D"] = _number(
        values.get("source_validity_min_radius_over_D"),
        "source_validity_min_radius_over_D", positive=True, allow_none=True)
    out["source_validity_provisional"] = _flag(
        values.get("source_validity_provisional", False),
        "source_validity_provisional")
    offset = values.get("source_offset_over_D", [0.0, 0.0, 0.0])
    if (not isinstance(offset, (list, tuple)) or len(offset) != 3):
        raise InterfaceConfigError(
            f"{where}: interface.source_offset_over_D must be three numbers")
    offset = tuple(_number(v, "source_offset_over_D") for v in offset)
    if any(v != 0.0 for v in offset):
        raise InterfaceConfigError(
            f"{where}: interface.source_offset_over_D is {list(offset)}. Only "
            f"the exit centre (0, 0, 0) is implemented in this increment "
            f"(assumption CP-1); a displaced source needs its own mapping.")
    out["source_offset_over_D"] = offset
    out["angular_exponent_offset"] = _number(
        values.get("angular_exponent_offset", 0.41), "angular_exponent_offset")
    out["tails"] = _choice(values.get("tails", "exclude"), "tails",
                           TAIL_TREATMENTS)
    for key, default, lo in (("radius_scan_points", 400, 16),
                             ("theta_scan_points", 91, 3)):
        value = values.get(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or value < lo:
            raise InterfaceConfigError(
                f"{where}: interface.{key} must be an integer >= {lo}, "
                f"got {value!r}")
        out[key] = int(value)
    out["domain_clearance_over_D"] = _number(
        values.get("domain_clearance_over_D", 0.5), "domain_clearance_over_D",
        minimum=0.0)
    out["diagnostic_radius_over_D"] = _number(
        values.get("diagnostic_radius_over_D", 10.0), "diagnostic_radius_over_D",
        positive=True)
    out["diagnostic_clearance_over_D"] = _number(
        values.get("diagnostic_clearance_over_D", 0.5),
        "diagnostic_clearance_over_D", minimum=0.0)
    out["meshed_flux_tolerance"] = _number(
        values.get("meshed_flux_tolerance", 0.02), "meshed_flux_tolerance",
        positive=True, maximum=1.0)
    evidence = values.get("validation_evidence")
    if evidence is not None and not isinstance(evidence, str):
        raise InterfaceConfigError(
            f"{where}: interface.validation_evidence must be a path or null")
    out["validation_evidence"] = evidence
    label = values.get("label", "")
    out["label"] = "" if label is None else str(label)

    cfg = InterfaceConfig(**out)
    if cfg.purpose == "development-smoke" and cfg.run_mode == "production":
        raise InterfaceConfigError(
            f"{where}: a development-smoke run cannot have run_mode production")
    if cfg.purpose != "development-smoke" and cfg.tails != "exclude":
        raise InterfaceConfigError(
            f"{where}: tails '{cfg.tails}' injects the continuum tail outside "
            f"the core, which is a different model; it is allowed only for "
            f"development-smoke runs, never for {cfg.purpose}")
    return cfg


def read_settings(case_dir: Path) -> dict | None:
    """The raw ``interface.yaml`` of a case, or ``None`` if there is none."""
    path = Path(case_dir) / SETTINGS_FILE
    if not path.is_file():
        return None
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise InterfaceConfigError(f"{path}: not valid YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise InterfaceConfigError(f"{path}: must be a mapping")
    unknown = sorted(set(raw) - set(SECTIONS))
    if unknown:
        raise InterfaceConfigError(
            f"{path}: unknown section(s) {unknown}; valid: {list(SECTIONS)}")
    if "interface" not in raw:
        raise InterfaceConfigError(f"{path}: no 'interface' section")
    return raw


def load_interface(case_dir: Path) -> InterfaceConfig | None:
    """The case's interface settings, or ``None`` for the direct exit.

    ``None`` when ``interface.yaml`` is absent or says ``enabled: false``; a
    malformed file raises. This is THE mode selection: every script and module
    that has to know whether a case injects at the exit or at a sphere asks
    here.
    """
    raw = read_settings(case_dir)
    if raw is None:
        return None
    cfg = parse_interface(raw["interface"], str(Path(case_dir) / SETTINGS_FILE))
    return cfg if cfg.enabled else None
