r"""``case.yaml`` for the nozzle-plume family: schema, loading, validation.

The file is meant to be read and edited by hand, so it is kept to what a person
would actually know about a case: the nozzle exit and the state on it, the gas,
the plate, the box, how fine the mesh is, and how the DSMC is run. Everything
else is derived and printed, never typed in.

Three rules shape the schema.

**No redundant keys.** An exit diameter is an input; the exit area is printed
next to it. A molecular diameter is an input; the collision cross-section is
printed. Two keys for one fact are two chances for them to disagree, and the
loader would then have to pick one.

**Alternatives are spelled as alternatives.** The exit density may be given as a
number density, a mass density, or a Knudsen number -- exactly one of them. The
exit speed as a speed or a speed ratio -- exactly one. Giving two is an error
rather than a silent precedence rule.

**``null`` means "derive it".** The particle weight, the time step and the run
times are all derived when left ``null``, and all of them can be pinned by
writing a number instead. ``dsmc.n_equivalent_particles`` is the one most worth
pinning: the weight is the cost of the run, and for a dense exit it is the first
thing that has to give.

Unknown keys are an error everywhere, so a typo in a physical parameter stops
the run instead of silently falling back to a default.
"""

from __future__ import annotations

import math
import re
from dataclasses import MISSING, dataclass, field, fields, replace
from pathlib import Path
from typing import Any

import yaml

from plumetools.nozzleplume import MODEL
from plumetools.nozzleplume.gas import (
    MEAN_FREE_PATH_CONVENTIONS,
    SPECIES_NAME,
    VhsSpecies,
)
from plumetools.rarefaction.settings import DiagnosticSettings, SettingsError


class ConfigError(ValueError):
    """``case.yaml`` is missing, malformed, or internally inconsistent."""


#: What a refinement region may be called. The name is a key in
#: ``snappyHexMeshDict``'s ``geometry`` dictionary.
_REGION_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")

#: Axes a symmetry plane may be normal to. ``y`` is the plane ``y = 0``; the
#: modelled domain is then ``y >= 0``. The plume axis ``x`` is never one.
SYMMETRY_AXES = ("y", "z")

#: Binary collision models this family writes coefficients for.
#:
#: ``VariableHardSphere``                 elastic VHS. The default.
#: ``LarsenBorgnakkeVariableHardSphere``  VHS with internal-energy exchange; the
#:                                        one to use with internal modes.
#: ``NoBinaryCollision``                  collisionless. A diagnostic, not a model
#:                                        of anything dense.
COLLISION_MODELS = ("VariableHardSphere", "LarsenBorgnakkeVariableHardSphere",
                    "NoBinaryCollision")

#: Gas-surface models for the plate.
#:
#: ``MaxwellianThermal``     fully diffuse, fully accommodated to the wall
#:                           temperature. The default.
#: ``SpecularReflection``    mirror reflection; no energy exchange.
#: ``MixedDiffuseSpecular``  ``dsmc.diffuse_fraction`` diffuse, the rest specular.
WALL_MODELS = ("MaxwellianThermal", "SpecularReflection", "MixedDiffuseSpecular")

#: ``decomposePar`` methods this family writes.
DECOMPOSITION_METHODS = ("scotch", "simple")

#: Shapes a refinement region may take.
#:
#: ``cone``  coaxial with the plume: ``radius1_m`` at ``x_min_m``, ``radius2_m``
#:           at ``x_max_m``, linear between. A ``searchableCone``.
#: ``box``   axis-aligned, all six bounds given.
REGION_SHAPES = ("cone", "box")

#: Keys that belong to one shape and are an error on the other.
_CONE_KEYS = ("radius1_m", "radius2_m", "inner_radius1_m", "inner_radius2_m")
_BOX_KEYS = ("y_min_m", "y_max_m", "z_min_m", "z_max_m")

#: Patch roles and their default names.
DEFAULT_PATCH_NAMES = {
    "nozzle": "nozzle",
    "upstream_vacuum": "upstreamVacuum",
    "outer": "vacuum",
    "plate": "plate",
    "symmetry": "symmetry",
}


@dataclass(frozen=True)
class NozzleConfig:
    """The circular exit.

    Attributes:
        exit_diameter_m: ``D`` [m]. The radius and the exit area are derived.
    """

    exit_diameter_m: float

    @property
    def radius_m(self) -> float:
        """``R0 = D/2`` [m]."""
        return 0.5 * float(self.exit_diameter_m)

    @property
    def area_m2(self) -> float:
        """``pi R0^2`` [m^2], the exact area the meshed patch approximates."""
        return math.pi * self.radius_m ** 2


@dataclass(frozen=True)
class ExitConfig:
    """The uniform drifting Maxwellian on the exit disk.

    Attributes:
        temperature_K: exit static temperature [K].
        speed_m_per_s: bulk speed along ``+x`` [m/s].
        speed_ratio: ``S = U / sqrt(2 R T)``, as an alternative to the speed.
        number_density_per_m3: exit number density [1/m^3].
        mass_density_kg_per_m3: exit mass density [kg/m^3], as an alternative.
        knudsen: exit Knudsen number on the exit diameter, as an alternative.

    **Exactly one** of the two speed keys and **exactly one** of the three
    density keys must be set; the others stay ``null``.
    """

    temperature_K: float
    speed_m_per_s: float | None = None
    speed_ratio: float | None = None
    number_density_per_m3: float | None = None
    mass_density_kg_per_m3: float | None = None
    knudsen: float | None = None


@dataclass(frozen=True)
class GasConfig:
    """The single VHS species. See :class:`plumetools.nozzleplume.gas.VhsSpecies`.

    Attributes:
        molecular_mass_kg: molecular mass [kg]. For a mixture treated as one
            species, the mean molecular mass.
        molecular_diameter_m: VHS reference diameter [m].
        species_name: the ``typeId`` in ``dsmcProperties``.
        omega: VHS exponent. 0.5 -- a hard sphere -- makes the collision
            cross-section ``pi d^2`` at every temperature.
        t_ref_K: reference temperature of the diameter [K]; inert at
            ``omega = 0.5``.
        internal_degrees_of_freedom: internal modes carried by each molecule.
        mean_free_path_convention: ``vhs`` or ``hard_sphere``; identical at
            ``omega = 0.5``.
    """

    molecular_mass_kg: float
    molecular_diameter_m: float
    species_name: str = "gas"
    omega: float = 0.5
    t_ref_K: float = 273.0
    internal_degrees_of_freedom: int = 0
    mean_free_path_convention: str = "vhs"

    def species(self) -> VhsSpecies:
        """This section as a :class:`~plumetools.nozzleplume.gas.VhsSpecies`."""
        return VhsSpecies(
            name=str(self.species_name),
            mass_kg=float(self.molecular_mass_kg),
            diameter_m=float(self.molecular_diameter_m),
            omega=float(self.omega),
            t_ref_K=float(self.t_ref_K),
            internal_degrees_of_freedom=int(self.internal_degrees_of_freedom),
        )


@dataclass(frozen=True)
class PlateConfig:
    """A flat rectangular plate, normal to the plume axis.

    The section is optional: leave it out, or write ``plate: null``, for a free
    plume into vacuum.

    Attributes:
        standoff_m: exit plane to the plate's **upstream face** [m].
        thickness_m: extent along ``x``, the plume axis [m].
        size_y_m, size_z_m: full extents along ``y`` and ``z`` [m].
        centre_y_m, centre_z_m: where the plate's centre sits off the plume
            axis [m]. Zero centres the plate on the plume.
        temperature_K: wall temperature [K]. A diffuse wall re-emits at it.
    """

    standoff_m: float
    thickness_m: float
    size_y_m: float
    size_z_m: float
    centre_y_m: float = 0.0
    centre_z_m: float = 0.0
    temperature_K: float = 300.0


@dataclass(frozen=True)
class DomainConfig:
    """The vacuum box. The exit plane ``x = 0`` is its upstream face.

    Attributes:
        x_max_m: downstream boundary [m].
        y_half_m, z_half_m: lateral half-extents [m].
        symmetry_planes: any of ``y`` and ``z``. ``[y]`` models ``y >= 0`` with
            a symmetry plane at ``y = 0``; ``[y, z]`` models a quarter. Empty --
            the default -- is the full 3-D box.

    Each extent is rounded **up** to a whole number of background cells, so the
    background mesh is exactly cubic and ``y = 0`` and ``z = 0`` are always
    lattice planes. The rounding is reported.
    """

    x_max_m: float
    y_half_m: float
    z_half_m: float
    symmetry_planes: tuple = ()


@dataclass(frozen=True)
class RegionConfig:
    """One volume refined to a fixed octree level.

    Attributes:
        name: the ``geometry`` entry name in ``snappyHexMeshDict``.
        level: octree level; the cell size is
            ``mesh.background_cell_size_m / 2**level``.
        x_max_m: downstream end [m].
        shape: ``cone`` (coaxial with the plume) or ``box``.
        x_min_m: upstream end [m]. Zero is the exit plane.
        radius1_m: cone radius at ``x_min_m`` [m].
        radius2_m: cone radius at ``x_max_m`` [m].
        inner_radius1_m, inner_radius2_m: inner radii at the two ends [m], for a
            hollow cone. ``null`` -- the default -- is solid.
        y_min_m, y_max_m, z_min_m, z_max_m: box bounds [m].

    A cone is written as a ``searchableCone`` with ``point1 = (x_min_m 0 0)``
    and ``point2 = (x_max_m 0 0)``, so the ``1`` and ``2`` of the radii are the
    ``1`` and ``2`` of that dictionary: upstream end, downstream end. The radius
    is linear in ``x`` between them, which is the shape of a plume: one that
    widens downstream follows it with a fraction of the cells a cylinder wide
    enough for its far end would spend on its near one.

    Equal radii are a cylinder, and that is the only way one is written.
    """

    name: str
    level: int
    x_max_m: float
    shape: str = "cone"
    x_min_m: float = 0.0
    radius1_m: float | None = None
    radius2_m: float | None = None
    inner_radius1_m: float | None = None
    inner_radius2_m: float | None = None
    y_min_m: float | None = None
    y_max_m: float | None = None
    z_min_m: float | None = None
    z_max_m: float | None = None


@dataclass(frozen=True)
class MeshConfig:
    """A uniform cubic background, refined by ``snappyHexMesh``.

    Attributes:
        background_cell_size_m: the level-0 cell [m]. Every other cell size in
            the mesh is this halved a whole number of times.
        plate_refinement_level: octree level at the plate surface.
        refinement_regions: volumes refined to a fixed level; see
            :class:`RegionConfig`. The exit disk has to sit inside one of them
            (or at level 0), entirely at one level.
        n_cells_between_levels: buffer cells between successive levels.
        max_cells: refuse to write a mesh *estimated* above this. An octree's
            cell count is not known until it is built, so this is a guard
            against a mistyped level, not a budget met to the cell.
        patch_names: role -> mesh patch name.
    """

    background_cell_size_m: float
    plate_refinement_level: int = 2
    refinement_regions: tuple = ()
    n_cells_between_levels: int = 3
    max_cells: int = 3_000_000
    patch_names: dict = field(default_factory=lambda: dict(DEFAULT_PATCH_NAMES))


@dataclass(frozen=True)
class DsmcConfig:
    """Solver settings written into ``dsmcProperties`` and ``controlDict``.

    Attributes:
        n_equivalent_particles: real molecules per simulated parcel. ``null`` ->
            derived so the exit cell holds ``target_particles_per_exit_cell``.
            **Written as given when set**: this is the knob that trades
            statistical quality for run time, and for a dense exit the derived
            value is rarely affordable.
        target_particles_per_exit_cell: the occupancy the derivation aims at.
            Unused when the weight is pinned.
        collision_model: one of :data:`COLLISION_MODELS`.
        relaxation_collision_number: Larsen-Borgnakke ``Z``; used only by that
            model.
        wall_interaction_model: one of :data:`WALL_MODELS`.
        diffuse_fraction: used only by ``MixedDiffuseSpecular``.
        delta_t_s: time step [s]. ``null`` -> derived from ``courant_target``.
        courant_target: fraction of the smallest cell a fast molecule may cross
            per step.
        end_time_s, average_start_s, write_interval_s: ``null`` -> derived from
            the two durations below.
        transient_domain_transits: run this many transits of the domain at the
            exit speed before averaging starts.
        sampling_domain_transits: then average over this many.
        initial_number_density_per_m3, initial_temperature_K: the near-vacuum
            fill ``dsmcInitialise`` writes. Not a physical state.
        n_subdomains: MPI ranks for ``decomposeParDict``.
        decomposition_method: ``scotch`` or ``simple``.
    """

    n_equivalent_particles: float | None = None
    target_particles_per_exit_cell: float = 20.0
    collision_model: str = "VariableHardSphere"
    relaxation_collision_number: float = 5.0
    wall_interaction_model: str = "MaxwellianThermal"
    diffuse_fraction: float = 1.0
    delta_t_s: float | None = None
    courant_target: float = 0.2
    end_time_s: float | None = None
    average_start_s: float | None = None
    write_interval_s: float | None = None
    transient_domain_transits: float = 4.0
    sampling_domain_transits: float = 2.0
    initial_number_density_per_m3: float = 1.0e+10
    initial_temperature_K: float = 300.0
    n_subdomains: int = 4
    decomposition_method: str = "scotch"


@dataclass(frozen=True)
class ChecksConfig:
    """DSMC quality thresholds. See :mod:`plumetools.nozzleplume.checks`.

    Attributes:
        enabled: run the checks at all.
        max_cell_over_mfp: exit cell over exit mean free path; a warning above.
        max_time_step_over_collision_time: ``deltaT`` over the exit mean
            collision time; a warning above.
        max_courant: a hard error above -- molecules would skip cells.
        warn_courant: a warning above.
        min_particles_per_cell: a warning below, in the exit cell.
        min_nozzle_faces: a hard error below -- the exit would be a square.
        max_nozzle_area_error: a hard error above; the mass flow scales with it.
        warn_parcels: a warning above this many estimated parcels.
    """

    enabled: bool = True
    max_cell_over_mfp: float = 1.0
    max_time_step_over_collision_time: float = 1.0
    max_courant: float = 1.0
    warn_courant: float = 0.5
    min_particles_per_cell: float = 5.0
    min_nozzle_faces: int = 12
    max_nozzle_area_error: float = 0.05
    warn_parcels: float = 3.0e+7


@dataclass(frozen=True)
class PostConfig:
    """Post-processing.

    Attributes:
        centreline_radius_over_D: cells within this radius of the axis are
            averaged into a centreline sample. A single cell column is noise; a
            fat tube is not the centreline. Where the cells are larger than the
            tube it is widened to the four cells round the axis, so the profile
            does not stop where the mesh coarsens.
        centreline_points: bins along the axis.
    """

    centreline_radius_over_D: float = 0.25
    centreline_points: int = 200


@dataclass(frozen=True)
class CaseConfig:
    """A whole validated ``case.yaml``."""

    nozzle: NozzleConfig
    exit: ExitConfig
    gas: GasConfig
    domain: DomainConfig
    mesh: MeshConfig
    plate: PlateConfig | None = None
    dsmc: DsmcConfig = field(default_factory=DsmcConfig)
    checks: ChecksConfig = field(default_factory=ChecksConfig)
    post: PostConfig = field(default_factory=PostConfig)
    #: Local mean free path, Knudsen numbers and mesh resolution after a run;
    #: see :mod:`plumetools.rarefaction`. The exit diameter is always one of
    #: its reference lengths, as ``D``.
    diagnostics: DiagnosticSettings = field(default_factory=DiagnosticSettings)
    model: str = MODEL
    #: Free-form provenance. Not validated: it records where a case came from,
    #: and over-constraining that makes recording an awkward fact harder than
    #: omitting it.
    meta: dict = field(default_factory=dict)

    def species(self) -> VhsSpecies:
        """The VHS species, as the gas section describes it."""
        return self.gas.species()


#: Sections every ``case.yaml`` must have, and the ones it may leave out.
_REQUIRED_SECTIONS = {
    "nozzle": NozzleConfig,
    "exit": ExitConfig,
    "gas": GasConfig,
    "domain": DomainConfig,
    "mesh": MeshConfig,
}
_OPTIONAL_SECTIONS = {
    "dsmc": DsmcConfig,
    "checks": ChecksConfig,
    "post": PostConfig,
}
_TOP_LEVEL = (set(_REQUIRED_SECTIONS) | set(_OPTIONAL_SECTIONS)
              | {"plate", "model", "meta", "diagnostics"})


def _build(cls, data, where: str):
    """Instantiate a section dataclass from a mapping.

    Unknown keys and missing required keys are both errors, each naming the
    section, so a typo in a physical parameter stops the run instead of silently
    falling back to a default.

    Numeric-looking strings are coerced, because YAML 1.1 parses ``1.0e14`` as a
    **string** (it requires a signed exponent) while ``1.0e+14`` is a float. The
    failure that causes without this is an unformattable value deep inside a
    dictionary renderer, a long way from the line responsible.
    """
    if not isinstance(data, dict):
        raise ConfigError(
            f"{where}: expected a mapping, got {type(data).__name__}")
    known = {f.name for f in fields(cls)}
    unknown = sorted(set(data) - known)
    if unknown:
        raise ConfigError(
            f"{where}: unknown key(s) {unknown}; valid keys are {sorted(known)}")
    required = sorted(
        f.name for f in fields(cls)
        if f.default is MISSING and f.default_factory is MISSING
        and data.get(f.name) is None)
    if required:
        raise ConfigError(f"{where}: missing required key(s) {required}")

    coerced: dict[str, Any] = {}
    for f in fields(cls):
        if f.name not in data:
            continue
        value = data[f.name]
        annotation = str(f.type)
        if value is None:
            coerced[f.name] = None
            continue
        if isinstance(value, bool) and annotation != "bool":
            raise ConfigError(
                f"{where}: {f.name} is {value!r}, a boolean, where a "
                f"{annotation.split(' ')[0]} is expected")
        if annotation.startswith("tuple") and isinstance(value, list):
            value = tuple(value)
        elif isinstance(value, str) and ("float" in annotation or "int" in annotation):
            try:
                value = float(value) if "float" in annotation else int(value)
            except ValueError:
                raise ConfigError(
                    f"{where}: {f.name} is {value!r}, which is not a number. If "
                    f"it looks like one, check the exponent: YAML requires a "
                    f"signed exponent, so 1.0e+14 is a float and 1.0e14 is a "
                    f"string.") from None
        coerced[f.name] = value
    return cls(**coerced)


def _build_mesh(data, where: str) -> MeshConfig:
    """The mesh section, with its region list and its patch-name overrides."""
    if not isinstance(data, dict):
        raise ConfigError(f"{where}: expected a mapping, got {type(data).__name__}")
    data = dict(data)

    regions = data.pop("refinement_regions", None) or []
    if not isinstance(regions, list):
        raise ConfigError(
            f"{where}: refinement_regions must be a list of mappings, got "
            f"{type(regions).__name__}")
    for i, region in enumerate(regions):
        # Named before the generic unknown-key error can be: a case.yaml from
        # before the cones would otherwise be told only that 'radius_m' is not
        # a key, and not what replaced it.
        if isinstance(region, dict) and (
                region.get("shape") == "cylinder"
                or ("radius_m" in region and "shape" not in region)):
            raise ConfigError(
                f"{where}:refinement_regions[{i}]: this region is written as a "
                f"cylinder (radius_m, or shape: cylinder). The axisymmetric "
                f"region is a cone now: give radius1_m, its radius at x_min_m, "
                f"and radius2_m, its radius at x_max_m. The same number for "
                f"both is the cylinder this was.")
    built = tuple(_build(RegionConfig, region, f"{where}:refinement_regions[{i}]")
                  for i, region in enumerate(regions))

    # Overrides are merged into the defaults, so a case naming one patch
    # differently does not have to restate the other four.
    names = data.pop("patch_names", None)
    if names is not None and not isinstance(names, dict):
        raise ConfigError(f"{where}: patch_names must be a mapping")
    unknown = sorted(set(names or {}) - set(DEFAULT_PATCH_NAMES))
    if unknown:
        raise ConfigError(
            f"{where}: unknown patch role(s) {unknown}; valid roles are "
            f"{sorted(DEFAULT_PATCH_NAMES)}")
    patch_names = {**DEFAULT_PATCH_NAMES, **(names or {})}

    return replace(_build(MeshConfig, data, where),
                   refinement_regions=built, patch_names=patch_names)


def load_case_config(case_dir: Path) -> CaseConfig:
    """Read and validate ``<case_dir>/case.yaml``.

    Args:
        case_dir: a case directory, or the ``case.yaml`` file itself.

    Returns:
        A validated :class:`CaseConfig`.

    Raises:
        ConfigError: for a missing file, an unknown or missing key, the wrong
            ``model``, or any failed cross-check in :func:`validate`.
    """
    path = Path(case_dir)
    if path.is_dir():
        path = path / "case.yaml"
    if not path.is_file():
        raise ConfigError(f"no case.yaml at {path}")

    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path}: not valid YAML: {exc}") from None
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: expected a mapping at the top level")

    unknown = sorted(set(raw) - _TOP_LEVEL)
    if unknown:
        raise ConfigError(
            f"{path}: unknown top-level key(s) {unknown}; valid: {sorted(_TOP_LEVEL)}")

    model = raw.get("model")
    if model != MODEL:
        raise ConfigError(
            f"{path}: model is {model!r}, expected {MODEL!r}. This loader reads "
            f"the nozzle-plume schema only; a cai2012_circular_plume case.yaml "
            f"belongs to plumetools.cai2012.config, and a source-flow one to "
            f"plumetools.config.")

    missing = sorted(name for name in _REQUIRED_SECTIONS if raw.get(name) is None)
    if missing:
        raise ConfigError(f"{path}: missing required section(s) {missing}")

    sections: dict[str, Any] = {}
    for name, cls in _REQUIRED_SECTIONS.items():
        where = f"{path}:{name}"
        sections[name] = (_build_mesh(raw[name], where) if name == "mesh"
                          else _build(cls, raw[name], where))
    for name, cls in _OPTIONAL_SECTIONS.items():
        sections[name] = _build(cls, raw.get(name) or {}, f"{path}:{name}")

    try:
        sections["diagnostics"] = DiagnosticSettings.from_mapping(
            raw.get("diagnostics"), f"{path}:diagnostics")
    except SettingsError as exc:
        raise ConfigError(str(exc)) from None

    plate = raw.get("plate")
    cfg = CaseConfig(
        model=model,
        meta=raw.get("meta") or {},
        plate=None if plate is None else _build(PlateConfig, plate, f"{path}:plate"),
        **sections,
    )
    validate(cfg, path)
    return cfg


def _is_whole(value) -> bool:
    return float(value).is_integer()


def validate(cfg: CaseConfig, path: Path | str = "case.yaml") -> None:
    """Cross-check a config. Called by :func:`load_case_config`.

    Every check here would otherwise surface as a meshing failure, a solver
    abort with no indication of which key caused it, or -- worse -- a case that
    runs cleanly at the wrong physical conditions. Checks that need the derived
    geometry or the mesh plan live with those, in
    :func:`plumetools.nozzleplume.geometry.from_config` and
    :func:`plumetools.nozzleplume.mesh.plan`.

    Raises:
        ConfigError: on the first inconsistency, naming the key and the fix.
    """
    def fail(message: str):
        raise ConfigError(f"{path}: {message}")

    # --- exactly-one-of -----------------------------------------------------
    for label, keys in (
            ("speed", ("speed_m_per_s", "speed_ratio")),
            ("density", ("number_density_per_m3", "mass_density_kg_per_m3",
                         "knudsen"))):
        given = [k for k in keys if getattr(cfg.exit, k) is not None]
        if len(given) != 1:
            fail(f"the exit {label} must be given by exactly one of "
                 f"{['exit.' + k for k in keys]}; "
                 + (f"{['exit.' + k for k in given]} are all set. Keep one and "
                    f"set the others to null -- the rest are printed as derived "
                    f"values, so nothing is lost." if given
                    else "none is set."))

    # --- positive physical inputs -------------------------------------------
    positives = {
        "nozzle.exit_diameter_m": cfg.nozzle.exit_diameter_m,
        "exit.temperature_K": cfg.exit.temperature_K,
        "gas.molecular_mass_kg": cfg.gas.molecular_mass_kg,
        "gas.molecular_diameter_m": cfg.gas.molecular_diameter_m,
        "gas.t_ref_K": cfg.gas.t_ref_K,
        "domain.x_max_m": cfg.domain.x_max_m,
        "domain.y_half_m": cfg.domain.y_half_m,
        "domain.z_half_m": cfg.domain.z_half_m,
        "mesh.background_cell_size_m": cfg.mesh.background_cell_size_m,
        "dsmc.target_particles_per_exit_cell":
            cfg.dsmc.target_particles_per_exit_cell,
        "dsmc.courant_target": cfg.dsmc.courant_target,
        "dsmc.transient_domain_transits": cfg.dsmc.transient_domain_transits,
        "dsmc.sampling_domain_transits": cfg.dsmc.sampling_domain_transits,
        "dsmc.initial_number_density_per_m3":
            cfg.dsmc.initial_number_density_per_m3,
        "dsmc.initial_temperature_K": cfg.dsmc.initial_temperature_K,
        "dsmc.relaxation_collision_number": cfg.dsmc.relaxation_collision_number,
        "checks.warn_parcels": cfg.checks.warn_parcels,
        "post.centreline_radius_over_D": cfg.post.centreline_radius_over_D,
    }
    if cfg.plate is not None:
        positives.update({
            "plate.standoff_m": cfg.plate.standoff_m,
            "plate.thickness_m": cfg.plate.thickness_m,
            "plate.size_y_m": cfg.plate.size_y_m,
            "plate.size_z_m": cfg.plate.size_z_m,
            "plate.temperature_K": cfg.plate.temperature_K,
        })
    bad = sorted(name for name, value in positives.items()
                 if value is None or not math.isfinite(float(value))
                 or float(value) <= 0.0)
    if bad:
        fail(f"{bad} must be finite and positive")

    # Optional numbers: null derives, anything else has to be usable.
    for name, value in (
            ("exit.speed_m_per_s", cfg.exit.speed_m_per_s),
            ("exit.speed_ratio", cfg.exit.speed_ratio),
            ("exit.number_density_per_m3", cfg.exit.number_density_per_m3),
            ("exit.mass_density_kg_per_m3", cfg.exit.mass_density_kg_per_m3),
            ("exit.knudsen", cfg.exit.knudsen),
            ("dsmc.n_equivalent_particles", cfg.dsmc.n_equivalent_particles),
            ("dsmc.delta_t_s", cfg.dsmc.delta_t_s),
            ("dsmc.end_time_s", cfg.dsmc.end_time_s),
            ("dsmc.write_interval_s", cfg.dsmc.write_interval_s)):
        if value is not None and not (math.isfinite(float(value))
                                      and float(value) > 0.0):
            fail(f"{name} is {value}; it must be a positive number"
                 + ("" if name.startswith("exit.") else
                    ", or null to derive it"))
    if cfg.dsmc.average_start_s is not None and not float(
            cfg.dsmc.average_start_s) >= 0.0:
        fail(f"dsmc.average_start_s is {cfg.dsmc.average_start_s}; use null to "
             f"derive it, or a time >= 0")

    # --- the gas -------------------------------------------------------------
    if not SPECIES_NAME.match(str(cfg.gas.species_name)):
        fail(f"gas.species_name is {cfg.gas.species_name!r}. It becomes the "
             f"dsmcFoam typeId and part of the file name "
             f"0/boundaryNumberDensity_<name>, so use letters, digits and "
             f"underscores, starting with a letter.")
    if not 0.5 <= float(cfg.gas.omega) <= 1.0:
        fail(f"gas.omega is {cfg.gas.omega}; the VHS exponent lies in [0.5, 1.0] "
             f"(0.5 is a hard sphere)")
    dof = cfg.gas.internal_degrees_of_freedom
    if not _is_whole(dof) or int(dof) < 0:
        fail(f"gas.internal_degrees_of_freedom is {dof}; it must be a whole "
             f"number >= 0")
    if cfg.gas.mean_free_path_convention not in MEAN_FREE_PATH_CONVENTIONS:
        fail(f"unknown gas.mean_free_path_convention "
             f"{cfg.gas.mean_free_path_convention!r}; valid: "
             f"{list(MEAN_FREE_PATH_CONVENTIONS)}")

    # --- enumerations --------------------------------------------------------
    if cfg.dsmc.collision_model not in COLLISION_MODELS:
        fail(f"unknown dsmc.collision_model {cfg.dsmc.collision_model!r}; "
             f"valid: {list(COLLISION_MODELS)}")
    if cfg.dsmc.wall_interaction_model not in WALL_MODELS:
        fail(f"unknown dsmc.wall_interaction_model "
             f"{cfg.dsmc.wall_interaction_model!r}; valid: {list(WALL_MODELS)}")
    if cfg.dsmc.decomposition_method not in DECOMPOSITION_METHODS:
        fail(f"unknown dsmc.decomposition_method "
             f"{cfg.dsmc.decomposition_method!r}; valid: "
             f"{list(DECOMPOSITION_METHODS)}")
    if not 0.0 <= float(cfg.dsmc.diffuse_fraction) <= 1.0:
        fail(f"dsmc.diffuse_fraction is {cfg.dsmc.diffuse_fraction}; it is a "
             f"fraction, so it lies in [0, 1]")
    if not _is_whole(cfg.dsmc.n_subdomains) or int(cfg.dsmc.n_subdomains) < 1:
        fail(f"dsmc.n_subdomains is {cfg.dsmc.n_subdomains}; it must be a whole "
             f"number >= 1")

    # --- the domain ----------------------------------------------------------
    planes = tuple(cfg.domain.symmetry_planes)
    unknown = sorted(set(planes) - set(SYMMETRY_AXES))
    if unknown:
        fail(f"domain.symmetry_planes has {unknown}; valid entries are "
             f"{list(SYMMETRY_AXES)}. 'y' is the plane y = 0. The plume axis x "
             f"is never a symmetry plane.")
    if len(set(planes)) != len(planes):
        fail(f"domain.symmetry_planes repeats an axis: {list(planes)}")
    if cfg.plate is not None:
        for axis in planes:
            centre = getattr(cfg.plate, f"centre_{axis}_m")
            if float(centre) != 0.0:
                fail(f"domain.symmetry_planes includes '{axis}', but "
                     f"plate.centre_{axis}_m is {centre}. A plate off the plume "
                     f"axis is not symmetric about {axis} = 0; centre it, or "
                     f"drop that symmetry plane.")

    # --- the mesh ------------------------------------------------------------
    m = cfg.mesh
    for name, value, floor in (
            ("mesh.plate_refinement_level", m.plate_refinement_level, 0),
            ("mesh.n_cells_between_levels", m.n_cells_between_levels, 1),
            ("mesh.max_cells", m.max_cells, 1)):
        if not _is_whole(value) or int(value) < floor:
            fail(f"{name} is {value}; it must be a whole number >= {floor}")

    names = [m.patch_names[role] for role in DEFAULT_PATCH_NAMES]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        fail(f"mesh.patch_names reuses {duplicates} for more than one role; two "
             f"roles sharing a patch name would merge two boundaries with "
             f"different physics")

    region_names = [r.name for r in m.refinement_regions]
    duplicates = sorted({n for n in region_names if region_names.count(n) > 1})
    if duplicates:
        fail(f"mesh.refinement_regions repeats the name(s) {duplicates}")
    for index, region in enumerate(m.refinement_regions):
        where = f"mesh.refinement_regions[{index}] ({region.name!r})"
        if not _REGION_NAME.match(str(region.name)):
            fail(f"{where}: the name becomes a snappyHexMeshDict entry, so use "
                 f"letters, digits and underscores, starting with a letter")
        if region.name in names:
            fail(f"{where}: the name is also a patch name; snappyHexMesh keeps "
                 f"both in one geometry dictionary, so one would replace the "
                 f"other")
        if region.shape not in REGION_SHAPES:
            fail(f"{where}: unknown shape {region.shape!r}; valid: "
                 f"{list(REGION_SHAPES)}")
        if not _is_whole(region.level) or int(region.level) < 1:
            fail(f"{where}: level is {region.level}; it must be a whole number "
                 f">= 1 (level 0 is the background, which needs no region)")
        if not float(region.x_max_m) > float(region.x_min_m):
            fail(f"{where}: x_max_m ({region.x_max_m}) must be greater than "
                 f"x_min_m ({region.x_min_m})")
        if region.shape == "cone":
            stray = [k for k in _BOX_KEYS if getattr(region, k) is not None]
            if stray:
                fail(f"{where}: {stray} belong to shape 'box', and this is a "
                     f"cone")
            missing = [k for k in ("radius1_m", "radius2_m")
                       if getattr(region, k) is None]
            if missing:
                fail(f"{where}: a cone needs {missing} -- radius1_m is its "
                     f"radius at x_min_m, radius2_m its radius at x_max_m. "
                     f"Equal, they make a cylinder.")
            thickness = 0.0
            for end in ("1", "2"):
                outer = float(getattr(region, f"radius{end}_m"))
                inner = getattr(region, f"inner_radius{end}_m")
                inner = 0.0 if inner is None else float(inner)
                if not (math.isfinite(outer) and outer >= 0.0):
                    fail(f"{where}: radius{end}_m is {outer}; a radius is a "
                         f"number >= 0 (zero puts the cone's apex at that end)")
                if not (math.isfinite(inner) and 0.0 <= inner <= outer):
                    fail(f"{where}: inner_radius{end}_m is {inner}; it lies "
                         f"between 0 (solid) and radius{end}_m ({outer})")
                thickness = max(thickness, outer - inner)
            if not thickness > 0.0:
                fail(f"{where}: the cone encloses no volume: at neither end is "
                     f"its radius greater than its inner radius")
        else:
            missing = [k for k in _BOX_KEYS if getattr(region, k) is None]
            if missing:
                fail(f"{where}: a box needs {missing}")
            stray = [k for k in _CONE_KEYS if getattr(region, k) is not None]
            if stray:
                fail(f"{where}: {stray} belong to shape 'cone', and this is a "
                     f"box")
            for axis in ("y", "z"):
                lo = getattr(region, f"{axis}_min_m")
                hi = getattr(region, f"{axis}_max_m")
                if not float(hi) > float(lo):
                    fail(f"{where}: {axis}_max_m ({hi}) must be greater than "
                         f"{axis}_min_m ({lo})")

    # --- diagnostics ---------------------------------------------------------
    if "D" in cfg.diagnostics.reference_lengths_m:
        fail("diagnostics.reference_lengths_m has 'D', which is reserved: it is "
             "always the nozzle exit diameter, nozzle.exit_diameter_m, and gives "
             "KnD. Name the other length something else.")

    # --- checks and post -----------------------------------------------------
    if not 0.0 < float(cfg.checks.warn_courant) <= float(cfg.checks.max_courant):
        fail(f"checks.warn_courant ({cfg.checks.warn_courant}) must be positive "
             f"and no greater than checks.max_courant ({cfg.checks.max_courant})")
    if int(cfg.post.centreline_points) < 2:
        fail("post.centreline_points must be >= 2")
