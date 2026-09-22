r"""``case.yaml`` for the Cai 2012 family: schema, loading, validation.

Why this is not :mod:`plumetools.config`
----------------------------------------
That module's schema is welded to the source-flow model. It requires a
``stagnation`` section with a throat radius, cross-checks
``mesh.sphere_radius_m == geometry.sphere_radius_m``, and validates a
hemispherical inflow patch. None of those exist in this case: the inlet is a
circular nozzle disk with no sphere, no throat and no reservoir. Bolting an
optional Cai section onto it would mean every source-flow case carries keys it
must never set, and every Cai case carries a sphere radius it must never read --
and the frozen regression golden pinned to that module is exactly the thing §8 of
the case specification says not to disturb.

So this family owns its schema. The two loaders share a *shape* -- dataclasses,
unknown keys rejected, units in every key name -- but no code, and they cannot
load each other's files: ``model:`` must be ``cai2012_circular_plume`` here and
is rejected otherwise.

Provenance
----------
Every key is tagged in ``cases/cai2012/baseCase/case.yaml`` with ``[PAPER]``,
``[DERIVED]`` or ``[ASSUMPTION]``. A ``null`` on a numeric key means **derived at
mesh time** -- the cell size, the particle weight and the time step are all
computed from the Knudsen number rather than typed in, and a case that sets them
explicitly is overriding a derivation, which the manifest records.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

import yaml

from plumetools.cai2012 import MODEL
from plumetools.cai2012.gas import MEAN_FREE_PATH_CONVENTIONS, VhsSpecies


class CaiConfigError(ValueError):
    """``case.yaml`` is missing, malformed, or internally inconsistent."""


#: Where the Knudsen number's characteristic length comes from.
#:
#: Cai states the nozzle **diameter**. ``radius`` exists so that a comparison
#: with a paper using the other convention is a config change rather than a
#: silent factor of two in every density.
CHARACTERISTIC_LENGTHS = ("diameter", "radius")

#: Symmetry reductions of the full 3-D domain.
#:
#: ``none`` is the requested baseline: a full 3-D Cartesian box. ``half_y``
#: models ``y >= 0`` with a symmetry plane at ``y = 0``, halving the cell count
#: and the particle count. It is **opt-in**: nothing switches to it
#: automatically, because a domain that quietly became half a domain is not a
#: result anybody can compare.
SYMMETRY_MODES = ("none", "half_y")

#: What sets the start of sampling.
#:
#: ``cai``      ``dsmc.transient_collision_times`` reference collision times,
#:              Cai's own statement of 10 000 t0.
#: ``transits`` ``dsmc.transient_domain_transits`` transits of the domain at the
#:              exit speed. Cheaper, and the honest choice when Cai's transient
#:              is longer than the machine can run -- the manifest records which
#:              was used and both durations are always printed.
TRANSIENT_BASES = ("cai", "transits")

#: Inflow boundary models this family can write coefficients for.
INFLOW_MODELS = ("plumeFieldInflow",)


@dataclass(frozen=True)
class NozzleConfig:
    """The circular exit.

    Attributes:
        diameter_m: ``D`` [m]. **[PAPER]** 0.2 m.
    """

    diameter_m: float = 0.2

    @property
    def radius_m(self) -> float:
        """``R0 = D/2`` [m]. [DERIVED] -- there is deliberately no radius key."""
        return 0.5 * float(self.diameter_m)

    @property
    def area_m2(self) -> float:
        """``pi R0^2`` [m^2], the exact disk area the meshed patch approximates."""
        return math.pi * self.radius_m ** 2


@dataclass(frozen=True)
class ExitConfig:
    """The uniform drifting Maxwellian on the exit disk.

    Attributes:
        speed_ratio: ``S0 = U0 / sqrt(2 R T0)``. **[PAPER]** 2.0.
        T0_K: exit temperature [K]. **[ASSUMPTION]** -- Cai does not appear to
            print an absolute ``T0``, and the normalised quantities he validates
            depend on ``S0`` and ``Kn``, not on it. 300 K, stated rather than
            invented as a paper value.
        knudsen: ``Kn`` for this case. **[PAPER]** one of 100, 0.1, 0.01;
            ``generate_cases.py`` writes it per case.
        characteristic_length: which nozzle dimension ``Kn`` is built on.
            **[PAPER]** the diameter.
    """

    speed_ratio: float = 2.0
    T0_K: float = 300.0
    knudsen: float = 100.0
    characteristic_length: str = "diameter"


@dataclass(frozen=True)
class GasConfig:
    """The VHS species. See :data:`plumetools.cai2012.gas.ARGON`.

    Attributes:
        species_name: the ``typeId`` in ``dsmcProperties``. **[PAPER]** argon.
        mass_kg, diameter_m, omega, t_ref_K: VHS constants.
            **[ASSUMPTION / existing repository model]** -- the paper names the
            gas and the collision model but not the coefficients.
        internal_degrees_of_freedom: 0; argon is monatomic.
        mean_free_path_convention: ``vhs`` or ``hard_sphere``.
            **[ASSUMPTION]** -- Cai does not state GRASP's convention.
    """

    species_name: str = "Ar"
    mass_kg: float = 6.63e-26
    diameter_m: float = 4.17e-10
    omega: float = 0.74
    t_ref_K: float = 273.0
    internal_degrees_of_freedom: int = 0
    mean_free_path_convention: str = "vhs"

    def species(self) -> VhsSpecies:
        """This section as a :class:`~plumetools.cai2012.gas.VhsSpecies`."""
        return VhsSpecies(
            name=self.species_name,
            mass_kg=float(self.mass_kg),
            diameter_m=float(self.diameter_m),
            omega=float(self.omega),
            t_ref_K=float(self.t_ref_K),
            internal_degrees_of_freedom=int(self.internal_degrees_of_freedom),
        )


@dataclass(frozen=True)
class GeometryConfig:
    """The physical domain, in nozzle diameters.

    Attributes:
        symmetry_mode: one of :data:`SYMMETRY_MODES`.
        x_max_over_D: downstream extent. **[PAPER]** the sampling range Cai
            plots is ``0 <= X/D <= 10``.
        y_half_over_D, z_half_over_D: lateral half-extents.
        plume_axis: recorded, not configurable in practice -- the mesh, the
            inflow and the post-processing all assume ``+x``, and a case that
            set anything else would generate a mesh that disagreed with them.
    """

    symmetry_mode: str = "none"
    x_max_over_D: float = 10.0
    y_half_over_D: float = 10.0
    z_half_over_D: float = 10.0
    plume_axis: tuple = (1.0, 0.0, 0.0)


@dataclass(frozen=True)
class MeshConfig:
    """The graded multi-block Cartesian mesh.

    A fine uniform **core** around the nozzle and the near plume, with geometric
    expansion outward in every direction. Cai's uniform axisymmetric grid does
    not translate: at his ``dx = lambda0(Kn = 0.01) = 2 mm`` the full 3-D box is
    4e9 cells, which :func:`plumetools.cai2012.mesh.uniform_cost` computes and
    prints before anything is written.

    Attributes:
        type: ``graded_cartesian``. The only generator here.
        core_cell_size_m: cell size in the core [m]. ``null`` -> derived as
            ``min(target_cell_over_mfp * lambda0, max_core_cell_over_D * D)``,
            then coarsened if the cell budget requires it (which is reported,
            never silent).
        target_cell_over_mfp: the DSMC criterion the derivation aims at.
            **[PAPER]** Cai states ``dx / lambda0 = 1``.
        max_core_cell_over_D: cap on the derived core cell, in nozzle diameters.
            The collision criterion alone is useless at ``Kn = 100``, where
            ``lambda0`` is a hundred nozzle diameters and would allow a cell
            larger than the nozzle. This is the **geometric** requirement that
            sits underneath it: the exit disk and the near plume have to be
            resolved whatever the collision rate is. 0.05 gives 20 cells across
            the diameter.
        core_x_over_D, core_half_over_D: extent of the uniform core.
        outer_x_cells, outer_lateral_cells: cells in each expanding segment.
        outer_expansion: ``simpleGrading`` ratio across an outer segment, last
            cell over first. Applied outward from the core in every direction.
        max_cells: refuse to write a mesh larger than this. The guard §5 asks
            for: a mis-set Knudsen number would otherwise silently request a
            multi-billion-cell mesh. It binds on ``core_cell_size_m`` too, not
            only on the derived cell -- a pinned value that does not fit is
            coarsened and the fact is reported, because the alternative is a
            mesh nobody asked for on a machine that cannot hold it.
        min_core_cell_size_m: floor on the derived core cell [m], or ``null``.
        patch_names: roles ``nozzle`` / ``upstream_vacuum`` / ``outer`` /
            ``symmetry`` to mesh patch names.
    """

    type: str = "graded_cartesian"
    core_cell_size_m: float | None = None
    target_cell_over_mfp: float = 1.0
    max_core_cell_over_D: float = 0.05
    core_x_over_D: float = 2.0
    core_half_over_D: float = 1.5
    outer_x_cells: int = 40
    outer_lateral_cells: int = 30
    outer_expansion: float = 40.0
    max_cells: int = 4_000_000
    min_core_cell_size_m: float | None = None
    patch_names: dict = field(default_factory=lambda: {
        "nozzle": "nozzle",
        "upstream_vacuum": "upstreamVacuum",
        "outer": "vacuum",
        "symmetry": "symmetry",
    })


@dataclass(frozen=True)
class DsmcConfig:
    """Solver settings written into ``dsmcProperties`` and ``controlDict``.

    Attributes:
        binary_collision_model: **[PAPER]** VHS. ``VariableHardSphere`` is
            standard ``dsmcFoam``'s name for it. NTC is not a setting: it is the
            only collision *selection* algorithm ``DSMCCloud`` implements, so
            Cai's "no-time-counter" is satisfied by construction.
        collisions_enabled: ``true`` for **every** case including ``Kn = 100``.
            §3 of the case specification is explicit: the point is that DSMC
            approaches the collisionless solution on its own, not that collisions
            were switched off to make it.
        n_equivalent_particles: particle weight. ``null`` -> derived from
            ``resolution``.
        numerical_particle_multiplier: how many times the baseline numerical
            **parcel** population this case carries. One parcel stands for
            ``nEquivalentParticles`` real molecules, so ``N_numerical`` goes as
            ``1 / nEquivalentParticles`` and the weight is *divided* by this:

            .. code-block:: text

                nEquivalentParticles = baseline / numerical_particle_multiplier

            ``1.0`` is the baseline and leaves the weight exactly as derived (or
            exactly as pinned). It changes no physics -- the density, the mesh
            and the time step are untouched -- only how finely the same physical
            gas is sampled, which is what a statistical-resolution sweep varies.
        run_time_multiplier: scale on ``end_time_s``, applied **after** the
            transient is fixed, so every extra second is sampling time and
            ``average_start_s`` does not move. ``1.0`` is the baseline.
        output_frequency_multiplier: how many times more often to write than the
            baseline schedule. ``writeControl`` is ``timeStep``, so the interval
            is an integer number of steps and is rounded to the nearest one;
            :class:`~plumetools.cai2012.inflow.RunSettings` records the ratio
            actually achieved. ``1.0`` is the baseline.
        delta_t_s: time step [s]. ``null`` -> derived from ``courant_target``.
            **Not** touched by any of the multipliers above: a particle-count
            sweep that also moved the time step would confound the two.
        courant_target: fraction of the smallest cell a fast molecule may cross
            per step, used for the derivation.
        end_time_s, average_start_s, write_interval_s: ``null`` -> derived from
            the transient and sampling durations below.
        transient_basis: ``cai`` or ``transits``; see :data:`TRANSIENT_BASES`.
        transient_collision_times: **[PAPER]** Cai runs at least 10 000
            reference collision times before sampling.
        transient_domain_transits, sampling_domain_transits: the alternative
            basis, in transits of the domain length at the exit speed.
        initial_number_density_per_m3, initial_temperature_K: the near-vacuum
            fill ``dsmcInitialise`` writes. Not a physical state.
        n_subdomains: MPI ranks for ``decomposeParDict``.
    """

    binary_collision_model: str = "VariableHardSphere"
    collisions_enabled: bool = True
    n_equivalent_particles: float | None = None
    numerical_particle_multiplier: float = 1.0
    run_time_multiplier: float = 1.0
    output_frequency_multiplier: float = 1.0
    delta_t_s: float | None = None
    courant_target: float = 0.2
    end_time_s: float | None = None
    average_start_s: float | None = None
    write_interval_s: float | None = None
    transient_basis: str = "cai"
    transient_collision_times: float = 10000.0
    transient_domain_transits: float = 3.0
    sampling_domain_transits: float = 3.0
    initial_number_density_per_m3: float = 1.0e+10
    initial_temperature_K: float = 300.0
    n_subdomains: int = 4


@dataclass(frozen=True)
class ResolutionConfig:
    """Particle weighting.

    Attributes:
        target_particles_per_cell: occupancy target in the **exit** cell, which
            is the densest cell in the domain. Standard ``dsmcFoam`` has one
            weight for the whole domain, so this is the only place the target can
            be met by construction; everywhere else follows from the density and
            cell-size ratios and is reported, not assumed.
        reference_knudsen: **[PAPER]** Cai's ``dx/lambda0`` and ``dt/t0`` are
            referred to the ``Kn = 0.01`` exit properties, whatever case is being
            run. Used only for reporting the reference lengths and times.
        report_only: keep an explicit ``n_equivalent_particles`` instead of
            re-deriving one.
    """

    target_particles_per_cell: float = 20.0
    reference_knudsen: float = 0.01
    report_only: bool = False


@dataclass(frozen=True)
class ChecksConfig:
    """DSMC quality thresholds. See :mod:`plumetools.cai2012.checks`."""

    enabled: bool = True
    max_cell_over_mfp: float = 1.0
    max_courant: float = 1.0
    warn_courant: float = 0.5
    min_particles_per_cell: float = 5.0
    min_nozzle_faces: int = 12
    max_nozzle_area_error: float = 0.05


@dataclass(frozen=True)
class OutputConfig:
    """How ``0/`` and ``constant/`` are written.

    Attributes:
        dialect: ``standard`` -- OpenFOAM's own ``dsmcFoam``, v2512. There is no
            ``mnf`` path in this family: the fork reads ``boundaryT`` as a
            vector and this writes a scalar.
        inflow_model: ``plumeFieldInflow``. Stock ``FreeStream`` cannot express
            this case at all -- it injects on **every** ``patch``-type boundary,
            so the vacuum boundaries would become inlets, and the only way to
            stop that is to make them walls, which turns the vacuum into a box.
    """

    dialect: str = "standard"
    inflow_model: str = "plumeFieldInflow"


@dataclass(frozen=True)
class PostConfig:
    """Post-processing.

    Attributes:
        centerline_x_over_D_max: end of the extracted centreline. **[PAPER]**
            Cai plots to ``X/D = 10``.
        centerline_points: samples along it.
        centerline_radius_over_D: cells within this radius of the axis are
            averaged into a centreline sample. A single cell column is noise; a
            fat tube is not the centreline.
        plane: ``xz`` or ``xy`` -- which centre plane the contour is taken in.
            **[PAPER]** Cai's Figs. 6-8 are ``X/D`` against ``Z/D``.
        plane_half_over_D: lateral half-extent of the contour plane.
        contour_levels: **[PAPER]** the ``n/n0`` levels Cai draws.
        plane_points_x, plane_points_lateral: sampling resolution of the plane.
    """

    centerline_x_over_D_max: float = 10.0
    centerline_points: int = 200
    centerline_radius_over_D: float = 0.25
    plane: str = "xz"
    plane_half_over_D: float = 5.0
    contour_levels: tuple = (0.1, 0.01, 0.001)
    plane_points_x: int = 240
    plane_points_lateral: int = 240


@dataclass(frozen=True)
class CaiCaseConfig:
    """A whole validated ``case.yaml``."""

    model: str = MODEL
    nozzle: NozzleConfig = field(default_factory=NozzleConfig)
    exit: ExitConfig = field(default_factory=ExitConfig)
    gas: GasConfig = field(default_factory=GasConfig)
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    mesh: MeshConfig = field(default_factory=MeshConfig)
    dsmc: DsmcConfig = field(default_factory=DsmcConfig)
    resolution: ResolutionConfig = field(default_factory=ResolutionConfig)
    checks: ChecksConfig = field(default_factory=ChecksConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    post: PostConfig = field(default_factory=PostConfig)
    #: Free-form provenance. Not validated: it records where a case came from,
    #: and over-constraining that makes recording an awkward fact harder than
    #: omitting it.
    meta: dict = field(default_factory=dict)

    def species(self) -> VhsSpecies:
        """The VHS species, as the gas section describes it."""
        return self.gas.species()

    @property
    def characteristic_length_m(self) -> float:
        """The length ``Kn`` is built on [m]. [DERIVED] from ``nozzle``."""
        if self.exit.characteristic_length == "radius":
            return self.nozzle.radius_m
        return float(self.nozzle.diameter_m)


_SECTIONS = {
    "nozzle": NozzleConfig,
    "exit": ExitConfig,
    "gas": GasConfig,
    "geometry": GeometryConfig,
    "mesh": MeshConfig,
    "dsmc": DsmcConfig,
    "resolution": ResolutionConfig,
    "checks": ChecksConfig,
    "output": OutputConfig,
    "post": PostConfig,
}

_PASSTHROUGH = ("model", "meta")


def _build(cls, data: dict, where: str):
    """Instantiate a section dataclass, rejecting unknown keys.

    Unknown keys are an error rather than being ignored, so a typo in a physical
    parameter stops the run instead of silently falling back to a default.

    Numeric-looking strings are coerced, because YAML 1.1 parses ``1.0e14`` as a
    **string** (it requires a signed exponent) while ``1.0e+14`` is a float. The
    failure that causes without this is an unformattable value deep inside a
    dictionary renderer, a long way from the line responsible.
    """
    if not isinstance(data, dict):
        raise CaiConfigError(
            f"{where}: expected a mapping, got {type(data).__name__}")
    known = {f.name for f in fields(cls)}
    unknown = sorted(set(data) - known)
    if unknown:
        raise CaiConfigError(
            f"{where}: unknown key(s) {unknown}; valid keys are {sorted(known)}")

    coerced: dict[str, Any] = {}
    for f in fields(cls):
        if f.name not in data:
            continue
        value = data[f.name]
        annotation = str(f.type)
        if value is None:
            coerced[f.name] = None
            continue
        if annotation in ("tuple", "tuple | None") and isinstance(value, list):
            value = tuple(value)
        elif isinstance(value, str) and ("float" in annotation or "int" in annotation):
            try:
                value = float(value) if "float" in annotation else int(value)
            except ValueError:
                raise CaiConfigError(
                    f"{where}: {f.name} is {value!r}, which is not a number. If "
                    f"it looks like one, check the exponent: YAML requires a "
                    f"signed exponent, so 1.0e+14 is a float and 1.0e14 is a "
                    f"string.") from None
        coerced[f.name] = value
    return cls(**coerced)


def load_case_config(case_dir: Path) -> CaiCaseConfig:
    """Read and validate ``<case_dir>/case.yaml``.

    Args:
        case_dir: a case directory, or the ``case.yaml`` file itself.

    Returns:
        A validated :class:`CaiCaseConfig`.

    Raises:
        CaiConfigError: for a missing file, an unknown key, the wrong ``model``,
            or any failed cross-check in :func:`validate`.
    """
    path = Path(case_dir)
    if path.is_dir():
        path = path / "case.yaml"
    if not path.is_file():
        raise CaiConfigError(f"no case.yaml at {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise CaiConfigError(f"{path}: expected a mapping at the top level")

    unknown = sorted(set(raw) - set(_SECTIONS) - set(_PASSTHROUGH))
    if unknown:
        raise CaiConfigError(
            f"{path}: unknown top-level key(s) {unknown}; "
            f"valid: {sorted(set(_SECTIONS) | set(_PASSTHROUGH))}")

    model = raw.get("model")
    if model != MODEL:
        raise CaiConfigError(
            f"{path}: model is {model!r}, expected {MODEL!r}. This loader reads "
            f"the Cai 2012 circular-plume schema only; a source-flow case.yaml "
            f"belongs to plumetools.config, and the two families share no keys.")

    cfg = CaiCaseConfig(
        model=model,
        meta=raw.get("meta") or {},
        **{name: _build(cls, raw.get(name) or {}, f"{path}:{name}")
           for name, cls in _SECTIONS.items()},
    )
    validate(cfg, path)
    return cfg


#: Patch roles a ``graded_cartesian`` mesh must name.
REQUIRED_PATCH_ROLES = ("nozzle", "upstream_vacuum", "outer")


def validate(cfg: CaiCaseConfig, path: Path | str = "case.yaml") -> None:
    """Cross-check a config. Called by :func:`load_case_config`.

    Every check here would otherwise surface as a meshing failure, a solver abort
    with no indication of which key caused it, or -- worse -- a case that runs
    cleanly at the wrong physical conditions.

    Raises:
        CaiConfigError: on the first inconsistency, naming the key and the fix.
    """
    def fail(message: str):
        raise CaiConfigError(f"{path}: {message}")

    # --- positive physical inputs -------------------------------------------
    positives = {
        "nozzle.diameter_m": cfg.nozzle.diameter_m,
        "exit.T0_K": cfg.exit.T0_K,
        "exit.knudsen": cfg.exit.knudsen,
        "gas.mass_kg": cfg.gas.mass_kg,
        "gas.diameter_m": cfg.gas.diameter_m,
        "gas.t_ref_K": cfg.gas.t_ref_K,
        "geometry.x_max_over_D": cfg.geometry.x_max_over_D,
        "geometry.y_half_over_D": cfg.geometry.y_half_over_D,
        "geometry.z_half_over_D": cfg.geometry.z_half_over_D,
        "mesh.core_x_over_D": cfg.mesh.core_x_over_D,
        "mesh.core_half_over_D": cfg.mesh.core_half_over_D,
        "mesh.target_cell_over_mfp": cfg.mesh.target_cell_over_mfp,
        "mesh.max_core_cell_over_D": cfg.mesh.max_core_cell_over_D,
        "mesh.outer_expansion": cfg.mesh.outer_expansion,
        "resolution.target_particles_per_cell":
            cfg.resolution.target_particles_per_cell,
        "resolution.reference_knudsen": cfg.resolution.reference_knudsen,
        "dsmc.courant_target": cfg.dsmc.courant_target,
        "dsmc.numerical_particle_multiplier":
            cfg.dsmc.numerical_particle_multiplier,
        "dsmc.run_time_multiplier": cfg.dsmc.run_time_multiplier,
        "dsmc.output_frequency_multiplier": cfg.dsmc.output_frequency_multiplier,
        "dsmc.initial_number_density_per_m3": cfg.dsmc.initial_number_density_per_m3,
        "dsmc.initial_temperature_K": cfg.dsmc.initial_temperature_K,
    }
    bad = sorted(name for name, value in positives.items()
                 if value is None or not math.isfinite(float(value))
                 or float(value) <= 0.0)
    if bad:
        fail(f"{bad} must be finite and positive")

    for name, value in (("mesh.core_cell_size_m", cfg.mesh.core_cell_size_m),
                        ("mesh.min_core_cell_size_m", cfg.mesh.min_core_cell_size_m),
                        ("dsmc.delta_t_s", cfg.dsmc.delta_t_s),
                        ("dsmc.n_equivalent_particles",
                         cfg.dsmc.n_equivalent_particles),
                        ("dsmc.end_time_s", cfg.dsmc.end_time_s)):
        if value is not None and not float(value) > 0.0:
            fail(f"{name} is {value}; use null to derive it, or a positive number")

    # --- enumerations --------------------------------------------------------
    if cfg.exit.characteristic_length not in CHARACTERISTIC_LENGTHS:
        fail(f"unknown exit.characteristic_length "
             f"{cfg.exit.characteristic_length!r}; valid: "
             f"{list(CHARACTERISTIC_LENGTHS)}. Cai uses the diameter.")
    if cfg.gas.mean_free_path_convention not in MEAN_FREE_PATH_CONVENTIONS:
        fail(f"unknown gas.mean_free_path_convention "
             f"{cfg.gas.mean_free_path_convention!r}; valid: "
             f"{list(MEAN_FREE_PATH_CONVENTIONS)}")
    if cfg.geometry.symmetry_mode not in SYMMETRY_MODES:
        fail(f"unknown geometry.symmetry_mode {cfg.geometry.symmetry_mode!r}; "
             f"valid: {list(SYMMETRY_MODES)}")
    if cfg.dsmc.transient_basis not in TRANSIENT_BASES:
        fail(f"unknown dsmc.transient_basis {cfg.dsmc.transient_basis!r}; "
             f"valid: {list(TRANSIENT_BASES)}")
    if cfg.mesh.type != "graded_cartesian":
        fail(f"unknown mesh.type {cfg.mesh.type!r}; this family has one "
             f"generator, 'graded_cartesian'")
    if cfg.output.dialect != "standard":
        fail(f"output.dialect is {cfg.output.dialect!r}. This family targets "
             f"standard OpenFOAM dsmcFoam (v2512) only -- it writes boundaryT as "
             f"a volScalarField, which the MNF fork cannot read.")
    if cfg.output.inflow_model not in INFLOW_MODELS:
        fail(f"unknown output.inflow_model {cfg.output.inflow_model!r}; valid: "
             f"{list(INFLOW_MODELS)}. Stock FreeStream cannot express this case: "
             f"it injects on every patch-type boundary, so the vacuum boundaries "
             f"would become inlets.")
    if cfg.post.plane not in ("xz", "xy"):
        fail(f"unknown post.plane {cfg.post.plane!r}; valid: ['xz', 'xy']")

    # --- structural ----------------------------------------------------------
    missing_roles = [r for r in REQUIRED_PATCH_ROLES
                     if r not in cfg.mesh.patch_names]
    if missing_roles:
        fail(f"mesh.patch_names is missing role(s) {missing_roles}; "
             f"'graded_cartesian' needs all of {list(REQUIRED_PATCH_ROLES)}")
    if cfg.geometry.symmetry_mode == "half_y" and "symmetry" not in cfg.mesh.patch_names:
        fail("geometry.symmetry_mode is 'half_y' but mesh.patch_names has no "
             "'symmetry' role; the y = 0 plane would be written as an open "
             "boundary and would delete every particle that reached it")
    names = [cfg.mesh.patch_names[r] for r in REQUIRED_PATCH_ROLES]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        fail(f"mesh.patch_names reuses {duplicates} for more than one role; two "
             f"roles sharing a patch name would merge two boundaries with "
             f"different physics -- and here that means the vacuum boundary "
             f"would inject")

    if int(cfg.gas.internal_degrees_of_freedom) != 0:
        fail(f"gas.internal_degrees_of_freedom is "
             f"{cfg.gas.internal_degrees_of_freedom}; argon is monatomic. A "
             f"non-zero value would need a Larsen-Borgnakke collision model and "
             f"an internal-energy field this case does not write.")

    if cfg.dsmc.binary_collision_model != "VariableHardSphere":
        fail(f"dsmc.binary_collision_model is "
             f"{cfg.dsmc.binary_collision_model!r}; Cai states VHS, and a "
             f"monatomic gas has no internal modes for Larsen-Borgnakke to "
             f"redistribute")

    # --- geometric consistency ----------------------------------------------
    if cfg.mesh.max_core_cell_over_D > 0.25:
        fail(f"mesh.max_core_cell_over_D is {cfg.mesh.max_core_cell_over_D}; "
             f"that allows a core cell of {cfg.mesh.max_core_cell_over_D:g} D, "
             f"so the exit disk would be fewer than four cells across and the "
             f"'circular' nozzle would be a square")
    if cfg.mesh.core_half_over_D <= 0.5:
        fail(f"mesh.core_half_over_D is {cfg.mesh.core_half_over_D}; the core "
             f"must be wider than the nozzle radius (0.5 D) or the exit disk "
             f"would straddle the grading boundary")
    if cfg.mesh.core_x_over_D >= cfg.geometry.x_max_over_D:
        fail(f"mesh.core_x_over_D ({cfg.mesh.core_x_over_D}) must be less than "
             f"geometry.x_max_over_D ({cfg.geometry.x_max_over_D})")
    for axis in ("y", "z"):
        half = getattr(cfg.geometry, f"{axis}_half_over_D")
        if cfg.mesh.core_half_over_D >= half:
            fail(f"mesh.core_half_over_D ({cfg.mesh.core_half_over_D}) must be "
                 f"less than geometry.{axis}_half_over_D ({half})")
    if cfg.post.centerline_x_over_D_max > cfg.geometry.x_max_over_D:
        fail(f"post.centerline_x_over_D_max "
             f"({cfg.post.centerline_x_over_D_max}) reaches past the domain "
             f"({cfg.geometry.x_max_over_D}); the extracted centreline would run "
             f"out of mesh")
    if cfg.mesh.outer_expansion < 1.0:
        fail(f"mesh.outer_expansion is {cfg.mesh.outer_expansion}; it is the "
             f"ratio of the outermost cell to the one next to the core, so it "
             f"must be >= 1 -- a value below 1 would refine towards the vacuum "
             f"boundary")
    if min(int(cfg.mesh.outer_x_cells), int(cfg.mesh.outer_lateral_cells)) < 1:
        fail("mesh.outer_x_cells and mesh.outer_lateral_cells must be >= 1")
    if int(cfg.mesh.max_cells) < 1:
        fail("mesh.max_cells must be >= 1")
    if int(cfg.dsmc.n_subdomains) < 1:
        fail("dsmc.n_subdomains must be >= 1")

    if not 0.0 < float(cfg.checks.warn_courant) <= float(cfg.checks.max_courant):
        fail(f"checks.warn_courant ({cfg.checks.warn_courant}) must be positive "
             f"and no greater than checks.max_courant ({cfg.checks.max_courant})")
