"""Load and validate ``case.yaml``.

``case.yaml`` carries only what OpenFOAM dictionaries cannot express: the
analytical inflow model's physical inputs and the mesh parameters. It never
duplicates or overrides ``controlDict`` or ``dsmcProperties``, which remain the
authority for everything the solver reads.

Every physical key carries its unit in its name (``p0_pa``, ``T0_K``,
``throat_radius_m``). The original code carried none, and its comments actively
misdescribed the values -- ``# 5psi in Pa`` was copy-pasted unchanged onto 34500,
6894 and 16894 Pa, which are 5.00, 1.00 and 2.45 psi (finding DOC-02).

Unknown keys are rejected rather than ignored, so a typo in a physical parameter
fails loudly instead of silently falling back to a default.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

import yaml

from plumetools.mesh.boundary import read_boundary


class ConfigError(ValueError):
    """Raised when ``case.yaml`` is missing, malformed, or inconsistent."""


class ConfigWarning(UserWarning):
    """Emitted for a frozen inconsistency that is preserved rather than fixed."""


#: Mesh generators, and which module renders each.
#:
#: ``snappy_markelov`` is the AIAA 99-3455 geometry -- a background box carved by
#: an inflow hemisphere, a finite cylinder and a plate -- rendered by
#: :mod:`plumetools.markelov1999.mesh`.
MESH_TYPES = ("block_mesh_ogrid", "snappy_hex_sphere", "snappy_markelov")

#: How the inflow surface is made curved, for ``block_mesh_ogrid``.
PROJECTIONS = ("arc", "searchable_sphere", "none")

#: Inflow boundary models a case can ask for.
#:
#: ``FreeStream`` is stock standard dsmcFoam: one number density per species, and
#: injection on every ``patch``-type boundary with no selection list.
#: ``plumeFieldInflow`` is the library under ``applications/dsmcBoundaryModels``,
#: which takes an explicit patch list and a per-face number density field.
#: ``dsmcFreeStreamInflowFieldPatch`` is the MNF fork's own equivalent.
INFLOW_MODELS = ("FreeStream", "plumeFieldInflow", "dsmcFreeStreamInflowFieldPatch")


@dataclass(frozen=True)
class MeshConfig:
    """Parameters for generating the mesh.

    Attributes:
        type: ``"snappy_hex_sphere"`` -- a plain background box from blockMesh,
            with the inflow cavity carved out of it by snappyHexMesh, rendered by
            :mod:`plumetools.foamio.snappy`; or ``"block_mesh_ogrid"`` -- a
            5-block hexahedral O-grid with the hemisphere built directly into
            ``blockMeshDict``, rendered by :mod:`plumetools.foamio.blockmesh`.
            The O-grid is exact and purely hexahedral but has a non-orthogonality
            floor set by its block topology; see docs/mesh.md.
        sphere_radius_m: inflow surface radius [m]. Must equal
            ``geometry.sphere_radius_m``.
        box_half_width_m: domain half-extent in y and z [m].
        box_length_m: domain extent in x [m], from 0.
        n_tangential: O-grid only. Cells across each block face; the inflow patch
            gets ``5 * n_tangential**2`` faces.
        n_radial: O-grid only. Cells from the sphere to the box.
        radial_grading: O-grid only. ``simpleGrading`` expansion outward.
        projection: O-grid only. ``"arc"`` curves only the twelve block edges, so
            the face interiors stay ruled surfaces and the patch is **not** a true
            hemisphere -- the deficit reaches 16.3% of R at the cap face centre and
            does **not** shrink with ``n_tangential``. ``"searchable_sphere"``
            projects the edges *and* faces onto a ``searchableSphere`` primitive,
            giving an exact surface at no extra cost. ``"none"`` is an accepted
            alias for ``"arc"``, kept because earlier configs used it.
        background_cell_size_m: snappy only. Background box cell size [m]. The
            background has to *find* the sphere, not resolve it -- surface
            resolution is ``background_cell_size_m / 2**refinement_level`` -- so
            this may be as large as ``R`` (two cells across the diameter). Above
            that the cavity falls between cells and the generator refuses.
        refinement_level: snappy only. Octree levels at the sphere. Raising it
            halves the surface cell size, so check it against ``controlDict``'s
            ``deltaT``: a particle must not cross a full cell in one step.
        n_cells_between_levels: snappy only. Buffer cells between refinement levels.
        outer_patch_type: geometric type for the outer boundary. ``"patch"`` is
            physically right -- particles leave a plume domain there.

            But standard dsmcFoam's FreeStream injects on **every**
            ``isType<polyPatch>`` boundary, with no selection list
            (``FreeStream.C:57-62``), so an outer ``patch`` also becomes an
            inflow. It aborts immediately with "Zero boundary temperature
            detected", because boundaryT is zeroGradient there and evaluates to
            the zero internal field.

            Setting ``"wall"`` excludes it -- ``isType<>`` is an exact match --
            and lets the case run. The cost is that the outer boundary then
            REFLECTS instead of absorbing: the plume expands into a closed box.
            That is a real physics change, so it is opt-in, not the default.
            See docs/solver-compatibility.md.
        patch_names: maps the roles ``inflow`` / ``outer`` / ``symmetry`` to the
            patch names written into the mesh.
    """

    type: str = "block_mesh_ogrid"
    sphere_radius_m: float = 0.5
    box_half_width_m: float = 2.5
    box_length_m: float = 5.0
    n_tangential: int = 20
    n_radial: int = 24
    radial_grading: float = 10.0
    projection: str = "arc"
    outer_patch_type: str = "patch"
    background_cell_size_m: float = 0.125
    refinement_level: int = 3
    n_cells_between_levels: int = 3
    patch_names: dict = field(
        default_factory=lambda: {"inflow": "inflow", "outer": "vacuum", "symmetry": "sym"}
    )
    # --- snappy_markelov only ------------------------------------------------
    # Six independent extents rather than a half-width and a length: the AIAA
    # 99-3455 domain is not a cube and is not symmetric about y, so no two of
    # these can be inferred from the others.
    x_min_m: float | None = None
    x_max_m: float | None = None
    y_min_m: float | None = None
    y_max_m: float | None = None
    z_min_m: float | None = None
    z_max_m: float | None = None
    #: Octree levels at each carved surface. ``None`` falls back to
    #: ``refinement_level``, so a case may set one number for all three.
    inflow_refinement_level: int | None = None
    cylinder_refinement_level: int | None = None
    plate_refinement_level: int | None = None


@dataclass(frozen=True)
class GeometryConfig:
    patch: str = "inflow"
    sphere_radius_m: float = 0.5
    plume_axis: tuple = (1.0, 0.0, 0.0)
    polar_axis: tuple = (0.0, 0.0, 1.0)


@dataclass(frozen=True)
class GasConfig:
    gamma: float = 1.4
    molar_mass_g_per_mol: float = 28.0134
    species_name: str = "Ar"


@dataclass(frozen=True)
class StagnationConfig:
    p0_pa: float = 3275009.71275  # 475 * PSI_TO_PA, to full precision
    T0_K: float = 300.0
    throat_radius_m: float = 0.0041275


@dataclass(frozen=True)
class AngularConfig:
    exponent_offset: float = 0.41
    form: str = "legacy_separable"
    normalization_gamma: float = 1.4
    quadrature_points: int = 500


#: Which dsmcFoam the case targets. The two read different boundary field types.
DIALECTS = ("standard", "mnf")


@dataclass(frozen=True)
class OutputConfig:
    """How the ``0/`` inflow fields are written.

    Attributes:
        inflow_model: which ``InflowBoundaryModel`` the case runs; one of
            :data:`INFLOW_MODELS`. This is not cosmetic. Stock ``FreeStream``
            injects on **every** ``patch``-type boundary with no selection list,
            which is what forces the ``outer_patch_type: wall`` compromise;
            ``plumeFieldInflow`` takes an explicit patch list, so the vacuum
            boundary can stay an open ``patch`` and actually absorb particles.
            The value is written into ``constant/dsmcProperties`` and drives the
            ``patch``-vs-``wall`` warning below.
        dialect: ``"standard"`` for OpenFOAM's own ``dsmcFoam`` (v2512 checked),
            ``"mnf"`` for the micro/nano-flow fork's ``dsmcFoam+``. They differ in
            one load-bearing way: standard reads ``0/boundaryT`` as a
            **volScalarField**, while the fork reads a **volVectorField** holding
            per-component translational temperature as ``(T 0 0)``. Writing the
            wrong one is a hard read failure, not a silent mis-run.

            Standard also has no per-face number density -- see
            :mod:`plumetools.foamio.fields` and docs/solver-compatibility.md.
        patches: non-inflow patch specs, in output order.
    """

    dialect: str = "standard"
    patches: dict = field(default_factory=dict)
    inflow_model: str = "FreeStream"


@dataclass(frozen=True)
class CylinderConfig:
    """The cylinder of AIAA 99-3455. Every value is printed by the paper.

    Attributes:
        radius_m: 3 in = 0.0762 m.
        length_m: extent along ``z``; 18 in = 0.4572 m.
        centre_x_m: 11.75 in = 0.29845 m.
        centre_z_m: 0 -- centred on the plume axis.
    """

    radius_m: float = 0.0762
    length_m: float = 0.4572
    centre_x_m: float = 0.29845
    centre_z_m: float = 0.0


@dataclass(frozen=True)
class PlateConfig:
    """The flat plate of AIAA 99-3455.

    Attributes:
        width_y_m: 6 in = 0.1524 m [PAPER].
        height_z_m: 15 in = 0.381 m [PAPER].
        thickness_m: extent in ``x``. **ASSUMPTION** -- the paper does not clearly
            establish it; baseline 0.5 in = 0.0127 m. See
            :mod:`plumetools.markelov1999.geometry`.

    There is deliberately no ``x`` here. The plate's position is *derived* from
    the cylinder and the gap, so a literal would be a second, silently divergent
    source of truth.
    """

    width_y_m: float = 0.1524
    height_z_m: float = 0.381
    thickness_m: float = 0.0127


@dataclass(frozen=True)
class BodiesConfig:
    """Solid bodies in the plume, and the symmetry plane.

    Attributes:
        gap_m: cylinder downstream surface to plate upstream face; 6 in = 0.1524 m
            [PAPER]. The plate's ``x`` follows from this.
        symmetry_plane_y_m: the modelled half is ``y >= this``. Only ``y`` is
            halved -- see :mod:`plumetools.markelov1999.geometry` for why ``z`` is
            not.
        cylinder: see :class:`CylinderConfig`.
        plate: see :class:`PlateConfig`.
    """

    gap_m: float = 0.1524
    symmetry_plane_y_m: float = 0.0
    cylinder: CylinderConfig = field(default_factory=CylinderConfig)
    plate: PlateConfig = field(default_factory=PlateConfig)

    _nested = {"cylinder": CylinderConfig, "plate": PlateConfig}


@dataclass(frozen=True)
class MoleculeConfig:
    """VHS / Larsen-Borgnakke molecular properties for one species.

    Attributes:
        name: the ``typeIdList`` entry, and the suffix of
            ``0/boundaryNumberDensity_<name>``.
        molar_mass_g_per_mol: 28.0134 for N2. **The authoritative mass input** --
            ``mass_kg`` is derived from it unless overridden.
        mass_kg: mass of one molecule. ``None`` derives it as
            ``M / (1000 * N_A)``, which is what the source-flow model uses, so the
            analytical model and the solver cannot disagree about ``m``. The
            OpenFOAM tutorials carry a rounded 46.5e-27 kg, 0.04% away; it is
            recorded in ``source`` but not used, because two masses that nearly
            agree are harder to debug than one.
        gamma: ratio of specific heats. 1.4 for a diatomic gas with 2 rotational
            degrees of freedom and no vibration -- consistent with
            ``internal_degrees_of_freedom`` below, and checked.
        diameter_m: VHS reference diameter.
        omega: VHS viscosity-temperature exponent.
        internal_degrees_of_freedom: 2 for N2 rotation. Vibration is not modelled
            by standard dsmcFoam's Larsen-Borgnakke implementation, so this is
            rotation only.
        source: where the non-paper coefficients came from. The paper prints the
            gas, the collision-model family and ZR, but not the VHS coefficients
            modern OpenFOAM needs, so these are a documented substitution and not
            a reproduction. Carried into every case summary.
    """

    name: str = "N2"
    molar_mass_g_per_mol: float = 28.0134
    mass_kg: float | None = None
    gamma: float = 1.4
    diameter_m: float = 4.17e-10
    omega: float = 0.74
    internal_degrees_of_freedom: int = 2
    source: str = ("OpenFOAM v2512 tutorials/discreteMethods/dsmcFoam/wedge15Ma5 "
                   "(diameter, omega, internalDegreesOfFreedom); mass derived from "
                   "the molar mass rather than the tutorial's rounded 46.5e-27 kg")


@dataclass(frozen=True)
class DsmcConfig:
    """Everything ``constant/dsmcProperties`` and ``system/controlDict`` need.

    Attributes:
        species: see :class:`MoleculeConfig`.
        n_equivalent_particles: the uniform particle weight. ``None`` means
            "derive from ``resolution``" -- :mod:`plumetools.markelov1999.resolution`
            picks it so the sizing region hits the target occupancy, and
            ``generate_cases.py`` writes the chosen number back into each
            case's ``case.yaml`` so it is auditable rather than implicit.

            Standard dsmcFoam has **one** weight for the whole domain
            (``DSMCCloud::nParticle_``, a single scalar). There is no radial or
            adaptive weighting to configure, so occupancy in the sparse regions
            follows from the choice made in the dense one.
        delta_t_s: solver time step.
        end_time_s: total simulated time.
        write_interval_s: how often fields are written.
        average_start_s: when ``fieldAverage`` starts accumulating. Everything
            before this is treated as transient and discarded.
        binary_collision_model: ``LarsenBorgnakkeVariableHardSphere`` gives VHS
            collisions plus internal-energy redistribution, which is the family
            the paper names.
        t_ref_K: VHS reference temperature.
        rotational_collision_number: ZR. The paper prints 5.
        wall_interaction_model: gas-surface model for the cylinder and plate.
        wall_temperature_K: written into ``0/boundaryT`` on the wall patches.
            ``MaxwellianThermal`` reads the wall temperature from there, not from
            its own coefficients, so a zero here would leave reflected particles
            with no thermal speed at all.
        initial_number_density_per_m3: what ``dsmcInitialise`` fills the domain
            with. A near-vacuum seed, not a physical state.
        initial_temperature_K: likewise.
        n_subdomains: default ``decomposePar`` count.
    """

    species: MoleculeConfig = field(default_factory=MoleculeConfig)
    n_equivalent_particles: float | None = None
    delta_t_s: float = 2.0e-7
    end_time_s: float = 4.0e-3
    write_interval_s: float = 5.0e-4
    average_start_s: float = 2.0e-3
    binary_collision_model: str = "LarsenBorgnakkeVariableHardSphere"
    t_ref_K: float = 273.0
    rotational_collision_number: float = 5.0
    wall_interaction_model: str = "MaxwellianThermal"
    wall_temperature_K: float = 300.0
    initial_number_density_per_m3: float = 1.0e14
    initial_temperature_K: float = 300.0
    n_subdomains: int = 4

    _nested = {"species": MoleculeConfig}


@dataclass(frozen=True)
class ResolutionConfig:
    """Particle-resolution targets and where they are measured.

    Attributes:
        target_particles_per_cell: the DSMC occupancy target. 20 for this study.
        sizing_region: which region of interest sets the uniform particle weight
            -- ``"cylinder"``, ``"plate"`` or ``"wake"``.

            With one weight for the whole domain the target can be met in exactly
            one region; the others follow from the density and cell-size ratios.
            Sizing on ``cylinder`` is the default because the windward/leeward
            pressure ratio is the primary result and the run stays affordable;
            sizing on ``plate`` meets the target everywhere at several times the
            particle count. Whichever is chosen, the estimator reports the
            occupancy in **all** regions, so the shortfall is never implicit.
        wake_offset_m: how far behind the cylinder base the wake sample sits.
        report_only: if true the estimator reports and never overrides an
            explicitly configured ``dsmc.n_equivalent_particles``.
    """

    target_particles_per_cell: float = 20.0
    sizing_region: str = "cylinder"
    wake_offset_m: float = 0.0254
    report_only: bool = False


@dataclass(frozen=True)
class ChecksConfig:
    """Thresholds for the DSMC mesh and time-step quality checks.

    Attributes:
        max_cell_over_mfp: warn when a cell exceeds this multiple of the local
            mean free path. DSMC theory wants a cell well under one mfp.
        max_courant: **hard error** above this. A particle crossing more than one
            of the smallest cells in a step invalidates the collision sampling,
            and the case specification lists it among the fail-fast conditions.
        warn_courant: warn above this.
        min_particles_per_cell: warn when a region of interest falls below this.
        enabled: run the checks at all. Off makes every check a no-op, which is
            occasionally useful when deliberately exploring a bad mesh.
    """

    max_cell_over_mfp: float = 1.0
    max_courant: float = 1.0
    warn_courant: float = 0.5
    min_particles_per_cell: float = 5.0
    enabled: bool = True


@dataclass(frozen=True)
class PressureWindowConfig:
    """One area-weighted surface-averaging window on the cylinder.

    Attributes:
        name: label used in the case summary and the study table.
        azimuth_deg: window centre, measured about the cylinder axis from ``+x``.
            ``180`` is windward (facing the source), ``0`` is leeward (the base).
        half_angle_deg: angular half-width of the window.
        axial_half_height_m: axial half-width about the cylinder mid-span.

    A window, never a single face: a one-face reading depends entirely on where
    snappyHexMesh happened to put that face, and would change with the refinement
    level rather than with the physics.
    """

    name: str = "windward"
    azimuth_deg: float = 180.0
    half_angle_deg: float = 15.0
    axial_half_height_m: float = 0.0381


@dataclass(frozen=True)
class PostConfig:
    """Post-processing configuration.

    Attributes:
        surface_pressure_field: the field whose wall-normal component is the
            surface pressure. ``fDMean`` is the time average of ``fD``, the force
            density dsmcFoam accumulates in ``DSMCParcel::hitWallPatch`` with
            dimensions ``[1 -1 -2 0 0 0 0]`` = Pa. A single-timestep ``fD`` is a
            one-step momentum tally and far too noisy to read as a pressure.
        windows: the averaging windows, in output order.
        cylinder_patch / plate_patch: patch names to measure on.
    """

    surface_pressure_field: str = "fDMean"
    cylinder_patch: str = "cylinder"
    plate_patch: str = "plate"
    windows: tuple = (
        {"name": "windward", "azimuth_deg": 180.0, "half_angle_deg": 15.0,
         "axial_half_height_m": 0.0381},
        {"name": "leeward", "azimuth_deg": 0.0, "half_angle_deg": 15.0,
         "axial_half_height_m": 0.0381},
        {"name": "midspan_side", "azimuth_deg": 90.0, "half_angle_deg": 15.0,
         "axial_half_height_m": 0.0381},
    )


@dataclass(frozen=True)
class SamplingConfig:
    radial_lines_deg: tuple = (0.0, 45.0, 90.0)
    r_start_m: float = 0.5
    r_end_m: float = 5.0
    n_points: int = 1000
    fields: tuple = ("U_Ar", "Ttra_Ar", "rhoN_Ar")


@dataclass(frozen=True)
class LegacyConfig:
    """Compatibility switches. Flipping one changes results -- do so deliberately.

    Attributes:
        rhoN_T0_K: stagnation temperature used by the density calculation only.
            The original hard-coded 300 K inside ``calculateRhoN`` while velocity
            and temperature used the configured value, so in the resolution-study
            cases density was computed at 300 K and velocity written for 3000 K
            (finding SM-02). ``None`` means "use ``stagnation.T0_K``", i.e. the
            corrected behaviour.
        clip_beyond_limiting_angle: zero the density past the limiting angle.
            ``False`` preserves the folding ``|cos|**n`` profile and the NaN path
            in ``f_phi`` (finding SM-06).
        emit_absent_patches: allow ``output.patches`` to name patches the mesh does
            not have, warning instead of failing. The original hard-coded
            ``cylinder`` and ``plate`` blocks into every generated field file, on
            meshes that have neither (finding AD-03); ``True`` reproduces that.
            Set ``False`` for new cases so a typo in a patch name is an error.
    """

    rhoN_T0_K: float | None = 300.0
    clip_beyond_limiting_angle: bool = False
    emit_absent_patches: bool = True


@dataclass(frozen=True)
class CaseConfig:
    model: str = "source_flow"
    mesh: MeshConfig = field(default_factory=MeshConfig)
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    gas: GasConfig = field(default_factory=GasConfig)
    stagnation: StagnationConfig = field(default_factory=StagnationConfig)
    angular: AngularConfig = field(default_factory=AngularConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    sampling: SamplingConfig = field(default_factory=SamplingConfig)
    legacy: LegacyConfig = field(default_factory=LegacyConfig)
    # --- AIAA 99-3455 sections ------------------------------------------------
    # Optional, and defaulted, so every existing case.yaml keeps loading unchanged.
    # A case that does not set them simply never reads them.
    bodies: BodiesConfig = field(default_factory=BodiesConfig)
    dsmc: DsmcConfig = field(default_factory=DsmcConfig)
    resolution: ResolutionConfig = field(default_factory=ResolutionConfig)
    checks: ChecksConfig = field(default_factory=ChecksConfig)
    post: PostConfig = field(default_factory=PostConfig)
    #: Free-form provenance carried into manifests and case summaries. Not
    #: validated: it records where a case came from, and over-constraining that
    #: would make recording an awkward fact harder than omitting it.
    meta: dict = field(default_factory=dict)


_SECTIONS = {
    "mesh": MeshConfig,
    "geometry": GeometryConfig,
    "gas": GasConfig,
    "stagnation": StagnationConfig,
    "angular": AngularConfig,
    "output": OutputConfig,
    "sampling": SamplingConfig,
    "legacy": LegacyConfig,
    "bodies": BodiesConfig,
    "dsmc": DsmcConfig,
    "resolution": ResolutionConfig,
    "checks": ChecksConfig,
    "post": PostConfig,
}

#: Top-level keys that are passed through rather than built into a dataclass.
_PASSTHROUGH = ("model", "meta")


def _build(cls, data: dict, where: str):
    """Instantiate a config dataclass, rejecting unknown keys.

    Nested sections -- ``bodies.cylinder``, ``dsmc.species`` -- are built
    recursively. A class declares them by carrying a ``_nested`` mapping of field
    name to class; the attribute is unannotated, so ``dataclasses`` ignores it and
    it never becomes a configurable key itself.

    Recursion is opt-in rather than inferred from the type annotation because
    ``from __future__ import annotations`` turns every annotation into a string,
    and resolving those would mean either ``eval`` or ``get_type_hints`` against a
    module namespace -- more machinery than an explicit two-entry mapping.
    """
    if not isinstance(data, dict):
        raise ConfigError(f"{where}: expected a mapping, got {type(data).__name__}")
    known = {f.name for f in fields(cls)}
    unknown = set(data) - known
    if unknown:
        raise ConfigError(
            f"{where}: unknown key(s) {sorted(unknown)}; valid keys are {sorted(known)}"
        )
    nested = getattr(cls, "_nested", {})
    coerced = {}
    for f in fields(cls):
        if f.name not in data:
            continue
        value = data[f.name]
        if f.name in nested:
            value = _build(nested[f.name], value or {}, f"{where}.{f.name}")
        elif f.type in ("tuple", tuple) and isinstance(value, list):
            value = tuple(value)
        elif isinstance(value, str):
            value = _coerce_number(value, f, where)
        coerced[f.name] = value
    return cls(**coerced)


def _coerce_number(value: str, f, where: str):
    """Turn a numeric-looking string into a number for a numeric field.

    YAML 1.1 requires the exponent of a float to carry a sign, so ``1.0e-14``
    parses as a float but ``1.0e14`` parses as a **string**. That is a genuine
    trap: the value looks right in the file, loads without complaint, and then
    fails deep inside a dictionary renderer with ``Unknown format code 'g' for
    object of type 'str'``, a long way from the line that caused it.

    So a string reaching a field declared ``float`` or ``int`` is converted here,
    or rejected with the key name and the fix.
    """
    annotation = str(f.type)
    if "float" not in annotation and "int" not in annotation:
        return value
    try:
        return float(value) if "float" in annotation else int(value)
    except ValueError:
        raise ConfigError(
            f"{where}: {f.name} is {value!r}, which is not a number. If it looks "
            f"like one, check the exponent: YAML requires a signed exponent, so "
            f"1.0e+14 is a float and 1.0e14 is a string."
        ) from None


def load_case_config(case_dir: Path) -> CaseConfig:
    """Read and validate ``<case_dir>/case.yaml``.

    Args:
        case_dir: an OpenFOAM case directory containing ``case.yaml``.

    Returns:
        A validated :class:`CaseConfig`.

    Raises:
        ConfigError: if the file is missing, contains unknown keys, or fails one
            of the cross-checks in :func:`validate_against_case`.
    """
    case_dir = Path(case_dir)
    path = case_dir / "case.yaml"
    if not path.is_file():
        raise ConfigError(f"no case.yaml at {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: expected a mapping at the top level")

    unknown = set(raw) - set(_SECTIONS) - set(_PASSTHROUGH)
    if unknown:
        raise ConfigError(
            f"{path}: unknown top-level key(s) {sorted(unknown)}; "
            f"valid: {sorted(set(_SECTIONS) | set(_PASSTHROUGH))}"
        )

    sections: dict[str, Any] = {
        name: _build(cls, raw.get(name, {}) or {}, f"{path}:{name}")
        for name, cls in _SECTIONS.items()
    }
    cfg = CaseConfig(
        model=raw.get("model", "source_flow"),
        meta=raw.get("meta") or {},
        **sections,
    )

    validate_against_case(cfg, case_dir, path)
    return cfg


#: Patch roles ``mesh.patch_names`` must fill for ``snappy_markelov``.
#:
#: ``upstream_vacuum`` is the ``x = 0`` plane outside the inflow cavity. It is an
#: open boundary, **not** a symmetry plane: the configuration is symmetric about
#: ``y = 0`` only, and calling this one "symmetry" would both misdescribe it and
#: reflect back any particle that scattered upstream.
MARKELOV_PATCH_ROLES = (
    "inflow", "cylinder", "plate", "outer", "symmetry", "upstream_vacuum",
)


def _validate_markelov(cfg: CaseConfig, path: Path) -> None:
    """Cross-checks specific to ``mesh.type: snappy_markelov``.

    Everything here would otherwise surface as a snappyHexMesh failure with no
    indication of which configuration key caused it, or -- worse -- as a mesh that
    builds cleanly around the wrong geometry.
    """
    m = cfg.mesh

    missing_extents = [
        name for name in ("x_min_m", "x_max_m", "y_min_m", "y_max_m", "z_min_m", "z_max_m")
        if getattr(m, name) is None
    ]
    if missing_extents:
        raise ConfigError(
            f"{path}: mesh.type is 'snappy_markelov' but mesh.{missing_extents} "
            f"{'are' if len(missing_extents) > 1 else 'is'} unset. This geometry is "
            f"neither cubic nor symmetric about y, so no extent can be inferred from "
            f"another; all six are required."
        )

    for lo_name, hi_name in (("x_min_m", "x_max_m"),
                             ("y_min_m", "y_max_m"),
                             ("z_min_m", "z_max_m")):
        lo, hi = getattr(m, lo_name), getattr(m, hi_name)
        if not hi > lo:
            raise ConfigError(
                f"{path}: mesh.{hi_name} ({hi}) must exceed mesh.{lo_name} ({lo})")

    missing_roles = [r for r in MARKELOV_PATCH_ROLES if r not in m.patch_names]
    if missing_roles:
        raise ConfigError(
            f"{path}: mesh.patch_names is missing role(s) {missing_roles}; "
            f"'snappy_markelov' needs all of {list(MARKELOV_PATCH_ROLES)}"
        )
    names = [m.patch_names[r] for r in MARKELOV_PATCH_ROLES]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ConfigError(
            f"{path}: mesh.patch_names reuses {duplicates} for more than one role; "
            f"two roles sharing a patch name would silently merge two boundaries "
            f"with different physics"
        )

    if cfg.model != "markelov1999_axisymmetric":
        raise ConfigError(
            f"{path}: mesh.type 'snappy_markelov' requires "
            f"model: markelov1999_axisymmetric, got {cfg.model!r}. The legacy "
            f"'source_flow' model measures its angle from +z and is frozen to "
            f"reproduce known defects; it must not be run on this geometry."
        )

    # The source-flow equations take the PHYSICAL orifice radius. Confusing it
    # with the hemispherical inflow radius scales the density by (R/r_e)**2 --
    # about 1.4e5 here -- so the two are required to differ by a wide margin.
    if cfg.stagnation.throat_radius_m >= cfg.geometry.sphere_radius_m:
        raise ConfigError(
            f"{path}: stagnation.throat_radius_m ({cfg.stagnation.throat_radius_m}) "
            f"must be smaller than geometry.sphere_radius_m "
            f"({cfg.geometry.sphere_radius_m}). The first is the physical orifice "
            f"the source-flow equations use; the second is the computational "
            f"surface the model is evaluated on. They are not interchangeable."
        )

    if cfg.resolution.sizing_region not in ("cylinder", "plate", "wake"):
        raise ConfigError(
            f"{path}: unknown resolution.sizing_region "
            f"{cfg.resolution.sizing_region!r}; valid: ['cylinder', 'plate', 'wake']"
        )

    # gamma and the internal degrees of freedom are two statements of the same
    # physics; a case that sets one without the other gets a silently wrong
    # limiting velocity, which is finding SM-03 in a new costume.
    dof = 3 + int(cfg.dsmc.species.internal_degrees_of_freedom)
    implied_gamma = (dof + 2) / dof
    if abs(implied_gamma - cfg.dsmc.species.gamma) > 1e-6:
        raise ConfigError(
            f"{path}: dsmc.species.gamma is {cfg.dsmc.species.gamma} but "
            f"internal_degrees_of_freedom {cfg.dsmc.species.internal_degrees_of_freedom} "
            f"implies {implied_gamma:.6f} (3 translational + "
            f"{cfg.dsmc.species.internal_degrees_of_freedom} internal). Set both "
            f"consistently: the model uses gamma, the solver uses the DoF count."
        )


def validate_against_case(cfg: CaseConfig, case_dir: Path, path: Path) -> None:
    """Cross-check the config against the mesh and the OpenFOAM dictionaries.

    Checks, in order:

    1. ``mesh.sphere_radius_m == geometry.sphere_radius_m``. These feed the mesh
        generator and the model respectively; if they disagree the model is
        evaluated on a different sphere than the one that was meshed (SM-09).
    2. ``geometry.patch`` exists in ``constant/polyMesh/boundary``.
    3. Every ``output.patches`` name exists in the mesh -- an error, unless
        ``legacy.emit_absent_patches`` is set, in which case it is a warning. The
        original hard-coded ``cylinder`` and ``plate`` entries into every generated
        field file, on meshes that have neither (finding AD-03).
    4. ``gas.species_name`` appears in ``dsmcProperties``' ``typeIdList``.

    Checks 2-4 are skipped when the corresponding file is absent, so a config can
    be loaded before the mesh has been generated.
    """
    if cfg.mesh.sphere_radius_m != cfg.geometry.sphere_radius_m:
        raise ConfigError(
            f"{path}: mesh.sphere_radius_m ({cfg.mesh.sphere_radius_m}) must equal "
            f"geometry.sphere_radius_m ({cfg.geometry.sphere_radius_m}); the mesh "
            f"generator and the source-flow model would otherwise use different spheres"
        )

    if cfg.mesh.type not in MESH_TYPES:
        raise ConfigError(
            f"{path}: unknown mesh.type {cfg.mesh.type!r}; valid: {list(MESH_TYPES)}"
        )
    if cfg.mesh.projection not in PROJECTIONS:
        raise ConfigError(
            f"{path}: unknown mesh.projection {cfg.mesh.projection!r}; "
            f"valid: {list(PROJECTIONS)}"
        )
    if cfg.output.inflow_model not in INFLOW_MODELS:
        raise ConfigError(
            f"{path}: unknown output.inflow_model {cfg.output.inflow_model!r}; "
            f"valid: {list(INFLOW_MODELS)}"
        )

    # The FreeStream warning below is specifically about FreeStream's missing
    # patch-selection list. A model that takes an explicit patch list does not
    # have the problem, and warning anyway would train readers to ignore it.
    if (cfg.output.dialect == "standard"
            and cfg.output.inflow_model == "FreeStream"
            and cfg.mesh.outer_patch_type == "patch"):
        warnings.warn(
            f"{path}: output.dialect is 'standard' and mesh.outer_patch_type is "
            f"'patch'. Standard dsmcFoam's FreeStream injects on EVERY patch-type "
            f"boundary (FreeStream.C:57-62), so the outer boundary becomes a second "
            f"inflow and the run aborts at the first timestep with 'Zero boundary "
            f"temperature detected'. Set mesh.outer_patch_type: wall to exclude it "
            f"-- at the cost of a reflecting, non-absorbing outer boundary. See "
            f"docs/solver-compatibility.md.",
            ConfigWarning,
            stacklevel=2,
        )

    if cfg.mesh.type == "snappy_markelov":
        _validate_markelov(cfg, path)

    if cfg.mesh.type == "block_mesh_ogrid" and cfg.mesh.projection in ("arc", "none"):
        warnings.warn(
            f"{path}: mesh.projection is {cfg.mesh.projection!r}, so only the twelve "
            f"block edges are curved and the inflow patch is NOT a true hemisphere -- "
            f"the radial deficit reaches 16.3% of R at the cap face centre and does "
            f"not shrink with n_tangential. Set mesh.projection: searchable_sphere "
            f"for an exact surface. See docs/mesh.md.",
            ConfigWarning,
            stacklevel=2,
        )

    boundary = case_dir / "constant" / "polyMesh" / "boundary"
    if boundary.is_file():
        patches = read_boundary(case_dir)
        if cfg.geometry.patch not in patches:
            raise ConfigError(
                f"{path}: geometry.patch {cfg.geometry.patch!r} is not in the mesh; "
                f"available patches: {sorted(patches)}"
            )
        missing = sorted(set(cfg.output.patches) - set(patches))
        if missing and not cfg.legacy.emit_absent_patches:
            raise ConfigError(
                f"{path}: output.patches names {missing} are not in the mesh; "
                f"available patches: {sorted(patches)}"
            )
        if missing:
            warnings.warn(
                f"{path}: output.patches names {missing} are not in the mesh "
                f"{sorted(patches)}, but legacy.emit_absent_patches is set so they "
                f"will still be written (finding AD-03)",
                ConfigWarning,
                stacklevel=2,
            )

    dsmc = case_dir / "constant" / "dsmcProperties"
    if dsmc.is_file():
        text = dsmc.read_text(encoding="utf-8", errors="replace")
        match = re.search(r"typeIdList\s*\(([^)]*)\)", text)
        if match:
            type_ids = match.group(1).split()
            if cfg.gas.species_name not in type_ids:
                raise ConfigError(
                    f"{path}: gas.species_name {cfg.gas.species_name!r} is not in "
                    f"dsmcProperties typeIdList {type_ids}"
                )
