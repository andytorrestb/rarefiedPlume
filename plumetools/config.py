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
MESH_TYPES = ("block_mesh_ogrid", "snappy_hex_sphere")

#: How the inflow surface is made curved, for ``block_mesh_ogrid``.
PROJECTIONS = ("arc", "searchable_sphere", "none")


@dataclass(frozen=True)
class MeshConfig:
    """Parameters for generating the mesh.

    Attributes:
        type: ``"block_mesh_ogrid"`` -- a 5-block hexahedral O-grid, rendered by
            :mod:`plumetools.foamio.blockmesh`; or ``"snappy_hex_sphere"`` -- a
            background box carved by snappyHexMesh, rendered by
            :mod:`plumetools.foamio.snappy`.
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
        background_cell_size_m: snappy only. Background box cell size [m]; must be
            at most about R/4 or the cavity falls between cells.
        refinement_level: snappy only. Octree levels at the sphere.
        n_cells_between_levels: snappy only. Buffer cells between refinement levels.
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
    background_cell_size_m: float = 0.125
    refinement_level: int = 3
    n_cells_between_levels: int = 3
    patch_names: dict = field(
        default_factory=lambda: {"inflow": "inflow", "outer": "vacuum", "symmetry": "sym"}
    )


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


_SECTIONS = {
    "mesh": MeshConfig,
    "geometry": GeometryConfig,
    "gas": GasConfig,
    "stagnation": StagnationConfig,
    "angular": AngularConfig,
    "output": OutputConfig,
    "sampling": SamplingConfig,
    "legacy": LegacyConfig,
}


def _build(cls, data: dict, where: str):
    """Instantiate a config dataclass, rejecting unknown keys."""
    if not isinstance(data, dict):
        raise ConfigError(f"{where}: expected a mapping, got {type(data).__name__}")
    known = {f.name for f in fields(cls)}
    unknown = set(data) - known
    if unknown:
        raise ConfigError(
            f"{where}: unknown key(s) {sorted(unknown)}; valid keys are {sorted(known)}"
        )
    coerced = {}
    for f in fields(cls):
        if f.name not in data:
            continue
        value = data[f.name]
        if f.type in ("tuple", tuple) and isinstance(value, list):
            value = tuple(value)
        coerced[f.name] = value
    return cls(**coerced)


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

    unknown = set(raw) - set(_SECTIONS) - {"model"}
    if unknown:
        raise ConfigError(
            f"{path}: unknown top-level key(s) {sorted(unknown)}; "
            f"valid: {sorted(set(_SECTIONS) | {'model'})}"
        )

    sections: dict[str, Any] = {
        name: _build(cls, raw.get(name, {}) or {}, f"{path}:{name}")
        for name, cls in _SECTIONS.items()
    }
    cfg = CaseConfig(model=raw.get("model", "source_flow"), **sections)

    validate_against_case(cfg, case_dir, path)
    return cfg


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
