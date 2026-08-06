"""The AIAA 99-3455 mesh setup: domain, patches, primitives, seed point.

Checked as text and as geometry, with no OpenFOAM. Whether snappyHexMesh accepts
the dictionaries is a ``needs_openfoam`` question, covered in
``tests/openfoam/``.
"""

from __future__ import annotations

import math
import re
from dataclasses import replace

import numpy as np
import pytest

from plumetools.config import (
    BodiesConfig,
    CaseConfig,
    ConfigError,
    MeshConfig,
    StagnationConfig,
)
from plumetools.foamio import mesh_pipeline, write_mesh_setup
from plumetools.markelov1999 import mesh as mm
from plumetools.markelov1999.geometry import from_config

PATCH_NAMES = {
    "inflow": "inflow", "cylinder": "cylinder", "plate": "plate",
    "outer": "vacuum", "symmetry": "symmetry", "upstream_vacuum": "upstreamVacuum",
}


def make_cfg(**mesh_overrides) -> CaseConfig:
    params = dict(
        type="snappy_markelov",
        sphere_radius_m=0.1524,
        x_min_m=0.0, x_max_m=0.9,
        y_min_m=0.0, y_max_m=0.3,
        z_min_m=-0.35, z_max_m=0.35,
        background_cell_size_m=0.025,
        inflow_refinement_level=2,
        cylinder_refinement_level=2,
        plate_refinement_level=3,
        n_cells_between_levels=2,
        outer_patch_type="patch",
        patch_names=dict(PATCH_NAMES),
    )
    params.update(mesh_overrides)
    mesh = MeshConfig(**params)
    from plumetools.config import GeometryConfig
    return CaseConfig(
        model="markelov1999_axisymmetric",
        mesh=mesh,
        geometry=GeometryConfig(patch="inflow", sphere_radius_m=0.1524),
        stagnation=StagnationConfig(p0_pa=34473.79, T0_K=300.0,
                                    throat_radius_m=0.00041275),
        bodies=BodiesConfig(),
    )


@pytest.fixture
def cfg():
    return make_cfg()


@pytest.fixture
def geom(cfg):
    return from_config(cfg)


# --------------------------------------------------------------------------- #
# domain
# --------------------------------------------------------------------------- #

def test_domain_bounds_come_from_all_six_keys(cfg):
    assert mm.domain_bounds(cfg) == (0.0, 0.9, 0.0, 0.3, -0.35, 0.35)


def test_background_divisions_and_cell_size(cfg):
    text = mm.render_background_block_mesh_dict(cfg)
    assert "(36 12 28)" in text
    assert mm.background_cell_size(cfg) == pytest.approx(0.7 / 28, rel=1e-9)


def test_surface_cell_sizes_halve_once_per_level(cfg):
    sizes = mm.surface_cell_sizes(cfg)
    bg = mm.background_cell_size(cfg)
    assert sizes["inflow"] == pytest.approx(bg / 4)
    assert sizes["cylinder"] == pytest.approx(bg / 4)
    assert sizes["plate"] == pytest.approx(bg / 8)


def test_plate_thickness_spans_at_least_four_cells(cfg, geom):
    """The plate is the thinnest body; below a few cells across, snappyHexMesh
    cannot carve an interior and it comes out as a dented cell layer."""
    across = geom.plate_thickness_m / mm.surface_cell_sizes(cfg)["plate"]
    assert across >= 4.0


def test_per_surface_levels_fall_back_to_the_shared_one():
    cfg = make_cfg(inflow_refinement_level=None, cylinder_refinement_level=None,
                   plate_refinement_level=None, refinement_level=3)
    assert mm.refinement_levels(cfg) == {"inflow": 3, "cylinder": 3, "plate": 3}


def test_domain_must_contain_both_bodies(cfg, geom):
    mm.check_domain_contains_geometry(cfg, geom)  # baseline is fine

    short = make_cfg(x_max_m=0.5)
    with pytest.raises(ValueError, match="the plate .* reaches or exceeds"):
        mm.check_domain_contains_geometry(short, from_config(short))


def test_cylinder_end_caps_must_be_strictly_inside_the_domain(geom):
    """A cap flush with the z boundary is meshed as though the cylinder ran
    forever, which removes exactly the 3D end effects this case studies."""
    flush = make_cfg(z_min_m=-0.2286, z_max_m=0.2286)
    with pytest.raises(ValueError, match="end caps must be in the domain"):
        mm.check_domain_contains_geometry(flush, from_config(flush))


def test_the_source_must_be_inside_the_domain():
    upstream = make_cfg(x_min_m=0.05, x_max_m=0.9)
    with pytest.raises(ValueError, match="source at x = 0 is outside"):
        mm.check_domain_contains_geometry(upstream, from_config(upstream))


def test_domain_report_states_the_downstream_buffer(cfg, geom):
    text = "\n".join(mm.check_domain_contains_geometry(cfg, geom))
    assert "downstream buffer past the plate" in text
    assert "cylinder diameters" in text


# --------------------------------------------------------------------------- #
# patches
# --------------------------------------------------------------------------- #

def test_all_six_patches_are_declared_across_the_two_dictionaries(cfg):
    block = mm.render_background_block_mesh_dict(cfg)
    snappy = mm.render_snappy_hex_mesh_dict(cfg)

    for name in ("symmetry", "upstreamVacuum", "vacuum"):
        assert re.search(rf"^\s*{name}$", block, re.M), f"{name} missing from blockMesh"
    for name in ("inflow", "cylinder", "plate"):
        assert re.search(rf"^\s*{name}$", snappy, re.M), f"{name} missing from snappy"


def test_the_five_required_patch_names_are_present(cfg):
    """Requirement 4 lists inflow, cylinder, plate, vacuum and symmetry."""
    combined = (mm.render_background_block_mesh_dict(cfg)
                + mm.render_snappy_hex_mesh_dict(cfg))
    for required in ("inflow", "cylinder", "plate", "vacuum", "symmetry"):
        assert required in combined


def test_symmetry_is_on_y_min_and_is_the_only_symmetry_patch(cfg):
    text = mm.render_background_block_mesh_dict(cfg)
    assert re.search(r"symmetry\s*\{\s*type symmetry;", text)
    assert text.count("type symmetry;") == 1
    entry = text[text.index("    symmetry"):text.index("    upstreamVacuum")]
    assert "y_min" in entry
    assert "x_min" not in entry and "z_min" not in entry


def test_upstream_plane_is_an_open_patch_not_a_symmetry_plane(cfg):
    """Requirement 4 names this explicitly. The configuration is symmetric about
    y = 0 and nothing else; calling x = 0 a symmetry plane would reflect back
    every particle that scattered upstream."""
    text = mm.render_background_block_mesh_dict(cfg)
    entry = text[text.index("    upstreamVacuum"):text.index("    vacuum")]
    assert "type patch;" in entry
    assert "symmetry" not in entry.replace("not a symmetry plane", "")
    assert "x_min" in entry


def test_vacuum_takes_the_remaining_four_faces(cfg):
    text = mm.render_background_block_mesh_dict(cfg)
    entry = text[text.index("    vacuum"):text.index("mergePatchPairs")]
    for face in ("x_max", "y_max", "z_min", "z_max"):
        assert face in entry
    assert "type patch;" in entry


def test_bodies_are_walls_and_the_inflow_is_a_patch(cfg):
    text = mm.render_snappy_hex_mesh_dict(cfg)
    refinement = text[text.index("refinementSurfaces"):text.index("resolveFeatureAngle")]
    found = dict(re.findall(
        r"(\w+)\s*\{\s*level \(\d+ \d+\);\s*patchInfo\s*\{\s*type (\w+);", refinement))
    assert found == {"inflow": "patch", "cylinder": "wall", "plate": "wall"}


# --------------------------------------------------------------------------- #
# searchable primitives
# --------------------------------------------------------------------------- #

def test_the_three_primitives_are_native_not_stl(cfg):
    text = mm.render_snappy_hex_mesh_dict(cfg)
    assert "searchableSphere" in text
    assert "searchableCylinder" in text
    assert "searchableBox" in text
    assert ".stl" not in text and "triSurfaceMesh" not in text


def test_surface_order_is_deterministic(cfg):
    """Patch creation order fixes the patch indices in constant/polyMesh/boundary,
    and every decomposed run and stored field ordering depends on them."""
    for _ in range(3):
        assert [s.name for s in mm.searchable_surfaces(cfg)] == [
            "inflow", "cylinder", "plate"]


def test_the_cylinder_primitive_matches_the_paper_geometry(cfg, geom):
    cylinder = mm.searchable_surfaces(cfg)[1]
    assert cylinder.radius == pytest.approx(0.0762, rel=1e-12)
    assert cylinder.length_m() == pytest.approx(0.4572, rel=1e-12)
    assert cylinder.point1 == (0.29845, 0.0, -0.2286)
    assert cylinder.point2 == (0.29845, 0.0, 0.2286)


def test_the_plate_primitive_sits_at_the_derived_x(cfg):
    plate = mm.searchable_surfaces(cfg)[2]
    assert plate.min[0] == pytest.approx(0.52705, rel=1e-9)
    assert plate.max[0] == pytest.approx(0.53975, rel=1e-9)
    assert plate.max[2] - plate.min[2] == pytest.approx(0.381, rel=1e-9)


def test_the_plate_box_straddles_the_symmetry_plane(cfg):
    """A half-width box stopping at y = 0 would leave a spurious internal face
    there; the full box is cut by the domain instead."""
    plate = mm.searchable_surfaces(cfg)[2]
    assert plate.min[1] < 0.0 < plate.max[1]
    assert plate.max[1] - plate.min[1] == pytest.approx(0.1524, rel=1e-9)


def test_the_inflow_sphere_is_centred_on_the_orifice(cfg):
    sphere = mm.searchable_surfaces(cfg)[0]
    assert sphere.centre == (0.0, 0.0, 0.0)
    assert sphere.radius == pytest.approx(0.1524, rel=1e-12)


# --------------------------------------------------------------------------- #
# locationInMesh
# --------------------------------------------------------------------------- #

def test_seed_point_is_outside_every_body(cfg, geom):
    x, y, z = mm.location_in_mesh(cfg, geom)

    assert math.sqrt(x * x + y * y + z * z) > geom.inflow_radius_m
    radial = math.hypot(x - geom.cylinder_centre_x_m, y)
    assert not (geom.cylinder_z_min_m <= z <= geom.cylinder_z_max_m
                and radial <= geom.cylinder_radius_m)
    assert not (geom.plate_upstream_x_m <= x <= geom.plate_downstream_x_m
                and abs(y) <= geom.plate_y_max_m
                and geom.plate_z_min_m <= z <= geom.plate_z_max_m)


def test_seed_point_is_strictly_inside_the_domain_and_off_the_symmetry_plane(cfg):
    x, y, z = mm.location_in_mesh(cfg)
    x0, x1, y0, y1, z0, z1 = mm.domain_bounds(cfg)
    assert x0 < x < x1 and y0 < y < y1 and z0 < z < z1
    assert y != 0.0


@pytest.mark.parametrize("size", [0.05, 0.035, 0.025, 0.02, 0.0125])
def test_seed_point_never_lands_on_a_background_cell_face(size):
    """A fixed fraction of the domain sits exactly on a cell face whenever it
    happens to be a multiple of the cell spacing. Snapping to the containing
    cell's centre makes it interior for any spacing."""
    cfg = make_cfg(background_cell_size_m=size)
    point = mm.location_in_mesh(cfg)
    x0, x1, y0, y1, z0, z1 = mm.domain_bounds(cfg)
    from plumetools.foamio.primitives import uniform_divisions
    nx, ny, nz = uniform_divisions(mm.domain_bounds(cfg), size)

    for value, lo, hi, n in ((point[0], x0, x1, nx),
                             (point[1], y0, y1, ny),
                             (point[2], z0, z1, nz)):
        h = (hi - lo) / n
        assert ((value - lo) / h) % 1.0 == pytest.approx(0.5), (
            "0.0 would put the seed on a cell face")


def test_seed_point_is_written_into_the_dictionary(cfg):
    text = mm.render_snappy_hex_mesh_dict(cfg)
    m = re.search(r"locationInMesh \(([-\d.eE+ ]+)\);", text)
    assert m
    np.testing.assert_allclose(
        [float(v) for v in m.group(1).split()], mm.location_in_mesh(cfg), rtol=1e-6)


# --------------------------------------------------------------------------- #
# writing and dispatch
# --------------------------------------------------------------------------- #

def test_writes_three_dictionaries(cfg, tmp_path):
    paths = mm.write_mesh_setup(tmp_path, cfg)
    assert [p.name for p in paths] == [
        "blockMeshDict", "snappyHexMeshDict", "meshQualityDict"]
    assert all(p.is_file() and p.stat().st_size > 0 for p in paths)
    assert all(p.parent.name == "system" for p in paths)


def test_written_dictionaries_are_lf_terminated(cfg, tmp_path):
    for p in mm.write_mesh_setup(tmp_path, cfg):
        assert b"\r\n" not in p.read_bytes(), p.name


def test_dispatch_selects_this_generator(cfg, tmp_path):
    assert [p.name for p in write_mesh_setup(tmp_path, cfg)] == [
        "blockMeshDict", "snappyHexMeshDict", "meshQualityDict"]
    assert mesh_pipeline(cfg) == ["blockMesh", "snappyHexMesh"]


def test_write_rejects_the_wrong_mesh_type(tmp_path):
    wrong = replace(make_cfg(), mesh=replace(make_cfg().mesh, type="snappy_hex_sphere"))
    with pytest.raises(ValueError, match="snappy_markelov"):
        mm.write_mesh_setup(tmp_path, wrong)


def test_no_wake_or_shock_refinement_regions_yet(cfg):
    """Requirement: refine only the three surfaces for now, but leave a clean
    configuration path."""
    text = mm.render_snappy_hex_mesh_dict(cfg)
    assert re.search(r"refinementRegions\s*\{\s*\}", text)


def test_describe_reports_the_numbers_a_reviewer_needs(cfg):
    text = "\n".join(mm.describe(cfg))
    assert "background" in text
    assert "plate thickness spans" in text
    assert "cylinder circumference spans" in text
    assert "locationInMesh" in text


# --------------------------------------------------------------------------- #
# config validation specific to this mesh type
# --------------------------------------------------------------------------- #

def write_case(tmp_path, text: str):
    case = tmp_path / "case"
    case.mkdir(exist_ok=True)
    (case / "case.yaml").write_text(text, encoding="utf-8")
    return case


BASE_YAML = """
model: markelov1999_axisymmetric
mesh:
  type: snappy_markelov
  sphere_radius_m: 0.1524
  x_min_m: 0.0
  x_max_m: 0.9
  y_min_m: 0.0
  y_max_m: 0.3
  z_min_m: -0.35
  z_max_m: 0.35
  patch_names: {inflow: inflow, cylinder: cylinder, plate: plate,
                outer: vacuum, symmetry: symmetry, upstream_vacuum: upstreamVacuum}
geometry:
  sphere_radius_m: 0.1524
stagnation:
  throat_radius_m: 0.00041275
output:
  inflow_model: plumeFieldInflow
"""


def test_a_complete_markelov_config_loads(tmp_path):
    from plumetools.config import load_case_config
    cfg = load_case_config(write_case(tmp_path, BASE_YAML))
    assert cfg.mesh.type == "snappy_markelov"
    assert cfg.model == "markelov1999_axisymmetric"
    assert cfg.bodies.cylinder.radius_m == pytest.approx(0.0762)
    assert cfg.bodies.plate.thickness_m == pytest.approx(0.0127)


def test_a_missing_domain_extent_is_an_error(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML.replace("  z_max_m: 0.35\n", "")
    with pytest.raises(ConfigError, match=r"mesh.\['z_max_m'\] is unset"):
        load_case_config(write_case(tmp_path, text))


def test_a_missing_patch_role_is_an_error(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML.replace(", upstream_vacuum: upstreamVacuum", "")
    with pytest.raises(ConfigError, match="missing role.*upstream_vacuum"):
        load_case_config(write_case(tmp_path, text))


def test_two_roles_sharing_a_patch_name_is_an_error(tmp_path):
    """Two boundaries with different physics would silently merge into one."""
    from plumetools.config import load_case_config
    text = BASE_YAML.replace("upstream_vacuum: upstreamVacuum",
                             "upstream_vacuum: vacuum")
    with pytest.raises(ConfigError, match="reuses"):
        load_case_config(write_case(tmp_path, text))


def test_the_legacy_model_is_refused_on_this_geometry(tmp_path):
    """The legacy source_flow measures its angle from +z and is frozen to
    reproduce known defects; running it here would look like a valid case."""
    from plumetools.config import load_case_config
    text = BASE_YAML.replace("model: markelov1999_axisymmetric", "model: source_flow")
    with pytest.raises(ConfigError, match="requires\\s+model: markelov1999_axisymmetric"):
        load_case_config(write_case(tmp_path, text))


def test_the_orifice_must_be_smaller_than_the_inflow_surface(tmp_path):
    """Confusing the physical orifice radius with the hemispherical inflow radius
    scales the density by (R/r_e)**2, about 1.4e5 here."""
    from plumetools.config import load_case_config
    text = BASE_YAML.replace("throat_radius_m: 0.00041275", "throat_radius_m: 0.1524")
    with pytest.raises(ConfigError, match="must be smaller than geometry.sphere_radius_m"):
        load_case_config(write_case(tmp_path, text))


def test_gamma_must_agree_with_the_internal_degrees_of_freedom(tmp_path):
    """SM-03 in a new costume: the model uses gamma, the solver uses the DoF
    count, and nothing otherwise ties them together."""
    from plumetools.config import load_case_config
    text = BASE_YAML + "dsmc:\n  species:\n    gamma: 1.4\n    internal_degrees_of_freedom: 0\n"
    with pytest.raises(ConfigError, match="implies 1.666667"):
        load_case_config(write_case(tmp_path, text))


def test_a_valid_monatomic_pairing_is_accepted(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML + ("dsmc:\n  species:\n    gamma: 1.6666666666666667\n"
                        "    internal_degrees_of_freedom: 0\n")
    load_case_config(write_case(tmp_path, text))


def test_the_free_stream_warning_is_not_raised_for_plume_field_inflow(tmp_path):
    """The warning is specifically about FreeStream's missing patch-selection
    list. Raising it for a model that has one would train readers to ignore it."""
    import warnings
    from plumetools.config import ConfigWarning, load_case_config
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConfigWarning)
        load_case_config(write_case(tmp_path, BASE_YAML))


def test_an_unknown_inflow_model_is_an_error(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML.replace("inflow_model: plumeFieldInflow", "inflow_model: madeUp")
    with pytest.raises(ConfigError, match="unknown output.inflow_model"):
        load_case_config(write_case(tmp_path, text))


def test_nested_sections_reject_unknown_keys(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML + "bodies:\n  cylinder:\n    diametre_m: 0.1524\n"
    with pytest.raises(ConfigError, match=r"bodies.cylinder: unknown key"):
        load_case_config(write_case(tmp_path, text))


def test_nested_sections_are_actually_nested(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML + "bodies:\n  gap_m: 0.3048\n  cylinder:\n    centre_x_m: 0.4\n"
    cfg = load_case_config(write_case(tmp_path, text))
    assert cfg.bodies.gap_m == pytest.approx(0.3048)
    assert cfg.bodies.cylinder.centre_x_m == pytest.approx(0.4)
    # untouched keys keep their defaults
    assert cfg.bodies.cylinder.radius_m == pytest.approx(0.0762)


def test_meta_is_carried_through_unvalidated(tmp_path):
    from plumetools.config import load_case_config
    text = BASE_YAML + "meta:\n  reference: AIAA 99-3455\n  pressure_psi: 25\n"
    cfg = load_case_config(write_case(tmp_path, text))
    assert cfg.meta["reference"] == "AIAA 99-3455"
    assert cfg.meta["pressure_psi"] == 25


def test_the_shipped_base_case_yaml_loads_and_meshes(tmp_path):
    """The committed baseCase must stay loadable: it is what generate_cases.py
    clones, so a broken one breaks every generated case at once."""
    from pathlib import Path

    from plumetools.config import load_case_config
    base = Path(__file__).resolve().parents[2] / "cases" / "markelov1999" / "baseCase"
    cfg = load_case_config(base)
    assert cfg.mesh.type == "snappy_markelov"
    paths = mm.write_mesh_setup(tmp_path, cfg)
    assert len(paths) == 3

