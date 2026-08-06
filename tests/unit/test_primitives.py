"""The shared blockMesh/snappyHexMesh building blocks.

These are exercised through both case families, but the properties below are
about the primitives themselves -- in particular the block-ordering fix, which is
invisible on a cubic domain and therefore needs a non-cubic one to state.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from plumetools.foamio.primitives import (
    BOX_FACES,
    BoxPatch,
    SearchableBox,
    SearchableCylinder,
    SearchableSphere,
    check_surfaces_are_resolvable,
    fmt,
    fmt_point,
    render_box_block_mesh_dict,
    render_geometry_block,
    render_refinement_surfaces_block,
    render_snappy_dict,
    uniform_divisions,
)

GEN = "tests.unit.test_primitives"
BOUNDS = (0.0, 0.9, 0.0, 0.3, -0.35, 0.35)


def box_dict(**kwargs):
    params = dict(
        bounds=BOUNDS,
        divisions=(36, 12, 28),
        patches=[
            BoxPatch("symmetry", "symmetry", ("y_min",)),
            BoxPatch("upstreamVacuum", "patch", ("x_min",)),
            BoxPatch("vacuum", "patch", ("x_max", "y_max", "z_min", "z_max")),
        ],
        generator=GEN,
        note="test box",
    )
    params.update(kwargs)
    return render_box_block_mesh_dict(**params)


# --------------------------------------------------------------------------- #
# formatting
# --------------------------------------------------------------------------- #

def test_lengths_keep_ten_significant_figures():
    """A 0.4 mm orifice radius in a 0.9 m domain has to survive the round trip."""
    assert fmt(0.00041275) == "0.00041275"
    assert fmt(0.52705) == "0.52705"
    assert float(fmt(0.29845 + 0.0762 + 0.1524)) == pytest.approx(0.52705, rel=1e-12)


def test_points_render_as_openfoam_vectors():
    assert fmt_point((0.0, 0.0, 0.0)) == "(0 0 0)"
    assert fmt_point((0.29845, 0.0, -0.2286)) == "(0.29845 0 -0.2286)"


# --------------------------------------------------------------------------- #
# the block ordering fix
# --------------------------------------------------------------------------- #

def test_block_local_axes_are_x_then_y_then_z():
    """The ordering this replaced ran its local x1 along +y and x3 along +x, so a
    (nx ny nz) triple was applied to the wrong axes. That was invisible on the
    only case using it, which is a cube with all three counts equal.

    blockMesh defines x1 as vertex 0 -> 1, x2 as 1 -> 2 and x3 as 0 -> 4. Check
    those three directions directly against the emitted vertex list.
    """
    text = box_dict()
    block = text[text.index("vertices"):text.index("blocks")]
    pts = np.array([[float(v) for v in m] for m in re.findall(
        r"\(\s*(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*\)", block)])
    assert len(pts) == 8

    np.testing.assert_allclose(pts[1] - pts[0], [0.9, 0.0, 0.0], atol=1e-12)   # x1 = +x
    np.testing.assert_allclose(pts[2] - pts[1], [0.0, 0.3, 0.0], atol=1e-12)   # x2 = +y
    np.testing.assert_allclose(pts[4] - pts[0], [0.0, 0.0, 0.7], atol=1e-12)   # x3 = +z


def test_divisions_land_on_the_axes_they_are_named_for():
    assert "(36 12 28)" in box_dict()


def test_uniform_divisions_give_near_cubic_cells():
    nx, ny, nz = uniform_divisions(BOUNDS, 0.025)
    assert (nx, ny, nz) == (36, 12, 28)
    for span, n in ((0.9, nx), (0.3, ny), (0.7, nz)):
        assert span / n == pytest.approx(0.025, rel=0.02)


def test_uniform_divisions_never_return_zero():
    assert uniform_divisions((0.0, 0.001, 0.0, 0.001, 0.0, 0.001), 1.0) == (1, 1, 1)


def test_uniform_divisions_rejects_a_nonpositive_size():
    with pytest.raises(ValueError, match="must be positive"):
        uniform_divisions(BOUNDS, 0.0)


# --------------------------------------------------------------------------- #
# box patches
# --------------------------------------------------------------------------- #

def test_every_box_face_must_be_claimed_exactly_once():
    """An unclaimed face is not a closed domain, and blockMesh's defaultFaces
    would absorb it under a name nothing configures."""
    with pytest.raises(ValueError, match="not assigned to any patch"):
        box_dict(patches=[BoxPatch("vacuum", "patch", ("x_min",))])


def test_a_face_claimed_twice_is_rejected():
    with pytest.raises(ValueError, match="claimed by both"):
        box_dict(patches=[
            BoxPatch("a", "patch", tuple(BOX_FACES)),
            BoxPatch("b", "patch", ("x_min",)),
        ])


def test_an_unknown_face_key_is_rejected():
    with pytest.raises(ValueError, match="unknown box face"):
        box_dict(patches=[BoxPatch("a", "patch", ("x_middle",))])


def test_patch_types_are_written_verbatim():
    text = box_dict()
    assert re.search(r"symmetry\s*\{\s*type symmetry;", text)
    assert re.search(r"upstreamVacuum\s*\{\s*type patch;", text)
    assert re.search(r"vacuum\s*\{\s*type patch;", text)


def test_the_box_has_exactly_six_faces_across_all_patches():
    body = box_dict()
    body = body[body.index("boundary"):body.index("mergePatchPairs")]
    assert len(re.findall(r"^\s*\(\d+ \d+ \d+ \d+\)\s*//", body, re.M)) == 6


def test_box_face_windings_all_point_outward():
    """OpenFOAM requires a boundary face wound so its normal leaves the domain.
    Reconstruct each face from the canonical vertex list and check the sign."""
    x0, x1, y0, y1, z0, z1 = BOUNDS
    verts = np.array([
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ])
    centre = verts.mean(axis=0)
    for key, quad in BOX_FACES.items():
        p = verts[list(quad)]
        normal = np.cross(p[1] - p[0], p[2] - p[0])
        outward = p.mean(axis=0) - centre
        assert np.dot(normal, outward) > 0, f"{key} is wound inward"


def test_inverted_bounds_are_rejected():
    with pytest.raises(ValueError, match="y bounds are inverted"):
        box_dict(bounds=(0.0, 0.9, 0.3, 0.0, -0.35, 0.35))


def test_nonpositive_divisions_are_rejected():
    with pytest.raises(ValueError, match="divisions must be >= 1"):
        box_dict(divisions=(36, 0, 28))


# --------------------------------------------------------------------------- #
# searchable surfaces
# --------------------------------------------------------------------------- #

def test_sphere_renders_as_a_primitive_with_no_stl():
    text = "\n".join(render_geometry_block(
        [SearchableSphere(name="inflow", level=2, radius=0.1524)]))
    assert "type    searchableSphere;" in text
    assert "radius  0.1524;" in text
    assert "centre  (0 0 0);" in text
    assert ".stl" not in text and "triSurfaceMesh" not in text


def test_cylinder_renders_with_both_end_caps():
    """Finite and capped: the paper's cylinder has ends inside the domain, and a
    cylinder running to the boundary would delete the 3D wake around them."""
    text = "\n".join(render_geometry_block([SearchableCylinder(
        name="cylinder", level=2, radius=0.0762,
        point1=(0.29845, 0.0, -0.2286), point2=(0.29845, 0.0, 0.2286))]))
    assert "type    searchableCylinder;" in text
    assert "point1  (0.29845 0 -0.2286);" in text
    assert "point2  (0.29845 0 0.2286);" in text
    assert "radius  0.0762;" in text


def test_cylinder_length_is_the_distance_between_caps():
    cyl = SearchableCylinder(name="c", point1=(0.3, 0.0, -0.2286),
                             point2=(0.3, 0.0, 0.2286), radius=0.0762)
    assert cyl.length_m() == pytest.approx(0.4572, rel=1e-12)


def test_box_renders_as_min_and_max_corners():
    text = "\n".join(render_geometry_block([SearchableBox(
        name="plate", level=3, min=(0.52705, -0.0762, -0.1905),
        max=(0.53975, 0.0762, 0.1905))]))
    assert "type    searchableBox;" in text
    assert "min     (0.52705 -0.0762 -0.1905);" in text
    assert "max     (0.53975 0.0762 0.1905);" in text


def test_box_min_feature_is_its_thinnest_dimension():
    """For a plate that is the thickness -- the dimension a background mesh is
    most likely to miss."""
    plate = SearchableBox(name="plate", min=(0.52705, -0.0762, -0.1905),
                          max=(0.53975, 0.0762, 0.1905))
    assert plate.min_feature_size_m() == pytest.approx(0.0127, rel=1e-9)


def test_degenerate_surfaces_are_rejected():
    with pytest.raises(ValueError, match="sphere radius must be positive"):
        render_geometry_block([SearchableSphere(name="s", radius=0.0)])
    with pytest.raises(ValueError, match="end caps coincide"):
        render_geometry_block([SearchableCylinder(
            name="c", point1=(0, 0, 0), point2=(0, 0, 0), radius=1.0)])
    with pytest.raises(ValueError, match="box axis z is degenerate"):
        render_geometry_block([SearchableBox(
            name="b", min=(0, 0, 0), max=(1, 1, 0))])
    with pytest.raises(ValueError, match="refinement level must be >= 0"):
        render_geometry_block([SearchableSphere(name="s", level=-1, radius=1.0)])


def test_duplicate_surface_names_are_rejected():
    """Dictionary keys are unique, so a duplicate would silently drop one body."""
    with pytest.raises(ValueError, match="duplicate searchable surface"):
        render_geometry_block([SearchableSphere(name="x", radius=1.0),
                               SearchableBox(name="x")])


# --------------------------------------------------------------------------- #
# background resolution guard
# --------------------------------------------------------------------------- #

def test_a_thin_plate_is_resolvable_when_its_own_level_is_high_enough():
    """A 12.7 mm plate under a 25 mm background is fine at level 3 (3.125 mm, four
    cells across) and not at level 1 (12.5 mm, one cell). Coarseness relative to
    the *thickness* is not by itself a problem: the plate is 400 mm tall, so it
    intersects hundreds of background cells and cannot be missed."""
    thin = dict(min=(0.5, -0.1, -0.2), max=(0.5127, 0.1, 0.2))
    check_surfaces_are_resolvable([SearchableBox(name="plate", level=3, **thin)], 0.025)
    with pytest.raises(ValueError, match="only 1.02 cells across"):
        check_surfaces_are_resolvable([SearchableBox(name="plate", level=1, **thin)], 0.025)


def test_the_resolution_error_names_the_level_that_would_fix_it():
    """25 mm background, 12.7 mm plate: two cells across needs 6.35 mm, so
    ceil(log2(25/6.35)) = 2. Checked by taking the suggestion and re-running."""
    thin = dict(min=(0.5, -0.1, -0.2), max=(0.5127, 0.1, 0.2))
    with pytest.raises(ValueError, match=r"Raise its refinement level to 2"):
        check_surfaces_are_resolvable([SearchableBox(name="plate", level=1, **thin)], 0.025)
    check_surfaces_are_resolvable([SearchableBox(name="plate", level=2, **thin)], 0.025)


def test_a_body_smaller_than_one_background_cell_is_rejected():
    """The genuine "cannot be found" case: every dimension inside one cell, so the
    body need never intersect a cell boundary."""
    speck = SearchableSphere(name="speck", level=4, radius=0.005)
    with pytest.raises(ValueError, match="whole body can sit inside one background cell"):
        check_surfaces_are_resolvable([speck], 0.025)


def test_the_guard_names_the_offending_surface():
    surfaces = [
        SearchableSphere(name="inflow", level=2, radius=0.1524),
        SearchableBox(name="plate", level=1, min=(0.5, -0.1, -0.2), max=(0.5127, 0.1, 0.2)),
    ]
    with pytest.raises(ValueError, match="'plate'"):
        check_surfaces_are_resolvable(surfaces, 0.025)


def test_the_shipped_baseline_combination_is_resolvable():
    """background 25 mm; inflow and cylinder at level 2, plate at level 3."""
    check_surfaces_are_resolvable([
        SearchableSphere(name="inflow", level=2, radius=0.1524),
        SearchableCylinder(name="cylinder", level=2, radius=0.0762,
                           point1=(0.29845, 0.0, -0.2286),
                           point2=(0.29845, 0.0, 0.2286)),
        SearchableBox(name="plate", level=3, min=(0.52705, -0.0762, -0.1905),
                      max=(0.53975, 0.0762, 0.1905)),
    ], 0.025)


# --------------------------------------------------------------------------- #
# the assembled snappyHexMeshDict
# --------------------------------------------------------------------------- #

@pytest.fixture
def surfaces():
    return [
        SearchableSphere(name="inflow", level=2, patch_type="patch", radius=0.1524),
        SearchableCylinder(name="cylinder", level=2, patch_type="wall", radius=0.0762,
                           point1=(0.29845, 0.0, -0.2286),
                           point2=(0.29845, 0.0, 0.2286)),
        SearchableBox(name="plate", level=3, patch_type="wall",
                      min=(0.52705, -0.0762, -0.1905), max=(0.53975, 0.0762, 0.1905)),
    ]


def test_every_surface_gets_a_geometry_and_a_refinement_entry(surfaces):
    text = render_snappy_dict(surfaces, location_in_mesh=(0.7, 0.15, 0.3),
                              n_cells_between_levels=2, generator=GEN, note="t")
    geometry = text[text.index("geometry"):text.index("castellatedMeshControls")]
    refinement = text[text.index("refinementSurfaces"):text.index("resolveFeatureAngle")]
    for name in ("inflow", "cylinder", "plate"):
        assert re.search(rf"^\s*{name}$", geometry, re.M), f"{name} missing from geometry"
        assert re.search(rf"^\s*{name}$", refinement, re.M), f"{name} missing from refinement"


def test_refinement_levels_are_per_surface(surfaces):
    text = "\n".join(render_refinement_surfaces_block(surfaces))
    assert text.count("level (2 2);") == 2   # inflow and cylinder
    assert text.count("level (3 3);") == 1   # plate


def test_bodies_are_walls_and_the_inflow_is_a_patch(surfaces):
    """Load-bearing for DSMC. A `patch` deletes every particle that reaches it, so
    a body declared as `patch` would absorb the plume and exert no force; a `wall`
    runs hitWallPatch, which is what records fD and applies the reflection model.
    An inflow declared as `wall` would do neither injection nor deletion."""
    text = "\n".join(render_refinement_surfaces_block(surfaces))
    blocks = re.findall(r"(\w+)\s*\{\s*level \(\d+ \d+\);\s*patchInfo\s*\{\s*type (\w+);",
                        text)
    assert dict(blocks) == {"inflow": "patch", "cylinder": "wall", "plate": "wall"}


def test_layers_are_off_and_snap_is_on(surfaces):
    text = render_snappy_dict(surfaces, location_in_mesh=(0.7, 0.15, 0.3),
                              n_cells_between_levels=2, generator=GEN, note="t")
    assert "castellatedMesh true;" in text
    assert "snap            true;" in text
    assert "addLayers       false;" in text


def test_refinement_regions_block_is_present_but_empty(surfaces):
    """Requirement: no wake/shock/gap refinement yet, but a clean place to add it."""
    text = render_snappy_dict(surfaces, location_in_mesh=(0.7, 0.15, 0.3),
                              n_cells_between_levels=2, generator=GEN, note="t")
    body = re.search(r"refinementRegions\s*\{\s*\}", text)
    assert body, "refinementRegions must exist and be empty"


def test_mesh_quality_dict_is_included(surfaces):
    text = render_snappy_dict(surfaces, location_in_mesh=(0.7, 0.15, 0.3),
                              n_cells_between_levels=2, generator=GEN, note="t")
    assert '#include "meshQualityDict"' in text


def test_negative_buffer_cells_are_rejected(surfaces):
    with pytest.raises(ValueError, match="n_cells_between_levels"):
        render_snappy_dict(surfaces, location_in_mesh=(0.7, 0.15, 0.3),
                           n_cells_between_levels=-1, generator=GEN, note="t")


def test_generated_dictionaries_name_their_generator(surfaces):
    text = render_snappy_dict(surfaces, location_in_mesh=(0.7, 0.15, 0.3),
                              n_cells_between_levels=2, generator=GEN, note="t")
    assert f"GENERATED by {GEN}" in text
    assert f"GENERATED by {GEN}" in box_dict()
