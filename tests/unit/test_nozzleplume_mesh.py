"""The mesh plan, the nozzle prediction, the cell estimate and the dictionaries.

Several numbers here are **measurements of a real mesh**, not expectations
worked out by hand: the shipped R-4D case was built with OpenFOAM v2512's
blockMesh, snappyHexMesh, topoSet and createPatch, and came out at 529 126
cells with a 112-face nozzle patch. The plan has to predict the second exactly
and the first closely, and these tests hold it to both without needing OpenFOAM.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from _nozzleplume import R4D_CASE, derive, derived, load_case
from plumetools.nozzleplume import geometry, mesh
from plumetools.nozzleplume.config import load_case_config

#: What OpenFOAM v2512 built from cases/r-4d/baseCase/case.yaml.
R4D_MEASURED_CELLS = 529_126
R4D_MEASURED_CELLS_BY_LEVEL = {0: 89_250, 1: 100_324, 2: 276_320, 3: 36_608,
                               4: 26_624}
R4D_MEASURED_NOZZLE_FACES = 112


@pytest.fixture(scope="module")
def r4d():
    cfg = load_case_config(R4D_CASE)
    return (cfg, *derive(cfg))


# --------------------------------------------------------------------------- #
# the shipped case against the mesh OpenFOAM built from it
# --------------------------------------------------------------------------- #

def test_r4d_nozzle_prediction_is_the_patch_openfoam_built(r4d):
    _, _, _, _, nozzle, _ = r4d
    assert nozzle.n_faces == R4D_MEASURED_NOZZLE_FACES
    assert nozzle.area_m2 == pytest.approx(0.109375)


def test_r4d_nozzle_area_is_within_the_configured_tolerance(r4d):
    cfg, _, _, _, nozzle, _ = r4d
    assert abs(nozzle.area_error) < cfg.checks.max_nozzle_area_error
    assert nozzle.area_error == pytest.approx(-0.0198, abs=2e-4)


def test_r4d_cell_estimate_is_close_to_the_mesh_openfoam_built(r4d):
    """Within 4% overall. An octree's size is not known until it is built, so
    this is a guard, not a count -- but a guard a third out is not one."""
    _, _, _, plan, _, _ = r4d
    assert plan.estimated_cells == pytest.approx(R4D_MEASURED_CELLS, rel=0.04)


def test_r4d_cell_estimate_is_close_at_every_level(r4d):
    """Level 2 is the plate's skin, and the level a naive estimate gets wrong:
    on 811 m^2 of plate, the cells kept at the surface level are half the
    mesh."""
    _, _, _, plan, _, _ = r4d
    for level, measured in R4D_MEASURED_CELLS_BY_LEVEL.items():
        assert plan.cells_by_level[level] == pytest.approx(measured, rel=0.06), level


def test_r4d_background_is_cubic_and_whole(r4d):
    _, geom, _, plan, _, _ = r4d
    assert plan.divisions == (44, 76, 32)
    for (lo, hi), n in zip(((0.0, geom.x_max_m), (geom.y_min_m, geom.y_max_m),
                            (geom.z_min_m, geom.z_max_m)), plan.divisions):
        assert (hi - lo) / n == pytest.approx(0.5)


def test_r4d_exit_and_plate_resolution(r4d):
    _, geom, _, plan, _, _ = r4d
    assert (plan.exit_level, plan.plate_level) == (4, 2)
    assert plan.exit_cell_size_m == pytest.approx(0.03125)
    assert plan.plate_cell_size_m == pytest.approx(0.125)
    assert geom.plate.thickness_m / plan.plate_cell_size_m == pytest.approx(2.0)


# --------------------------------------------------------------------------- #
# cell sizes
# --------------------------------------------------------------------------- #

def test_each_level_halves_the_cell(tmp_path):
    _, _, _, plan, _, _ = derived(tmp_path)
    assert [plan.cell_size_m(n) for n in range(4)] == [0.5, 0.25, 0.125, 0.0625]


def test_the_smallest_cell_is_the_finest_level_anywhere(tmp_path):
    """Region level 3, plate level 2: the region sets the time step."""
    _, _, _, plan, _, _ = derived(tmp_path)
    assert plan.max_level == 3
    assert plan.min_cell_size_m == pytest.approx(0.0625)

    _, _, _, finer_plate, _, _ = derived(
        tmp_path / "b", mesh={"plate_refinement_level": 4})
    assert finer_plate.min_cell_size_m == pytest.approx(0.03125)


def test_the_exit_cell_volume_is_a_cube(tmp_path):
    _, _, _, plan, _, _ = derived(tmp_path)
    assert plan.exit_cell_volume_m3 == pytest.approx(0.0625 ** 3)


# --------------------------------------------------------------------------- #
# the exit sits at one level
# --------------------------------------------------------------------------- #

def test_with_no_region_the_exit_is_at_the_background_level(tmp_path):
    cfg = load_case(tmp_path, mesh={"refinement_regions": [],
                                    "background_cell_size_m": 0.05},
                    domain={"x_max_m": 5.0, "y_half_m": 2.5, "z_half_m": 1.5},
                    plate={"thickness_m": 0.2})
    geom = geometry.from_config(cfg)
    assert mesh.exit_level([], geom, 0.05) == 0


def test_the_exit_takes_the_level_of_the_region_covering_it(tmp_path):
    _, _, _, plan, _, _ = derived(tmp_path)
    assert plan.exit_level == 3


def test_a_region_clear_of_the_exit_does_not_set_its_level(tmp_path):
    regions = [
        {"name": "exitCore", "level": 3, "x_max_m": 1.0, "radius_m": 0.5},
        {"name": "front", "level": 4, "shape": "box", "x_min_m": 3.0,
         "x_max_m": 3.9, "y_min_m": -1.0, "y_max_m": 1.0, "z_min_m": -0.5,
         "z_max_m": 0.5},
    ]
    _, _, _, plan, _, _ = derived(tmp_path, mesh={"refinement_regions": regions})
    assert plan.exit_level == 3
    assert plan.max_level == 4


def test_a_region_covering_only_part_of_the_exit_is_rejected(tmp_path):
    """The nozzle would be meshed at two cell sizes: the staircase would not be
    the one predicted, and 'the exit cell' would be right for only part of the
    inlet."""
    regions = [
        {"name": "exitCore", "level": 3, "x_max_m": 1.0, "radius_m": 0.5},
        {"name": "tight", "level": 4, "x_max_m": 1.0, "radius_m": 0.1},
    ]
    with pytest.raises(ValueError, match=r"levels \[3, 4\], not at one"):
        derived(tmp_path, mesh={"refinement_regions": regions})


def test_a_region_smaller_than_the_cells_around_it_is_never_seen(tmp_path):
    """snappyHexMesh refines by CELL CENTRE. The four 0.5 m cells round the axis
    have their centres 0.354 m from it, so a 0.3 m cylinder contains none of
    them and refines nothing -- whatever level it asks for. Only its buffer
    layers arrive, and they split the exit."""
    regions = [{"name": "small", "level": 3, "x_max_m": 1.0, "radius_m": 0.3}]
    with pytest.raises(ValueError, match="never seen"):
        derived(tmp_path, mesh={"refinement_regions": regions})


def test_a_region_just_reaching_the_surrounding_cell_centres_is_seen(tmp_path):
    """0.36 m is past the 0.354 m of the cells round the axis; from there every
    child over the disk has its centre inside too."""
    regions = [{"name": "enough", "level": 3, "x_max_m": 1.0, "radius_m": 0.36}]
    _, _, _, plan, _, _ = derived(tmp_path, mesh={"refinement_regions": regions})
    assert plan.exit_level == 3


def test_a_region_starting_past_the_first_cell_centres_does_not_refine_the_exit(
        tmp_path):
    """The first layer of 0.5 m cells has its centres at x = 0.25 m."""
    regions = [{"name": "late", "level": 3, "x_min_m": 0.3, "x_max_m": 1.0,
                "radius_m": 0.5}]
    with pytest.raises(ValueError, match="buffer layers of region 'late'"):
        derived(tmp_path, mesh={"refinement_regions": regions})


def test_a_plate_whose_buffer_layers_reach_the_exit_is_rejected(tmp_path):
    """A plate one metre from an unrefined exit refines part of it by itself."""
    with pytest.raises(ValueError, match="buffer layers of the plate"):
        derived(tmp_path, mesh={"refinement_regions": []},
                plate={"standoff_m": 1.0})


def test_the_exit_plane_cells_are_those_snappy_would_leave(tmp_path):
    """Level 3 over a 0.4 m disk: every 62.5 mm cell whose square reaches
    inside r = 0.2 m, at x = h/2."""
    cfg = load_case(tmp_path)
    geom = geometry.from_config(cfg)
    cells = mesh.exit_plane_cells(mesh.regions_from_config(cfg), geom, 0.5)
    assert {level for level, _ in cells} == {3}
    assert all(x == pytest.approx(0.03125) for _, (x, _, _) in cells)
    # Per quadrant, row by row out from the axis: 4 + 4 + 3 + 2 squares have
    # their nearest corner inside the circle.
    assert len(cells) == 4 * 13


def test_a_region_outside_the_domain_is_rejected(tmp_path):
    regions = [{"name": "exitCore", "level": 3, "x_max_m": 1.0, "radius_m": 0.5},
               {"name": "lost", "level": 2, "x_min_m": 7.0, "x_max_m": 8.0,
                "radius_m": 0.5}]
    with pytest.raises(ValueError, match="outside the domain"):
        derived(tmp_path, mesh={"refinement_regions": regions})


# --------------------------------------------------------------------------- #
# the nozzle patch
# --------------------------------------------------------------------------- #

def test_the_staircase_nozzle_is_four_fold_symmetric(tmp_path):
    """The axis is a lattice line, so the face count is a multiple of four."""
    _, _, _, _, nozzle, _ = derived(tmp_path)
    assert nozzle.n_faces % 4 == 0
    assert nozzle.n_faces > 0


def test_the_nozzle_area_is_whole_exit_faces(tmp_path):
    _, _, _, plan, nozzle, _ = derived(tmp_path)
    assert nozzle.area_m2 == pytest.approx(nozzle.n_faces * plan.exit_cell_size_m ** 2)


def test_the_nozzle_prediction_applies_the_face_centre_rule(tmp_path):
    """Counted here the slow way: every lattice face, centre inside R0."""
    _, geom, _, plan, nozzle, _ = derived(tmp_path)
    h = plan.exit_cell_size_m
    count = sum(
        1
        for j in range(-8, 8) for k in range(-8, 8)
        if ((j + 0.5) * h) ** 2 + ((k + 0.5) * h) ** 2 <= geom.radius_m ** 2)
    assert nozzle.n_faces == count


@pytest.mark.parametrize("planes,fraction", [(["y"], 2), (["z"], 2), (["y", "z"], 4)])
def test_a_symmetry_plane_keeps_its_share_of_the_nozzle_faces(tmp_path, planes,
                                                              fraction):
    _, _, _, _, full, _ = derived(tmp_path / "full")
    _, _, _, _, part, _ = derived(tmp_path / "part",
                                  domain={"symmetry_planes": planes})
    assert part.n_faces * fraction == full.n_faces
    assert part.area_error == pytest.approx(full.area_error)


def test_the_staircase_error_stays_inside_an_envelope_that_shrinks(tmp_path):
    """NOT monotone: counting lattice squares in a circle is the Gauss circle
    problem, and the error wanders in sign and size as the cell shrinks. For
    this 0.4 m exit it is -0.5% at 62.5 mm, -3.6% at half that and +1.8% at a
    quarter. What falls is the bound, in proportion to cell over radius."""
    def error(level: int) -> float:
        regions = [{"name": "exitCore", "level": level, "x_max_m": 1.0,
                    "radius_m": 0.6}]
        *_, nozzle, _ = derived(tmp_path / str(level),
                                mesh={"refinement_regions": regions})
        return nozzle.area_error

    errors = {level: error(level) for level in (3, 4, 5)}
    for level, value in errors.items():
        assert abs(value) < 0.5 * (0.5 / 2 ** level) / 0.2, level
    assert abs(errors[4]) > abs(errors[3])          # finer, and worse


@pytest.mark.parametrize("level,faces,error", [
    (4, 112, -0.01985),     # the shipped case
    (5, 460, +0.00640),
    (6, 1828, -0.00016),
])
def test_r4d_staircase_at_each_exit_level(tmp_path, level, faces, error):
    """What raising the exit region's level buys the R-4D nozzle: the injected
    mass flow is off by exactly this fraction."""
    import yaml

    data = yaml.safe_load((R4D_CASE / "case.yaml").read_text(encoding="utf-8"))
    data["mesh"]["refinement_regions"][0]["level"] = level
    data["mesh"]["max_cells"] = 50_000_000
    (tmp_path / "case.yaml").write_text(yaml.safe_dump(data), encoding="utf-8")

    *_, nozzle, _ = derive(load_case_config(tmp_path))
    assert nozzle.n_faces == faces
    assert nozzle.area_error == pytest.approx(error, abs=2e-5)


# --------------------------------------------------------------------------- #
# the cell estimate
# --------------------------------------------------------------------------- #

def test_with_nothing_to_refine_the_estimate_is_the_background(tmp_path):
    _, _, _, plan, _, _ = derived(
        tmp_path, plate=None,
        mesh={"refinement_regions": [], "background_cell_size_m": 0.05},
        domain={"x_max_m": 2.0, "y_half_m": 1.0, "z_half_m": 1.0})
    assert plan.cells_by_level == {0: pytest.approx(plan.n_background_cells)}
    assert not plan.needs_snappy


def test_a_region_adds_its_volume_at_its_level(tmp_path):
    """A lone cylinder of r = 0.5 m, 1 m long, at level 3 (62.5 mm)."""
    _, _, _, plan, _, _ = derived(tmp_path, plate=None)
    expected = math.pi * 0.5 ** 2 * 1.0 / 0.0625 ** 3
    assert plan.cells_by_level[3] == pytest.approx(expected, rel=0.03)


def test_buffer_layers_surround_a_region_at_each_coarser_level(tmp_path):
    """Levels change one at a time, so a level-3 region is wrapped in level 2
    and level 1 whether or not anything asked for them."""
    _, _, _, plan, _, _ = derived(tmp_path, plate=None)
    assert plan.cells_by_level[2] > 0
    assert plan.cells_by_level[1] > 0


def test_more_buffer_cells_make_a_larger_mesh(tmp_path):
    _, _, _, thin, _, _ = derived(tmp_path / "a", mesh={"n_cells_between_levels": 1})
    _, _, _, thick, _, _ = derived(tmp_path / "b", mesh={"n_cells_between_levels": 4})
    assert thick.estimated_cells > thin.estimated_cells


def test_a_symmetric_domain_is_estimated_at_its_share(tmp_path):
    _, _, _, full, _, _ = derived(tmp_path / "full")
    _, _, _, quarter, _, _ = derived(tmp_path / "quarter",
                                     domain={"symmetry_planes": ["y", "z"]})
    assert quarter.estimated_cells == pytest.approx(full.estimated_cells / 4, rel=0.03)


def test_an_estimate_over_the_ceiling_is_refused_with_a_breakdown(tmp_path):
    with pytest.raises(ValueError, match="mesh.max_cells") as excinfo:
        derived(tmp_path, mesh={"max_cells": 1000})
    assert "level 3" in str(excinfo.value)


def test_a_plate_too_thin_for_its_level_is_refused(tmp_path):
    """One cell through the thickness is a dented layer, not a body."""
    with pytest.raises(ValueError, match="plate_refinement_level"):
        derived(tmp_path, mesh={"plate_refinement_level": 1})


# --------------------------------------------------------------------------- #
# buffer layers
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("n_between,shell,skin", [
    (1, 1.5, 2), (2, 2.5, 2),               # hexRef8::consistentSlowRefinement
    (3, 2.5, 4), (4, 3.5, 4), (5, 4.5, 6),  # consistentSlowRefinement2
    (6, 5.5, 6), (7, 6.5, 8),
])
def test_buffer_layers_are_whole_cells_of_the_level_below(n_between, shell, skin):
    """Two rules, changing at three, each rounded to pairs of cells: a shell is
    n or n + 1 cells thick below three and n - 1 or n from there on; a surface
    keeps two cells a side up to two, then two more for every second step."""
    assert mesh.buffer_shell_cells(n_between) == shell
    assert mesh.surface_skin_cells(n_between) == skin


@pytest.mark.parametrize("n_between,skin_m", [(2, 0.25), (3, 0.5), (4, 0.5),
                                             (5, 0.75)])
def test_the_plate_s_skin_is_as_thick_as_its_cells_say(tmp_path, n_between, skin_m):
    """The minimal plate's upstream face is at x = 4 m and its cells are 0.125 m."""
    cfg = load_case(tmp_path)
    geom = geometry.from_config(cfg)
    level, inside, _, source = mesh.refinement_layers(geom, 0.5, [], 2, n_between)[0]
    assert (level, source) == (2, "the plate")
    assert bool(inside(4.0 - skin_m + 0.01, 0.0, 0.0))
    assert not bool(inside(4.0 - skin_m - 0.01, 0.0, 0.0))


def test_the_plate_s_skin_is_rounded_at_its_edges(tmp_path):
    """Within reach in a straight line, as a rule about distance makes it. The
    plate's edge is at y = 2 m and its skin is 0.5 m: a point 0.4 m ahead of the
    face AND 0.4 m beyond the edge is 0.57 m from the plate, and outside -- a
    box grown by 0.5 m would take it in."""
    cfg = load_case(tmp_path)
    geom = geometry.from_config(cfg)
    _, inside, bounds, _ = mesh.refinement_layers(geom, 0.5, [], 2, 3)[0]
    assert bool(inside(3.6, 2.0, 0.0))                  # 0.4 m, straight ahead
    assert bool(inside(4.1, 2.4, 0.0))                  # 0.4 m, straight out
    assert bool(inside(3.7, 2.3, 0.0))                  # 0.42 m, diagonally
    assert not bool(inside(3.6, 2.4, 0.0))              # 0.57 m, diagonally
    # The bounding box is still the square-cornered one that contains it.
    assert bounds == pytest.approx((3.5, 4.75, -2.5, 2.5, -1.5, 1.5))


def test_each_buffer_shell_adds_its_thickness_to_the_reach(tmp_path):
    """Level 2 skin 0.5 m, then one shell of level 1: 2.5 cells of 0.25 m."""
    cfg = load_case(tmp_path)
    geom = geometry.from_config(cfg)
    layers = mesh.refinement_layers(geom, 0.5, [], 2, 3)
    assert [level for level, *_ in layers] == [2, 1]
    _, shell, _, _ = layers[1]
    assert bool(shell(4.0 - 0.5 - 0.625 + 0.01, 0.0, 0.0))
    assert not bool(shell(4.0 - 0.5 - 0.625 - 0.01, 0.0, 0.0))


def test_three_cells_between_levels_is_not_three_cells():
    """The case's own setting: 2.5 on average, measured as two ahead of the
    plate and three behind it."""
    assert mesh.buffer_shell_cells(3) == pytest.approx((2 + 3) / 2)
    assert mesh.buffer_shell_cells(2) == mesh.buffer_shell_cells(3)


# --------------------------------------------------------------------------- #
# the seed point
# --------------------------------------------------------------------------- #

def test_the_seed_is_in_the_fluid_and_off_every_lattice_plane(tmp_path):
    _, geom, _, plan, _, _ = derived(tmp_path)
    x, y, z = plan.seed_point
    assert bool(geom.contains(x, y, z))
    assert not bool(geom.plate.contains(x, y, z, margin_m=plan.min_cell_size_m))
    # Not on a face of ANY level: its position in a finest cell is not dyadic.
    for value, low in ((x, 0.0), (y, geom.y_min_m), (z, geom.z_min_m)):
        fraction = ((value - low) / plan.min_cell_size_m) % 1.0
        assert 0.3 < fraction < 0.5


def test_the_seed_is_upstream_of_the_plate(tmp_path):
    _, geom, _, plan, _, _ = derived(tmp_path)
    assert plan.seed_point[0] < geom.plate.x_upstream_m


def test_the_seed_of_a_symmetric_domain_is_off_the_symmetry_planes(tmp_path):
    _, _, _, plan, _, _ = derived(tmp_path, domain={"symmetry_planes": ["y", "z"]})
    assert plan.seed_point[1] > 0.0 and plan.seed_point[2] > 0.0


# --------------------------------------------------------------------------- #
# dictionaries
# --------------------------------------------------------------------------- #

def test_block_mesh_dict_names_the_exit_plane_as_one_open_patch(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    text = mesh.render_block_mesh_dict(cfg, geom, plan)
    assert "hex (0 1 2 3 4 5 6 7) (12 12 8) simpleGrading (1 1 1)" in text
    upstream = text.split("upstreamVacuum")[1].split("}")[0]
    assert "type patch;" in upstream and "x_min" in upstream
    assert "type symmetry;" not in text


def test_block_mesh_dict_never_declares_a_wall(tmp_path):
    """The vacuum is 'patch'; the only wall comes from snappyHexMesh."""
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    assert "type wall;" not in mesh.render_block_mesh_dict(cfg, geom, plan)


def test_block_mesh_dict_gives_symmetry_planes_their_own_patch(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path,
                                       domain={"symmetry_planes": ["y", "z"]})
    text = mesh.render_block_mesh_dict(cfg, geom, plan)
    block = text.split("    symmetry\n")[1].split("}")[0]
    assert "type symmetry;" in block
    assert "y_min" in block and "z_min" in block
    outer = text.split("    vacuum\n")[1].split("}")[0]
    assert "y_min" not in outer and "z_min" not in outer and "y_max" in outer


def test_snappy_dict_carves_the_plate_as_a_wall(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    text = mesh.render_snappy_hex_mesh_dict(cfg, geom, plan)
    assert "type    searchableBox;" in text
    assert "min     (4 -2 -1);" in text and "max     (4.25 2 1);" in text
    surfaces = text.split("refinementSurfaces")[1].split("resolveFeatureAngle")[0]
    assert "level (2 2);" in surfaces and "type wall;" in surfaces


def test_snappy_dict_snaps_the_plate_to_its_edges(tmp_path):
    """Without it the plate's edges are chamfered and its face is 1.4% short."""
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    assert "implicitFeatureSnap     true;" in mesh.render_snappy_hex_mesh_dict(
        cfg, geom, plan)


def test_snappy_dict_carves_the_whole_plate_under_symmetry(tmp_path):
    """The domain cuts it at the plane; a half box would leave a face there."""
    cfg, geom, _, plan, _, _ = derived(tmp_path,
                                       domain={"symmetry_planes": ["y", "z"]})
    text = mesh.render_snappy_hex_mesh_dict(cfg, geom, plan)
    assert "min     (4 -2 -1);" in text


def test_snappy_dict_holds_each_region_at_its_level(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    text = mesh.render_snappy_hex_mesh_dict(cfg, geom, plan)
    assert "type    searchableCylinder;" in text
    assert "point1  (0 0 0);" in text and "point2  (1 0 0);" in text
    regions = text.split("refinementRegions")[1].split("locationInMesh")[0]
    assert "exitCore" in regions and "mode    inside;" in regions
    assert "levels  ((1e15 3));" in regions


def test_snappy_dict_seeds_the_plan_s_point(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    x, y, z = plan.seed_point
    assert (f"locationInMesh ({x:.10g} {y:.10g} {z:.10g});"
            in mesh.render_snappy_hex_mesh_dict(cfg, geom, plan))


def test_topo_set_selects_the_disk_out_of_the_exit_plane_only(tmp_path):
    """patchToFace BEFORE cylinderToFace, or internal faces are selected too."""
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    text = mesh.render_topo_set_dict(cfg, geom, plan)
    assert text.index("patchToFace") < text.index("cylinderToFace")
    assert "patch   upstreamVacuum;" in text
    assert "radius  0.2;" in text
    half = 0.25 * plan.exit_cell_size_m
    assert f"point1  ({-half:.10g} 0 0);" in text


def test_create_patch_makes_the_nozzle_a_patch_never_a_wall(tmp_path):
    cfg, *_ = derived(tmp_path)
    text = mesh.render_create_patch_dict(cfg)
    assert "name            nozzle;" in text
    assert "type        patch;" in text and "type        wall;" not in text
    assert "set             nozzleFaces;" in text


def test_write_mesh_setup_writes_the_whole_pipeline(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    written = [p.name for p in mesh.write_mesh_setup(tmp_path, cfg, geom, plan)]
    assert written == ["blockMeshDict", "snappyHexMeshDict", "meshQualityDict",
                       "topoSetDict", "createPatchDict"]
    for name in written:
        text = (tmp_path / "system" / name).read_text(encoding="utf-8")
        assert "\r" not in text and f"object      {name};" in text


def test_with_nothing_to_carve_a_stale_snappy_dict_is_removed(tmp_path):
    """./Allmesh runs snappyHexMesh if the file is there, so a leftover from an
    earlier configuration would carve a plate case.yaml no longer has."""
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    mesh.write_mesh_setup(tmp_path, cfg, geom, plan)
    assert (tmp_path / "system" / "snappyHexMeshDict").is_file()

    cfg, geom, _, plan, _, _ = derived(
        tmp_path, plate=None,
        mesh={"refinement_regions": [], "background_cell_size_m": 0.05},
        domain={"x_max_m": 2.0, "y_half_m": 1.0, "z_half_m": 1.0})
    written = [p.name for p in mesh.write_mesh_setup(tmp_path, cfg, geom, plan)]
    assert "snappyHexMeshDict" not in written
    assert not (tmp_path / "system" / "snappyHexMeshDict").exists()


def test_region_membership_grows_with_the_margin():
    region = mesh.Region(name="r", shape="cylinder", level=2, x_min_m=0.0,
                         x_max_m=1.0, radius_m=0.5)
    assert bool(region.contains(0.5, 0.3, 0.3))
    assert not bool(region.contains(0.5, 0.4, 0.4))
    assert bool(region.contains(0.5, 0.4, 0.4, grow_m=0.1))
    assert not bool(region.contains(1.05, 0.0, 0.0))
    assert bool(region.contains(1.05, 0.0, 0.0, grow_m=0.1))


def test_a_box_region_becomes_a_searchable_box():
    region = mesh.Region(name="b", shape="box", level=2, x_min_m=1.0, x_max_m=2.0,
                         y_min_m=-1.0, y_max_m=1.0, z_min_m=0.0, z_max_m=0.5)
    assert region.searchable().geometry_lines()[0].strip() == "type    searchableBox;"
    assert np.all(region.contains([1.5, 1.5], [0.0, 0.0], [0.25, 0.75])
                  == [True, False])
