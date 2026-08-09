"""The graded multi-block Cartesian mesh, the nozzle prediction, and the cost.

The mesh is where Cai's axisymmetric setup and this 3-D Cartesian one part
company, so these tests cover both halves of that: the cost of the literal
translation (which must be computed, not asserted), and the properties the
graded substitute has to preserve.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from _cai2012 import derived, load_case
from plumetools.cai2012 import mesh
from plumetools.cai2012.geometry import CaiGeometry, from_config
from plumetools.foamio.primitives import BOX_FACES


# --------------------------------------------------------------------------- #
# segments
# --------------------------------------------------------------------------- #

def test_uniform_segment_rounds_the_count_up():
    """Never coarser than asked for: a mesh silently 3% coarser than the DSMC
    criterion demanded is the kind of thing nobody notices."""
    segment = mesh.uniform_segment(0.0, 1.0, 0.3)
    assert segment.n_cells == 4
    assert segment.first_cell_m == pytest.approx(0.25)
    assert segment.first_cell_m <= 0.3


def test_uniform_segment_that_divides_exactly_is_not_over_refined():
    assert mesh.uniform_segment(0.0, 1.0, 0.25).n_cells == 4


def test_uniform_segment_edges_span_the_interval():
    edges = mesh.uniform_segment(0.0, 0.6, 0.1).cell_edges()
    assert edges[0] == pytest.approx(0.0)
    assert edges[-1] == pytest.approx(0.6)
    assert len(edges) == 7


def test_inverted_segment_is_rejected():
    with pytest.raises(ValueError, match="empty or inverted"):
        mesh.uniform_segment(1.0, 0.0, 0.1)


def test_geometric_segment_matches_the_core_cell_at_its_first_cell():
    """The expansion ratio is SOLVED FOR, not typed in. A step in cell size at
    the core boundary would show up in the sampled density as a feature that
    looks physical."""
    segment = mesh.geometric_segment(0.4, 2.0, 40, first_cell_m=0.01,
                                     max_expansion=100.0)
    assert segment.first_cell_m == pytest.approx(0.01, rel=1e-6)


def test_geometric_segment_widths_sum_to_the_length():
    segment = mesh.geometric_segment(0.3, 2.0, 30, 0.01, 100.0)
    edges = segment.cell_edges()
    assert edges[-1] == pytest.approx(2.0)
    assert np.all(np.diff(edges) > 0.0)


def test_geometric_segment_expands_monotonically():
    segment = mesh.geometric_segment(0.3, 2.0, 30, 0.01, 100.0)
    widths = np.diff(segment.cell_edges())
    assert np.all(np.diff(widths) > 0.0)
    assert segment.last_cell_m / segment.first_cell_m == pytest.approx(
        segment.expansion, rel=1e-6)


def test_geometric_segment_honours_the_expansion_cap():
    """When the cap binds, the first cell grows instead -- a visible jump, which
    MeshPlan.max_grading_jump reports."""
    segment = mesh.geometric_segment(0.3, 20.0, 5, 0.001, max_expansion=2.0)
    assert segment.expansion == pytest.approx(2.0)
    assert segment.first_cell_m > 0.001


def test_geometric_segment_falls_back_to_uniform_when_the_core_cell_is_ample():
    segment = mesh.geometric_segment(0.0, 1.0, 10, first_cell_m=0.5,
                                     max_expansion=40.0)
    assert segment.expansion == pytest.approx(1.0)
    assert segment.first_cell_m == pytest.approx(0.1)


def test_single_cell_segment_is_uniform():
    segment = mesh.geometric_segment(0.0, 1.0, 1, 0.1, 40.0)
    assert (segment.n_cells, segment.expansion) == (1, 1.0)
    assert segment.first_cell_m == pytest.approx(1.0)


def test_zero_cells_is_rejected():
    with pytest.raises(ValueError, match="n_cells"):
        mesh.geometric_segment(0.0, 1.0, 0, 0.1, 40.0)


# --------------------------------------------------------------------------- #
# the plan
# --------------------------------------------------------------------------- #

def test_plan_builds_eighteen_blocks_for_a_full_domain(tmp_path):
    """2 segments in x, 3 in each lateral direction."""
    _, _, _, plan, _ = derived(tmp_path)
    assert plan.n_blocks == 18
    assert (len(plan.x_segments), len(plan.y_segments), len(plan.z_segments)) == (
        2, 3, 3)


def test_half_domain_drops_the_negative_lateral_segment(tmp_path):
    _, _, _, plan, _ = derived(tmp_path, geometry={"symmetry_mode": "half_y"})
    assert len(plan.y_segments) == 2
    assert plan.n_blocks == 12


def test_the_core_cell_is_capped_by_geometry_at_high_knudsen(tmp_path):
    """At Kn = 100 lambda0 is 20 m -- a hundred nozzle diameters. The collision
    criterion alone would allow a cell larger than the nozzle, so
    max_core_cell_over_D is what actually sets the mesh."""
    _, geom, exit_state, plan, _ = derived(
        tmp_path, mesh={"core_cell_size_m": None})
    assert exit_state.mean_free_path_m == pytest.approx(20.0)
    assert plan.core_cell_size_m == pytest.approx(0.05 * 0.2)


def test_the_collision_criterion_wins_at_low_knudsen(tmp_path):
    """At Kn = 0.001 lambda0 = 0.2 mm, far below the geometric cap."""
    _, _, _, plan, _ = derived(
        tmp_path, exit={"knudsen": 0.001},
        mesh={"core_cell_size_m": None, "max_cells": 10 ** 12})
    assert plan.requested_cell_size_m == pytest.approx(0.0002)


def test_cell_budget_coarsens_and_says_so(tmp_path):
    """A mis-set Knudsen number must not silently request a billion cells."""
    _, _, _, plan, _ = derived(
        tmp_path, exit={"knudsen": 0.001},
        mesh={"core_cell_size_m": None, "max_cells": 200000})
    assert plan.coarsened is True
    assert plan.n_cells <= 200000
    assert plan.core_cell_size_m > plan.requested_cell_size_m
    assert "COARSENED" in "\n".join(plan.describe())


def test_an_impossible_budget_is_an_error_not_a_square_nozzle(tmp_path):
    with pytest.raises(ValueError, match="coarser than the nozzle radius"):
        derived(tmp_path, mesh={"core_cell_size_m": None, "max_cells": 8})


def test_cell_counts_multiply_out(tmp_path):
    _, _, _, plan, _ = derived(tmp_path)
    nx, ny, nz = plan.divisions
    assert plan.n_cells == nx * ny * nz


def test_the_smallest_cell_is_the_core_cell(tmp_path):
    _, _, _, plan, _ = derived(tmp_path)
    assert plan.min_cell_size_m == pytest.approx(
        min(plan.core_cell_sizes_m), rel=1e-9)


def test_the_largest_cell_is_out_at_the_vacuum_boundary(tmp_path):
    _, _, _, plan, _ = derived(tmp_path)
    assert plan.max_cell_size_m > plan.core_cell_size_m


def test_core_cell_sizes_are_the_realised_ones_not_the_requested(tmp_path):
    """Each axis rounds its count up, so the realised cells differ slightly and
    between axes. The particle weight is proportional to their product."""
    _, _, _, plan, _ = derived(tmp_path)
    dx, dy, dz = plan.core_cell_sizes_m
    assert plan.exit_cell_volume_m3 == pytest.approx(dx * dy * dz)
    for size in (dx, dy, dz):
        assert size <= plan.core_cell_size_m * (1.0 + 1e-9)


def test_cell_over_mean_free_path_is_reported(tmp_path):
    _, _, exit_state, plan, _ = derived(tmp_path, exit={"knudsen": 0.1})
    assert plan.cell_over_mfp == pytest.approx(
        plan.core_cell_size_m / exit_state.mean_free_path_m)


def test_grading_is_continuous_when_the_cap_does_not_bind(tmp_path):
    _, _, _, plan, _ = derived(tmp_path, mesh={"outer_expansion": 200.0})
    assert plan.max_grading_jump == pytest.approx(1.0, rel=1e-3)


def test_grading_jump_uses_the_core_side_cell_of_each_segment(tmp_path):
    """The low-side lateral segments expand towards decreasing y, so their small
    cell is the LAST one. Using first_cell_m would report the outer cell."""
    _, _, _, plan, _ = derived(tmp_path, mesh={"outer_expansion": 1.5})
    low = plan.y_segments[0]
    assert low.expansion < 1.0
    assert low.last_cell_m < low.first_cell_m
    assert plan.max_grading_jump == pytest.approx(
        max(s.min_cell_m for s in plan.x_segments + plan.y_segments
            + plan.z_segments if s.expansion != 1.0) / plan.core_cell_size_m)


def test_plan_is_serialisable(tmp_path):
    _, _, _, plan, _ = derived(tmp_path)
    document = plan.as_dict()
    assert document["n_cells"] == plan.n_cells
    assert document["coarsened_for_budget"] is False


# --------------------------------------------------------------------------- #
# the cost of Cai's own grid
# --------------------------------------------------------------------------- #

def test_uniform_cost_reproduces_the_four_billion_cell_figure(tmp_path):
    """dx = lambda0(Kn = 0.01) = 2 mm over 2 m x 4 m x 4 m."""
    cfg, geom, _, _, run = derived(tmp_path)
    cost = mesh.uniform_cost(cfg, geom, run.reference_cell_size_m)
    assert cost.cell_size_m == pytest.approx(0.002)
    assert cost.divisions == (1000, 2000, 2000)
    assert cost.n_cells == pytest.approx(4.0e9)


def test_uniform_cost_particles_use_the_occupancy_target(tmp_path):
    cfg, geom, _, _, run = derived(tmp_path)
    cost = mesh.uniform_cost(cfg, geom, run.reference_cell_size_m)
    assert cost.particles_at_target == pytest.approx(cost.n_cells * 20.0)


def test_uniform_cost_explains_why_it_is_not_built(tmp_path):
    cfg, geom, _, _, run = derived(tmp_path)
    lines = "\n".join(mesh.uniform_cost(cfg, geom, run.reference_cell_size_m)
                      .describe())
    assert "AXISYMMETRIC" in lines
    assert "NOT built" in lines


def test_uniform_cost_rejects_a_non_positive_cell(tmp_path):
    cfg, geom, _, _, _ = derived(tmp_path)
    with pytest.raises(ValueError, match="positive"):
        mesh.uniform_cost(cfg, geom, 0.0)


# --------------------------------------------------------------------------- #
# the nozzle patch
# --------------------------------------------------------------------------- #

def test_predicted_nozzle_area_approaches_pi_r_squared_as_the_cell_shrinks(tmp_path):
    """A circle has no exact Cartesian representation. The staircase converges."""
    errors = []
    for cell in (0.05, 0.02, 0.01, 0.005):
        _, geom, _, plan, _ = derived(tmp_path, mesh={"core_cell_size_m": cell})
        errors.append(abs(mesh.predicted_nozzle_faces(plan, geom).area_error))
    assert errors[-1] < errors[0]
    assert errors[-1] < 0.01


def test_predicted_nozzle_face_count_scales_as_the_area(tmp_path):
    """~pi/4 * (D/cell)^2 faces."""
    _, geom, _, plan, _ = derived(tmp_path, mesh={"core_cell_size_m": 0.01})
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    assert nozzle.n_faces == pytest.approx(math.pi / 4.0 * 400.0, rel=0.05)


def test_predicted_nozzle_area_is_the_sum_of_the_selected_face_areas(tmp_path):
    _, geom, _, plan, _ = derived(tmp_path, mesh={"core_cell_size_m": 0.01})
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    _, dy, dz = plan.core_cell_sizes_m
    assert nozzle.area_m2 == pytest.approx(nozzle.n_faces * dy * dz, rel=1e-9)


def test_nozzle_prediction_reports_cells_across_the_diameter(tmp_path):
    _, geom, _, plan, _ = derived(tmp_path, mesh={"core_cell_size_m": 0.01})
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    assert nozzle.cells_across_diameter == pytest.approx(20.0, rel=0.02)


def test_half_domain_nozzle_is_half_the_faces(tmp_path):
    _, geom_full, _, plan_full, _ = derived(
        tmp_path / "full", mesh={"core_cell_size_m": 0.01})
    _, geom_half, _, plan_half, _ = derived(
        tmp_path / "half", mesh={"core_cell_size_m": 0.01},
        geometry={"symmetry_mode": "half_y"})
    full = mesh.predicted_nozzle_faces(plan_full, geom_full)
    half = mesh.predicted_nozzle_faces(plan_half, geom_half)
    assert half.n_faces == pytest.approx(0.5 * full.n_faces, rel=0.1)


# --------------------------------------------------------------------------- #
# particle estimate
# --------------------------------------------------------------------------- #

def test_particle_estimate_scales_inversely_with_the_weight(tmp_path):
    _, geom, exit_state, plan, _ = derived(tmp_path)
    a = mesh.estimate_particles(plan, geom, exit_state, 1.0e9, samples_per_axis=8)
    b = mesh.estimate_particles(plan, geom, exit_state, 2.0e9, samples_per_axis=8)
    assert b["total_particles"] == pytest.approx(0.5 * a["total_particles"])


def test_particle_estimate_is_positive_and_finite(tmp_path):
    _, geom, exit_state, plan, run = derived(tmp_path)
    estimate = mesh.estimate_particles(
        plan, geom, exit_state, run.n_equivalent_particles, samples_per_axis=8)
    assert 0.0 < estimate["total_particles"] < float("inf")
    assert estimate["mean_particles_per_cell"] > 0.0


# --------------------------------------------------------------------------- #
# the dictionaries
# --------------------------------------------------------------------------- #

@pytest.fixture
def dicts(tmp_path):
    cfg, geom, _, plan, _ = derived(tmp_path)
    return cfg, geom, plan, mesh.render_block_mesh_dict(cfg, geom, plan)


def test_block_mesh_dict_has_one_block_per_lattice_cell(dicts):
    _, _, plan, text = dicts
    assert text.count("    hex (") == plan.n_blocks


def test_block_mesh_dict_vertex_count_is_the_lattice(dicts):
    _, _, plan, text = dicts
    body = text.split("vertices\n(", 1)[1].split(");", 1)[0]
    expected = ((len(plan.x_segments) + 1) * (len(plan.y_segments) + 1)
                * (len(plan.z_segments) + 1))
    assert len([line for line in body.splitlines() if line.strip()]) == expected


def test_block_mesh_dict_divisions_match_the_plan(dicts):
    _, _, plan, text = dicts
    counts = {s.n_cells for s in plan.x_segments}
    for n in counts:
        assert f"({n} " in text


def _uncommented(text: str) -> str:
    """The dictionary with its ``//`` comments stripped.

    The comments deliberately mention 'nozzle' and 'wall' to explain why neither
    is written, so an assertion about the dictionary's CONTENT has to look past
    them.
    """
    return "\n".join(line.split("//", 1)[0] for line in text.splitlines())


def test_x_zero_is_one_patch_the_nozzle_is_carved_later(dicts):
    """blockMesh can only name block faces, and a disk is not a block face."""
    cfg, _, _, text = dicts
    boundary = _uncommented(text).split("boundary\n(", 1)[1]
    assert cfg.mesh.patch_names["upstream_vacuum"] in boundary
    assert cfg.mesh.patch_names["nozzle"] not in boundary


def test_every_outer_boundary_is_type_patch_not_wall(dicts):
    """A wall would REFLECT: the plume would expand into a closed box."""
    _, _, _, text = dicts
    boundary = _uncommented(text).split("boundary\n(", 1)[1]
    assert "type wall" not in boundary
    assert boundary.count("type patch;") == 2       # upstreamVacuum + vacuum


def test_half_domain_writes_a_symmetry_patch(tmp_path):
    cfg, geom, _, plan, _ = derived(tmp_path, geometry={"symmetry_mode": "half_y"})
    text = mesh.render_block_mesh_dict(cfg, geom, plan)
    assert "type symmetry;" in text


def test_boundary_faces_cover_every_exterior_block_face(dicts):
    """Every face of the lattice's outer surface is claimed exactly once, or
    blockMesh's defaultFaces patch would silently absorb the rest."""
    _, _, plan, text = dicts
    ni, nj, nk = (len(plan.x_segments), len(plan.y_segments),
                  len(plan.z_segments))
    expected = 2 * (nj * nk + ni * nk + ni * nj)
    boundary = text.split("boundary\n(", 1)[1]
    quads = [line for line in boundary.splitlines()
             if line.strip().startswith("(") and line.strip().endswith(")")]
    assert len(quads) == expected


def test_boundary_faces_are_unique(dicts):
    _, _, _, text = dicts
    boundary = text.split("boundary\n(", 1)[1]
    quads = [line.strip() for line in boundary.splitlines()
             if line.strip().startswith("(") and line.strip().endswith(")")]
    assert len(set(quads)) == len(quads)


def test_topo_set_dict_restricts_to_the_exit_plane_before_the_radius(tmp_path):
    """Reversed, cylinderToFace would also select every internal face inside it."""
    cfg, geom, _, plan, _ = derived(tmp_path)
    text = mesh.render_topo_set_dict(cfg, geom, plan)
    assert text.index("patchToFace") < text.index("cylinderToFace")
    assert "upstreamVacuum" in text


def test_topo_set_selector_radius_is_the_nozzle_radius(tmp_path):
    cfg, geom, _, plan, _ = derived(tmp_path)
    text = mesh.render_topo_set_dict(cfg, geom, plan)
    assert "radius  0.1;" in text


def test_topo_set_selector_is_a_thin_disc_around_the_exit_plane(tmp_path):
    """A long cylinder could reach a boundary face somewhere else."""
    cfg, geom, _, plan, _ = derived(tmp_path)
    text = mesh.render_topo_set_dict(cfg, geom, plan)
    half = mesh.SELECTOR_HALF_THICKNESS_FRACTION * plan.min_cell_size_m
    assert f"point1  (-{half:.10g} 0 0);" in text
    assert f"point2  ({half:.10g} 0 0);" in text


def test_create_patch_dict_makes_the_nozzle_a_patch_not_a_wall(tmp_path):
    """plumeFieldInflow refuses to inject across a wallPolyPatch."""
    cfg, geom, _, _, _ = derived(tmp_path)
    text = _uncommented(mesh.render_create_patch_dict(cfg, geom))
    assert "name            nozzle;" in text
    assert "type        patch;" in text
    assert "wall" not in text
    assert "constructFrom   set;" in text


def test_write_mesh_setup_writes_the_four_dictionaries(tmp_path):
    cfg, geom, _, plan, _ = derived(tmp_path)
    paths = mesh.write_mesh_setup(tmp_path, cfg, geom, plan)
    assert [p.name for p in paths] == [
        "blockMeshDict", "topoSetDict", "createPatchDict", "meshQualityDict"]
    for path in paths:
        assert path.is_file() and path.read_text(encoding="utf-8").strip()


def test_written_dictionaries_use_lf_line_endings(tmp_path):
    """A dictionary whose line endings depend on the generating platform makes
    every diff unreadable."""
    cfg, geom, _, plan, _ = derived(tmp_path)
    for path in mesh.write_mesh_setup(tmp_path, cfg, geom, plan):
        assert b"\r\n" not in path.read_bytes()


def test_written_dictionaries_name_their_generator(tmp_path):
    cfg, geom, _, plan, _ = derived(tmp_path)
    for path in mesh.write_mesh_setup(tmp_path, cfg, geom, plan):
        assert mesh.GENERATOR in path.read_text(encoding="utf-8")
