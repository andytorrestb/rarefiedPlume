"""The 5-block hemispherical O-grid generator.

Runs without OpenFOAM: the dictionary is checked as a geometric object, not by
meshing it. `blockMesh` itself is exercised by the needs_openfoam tier.
"""

from __future__ import annotations

import math
import re
from dataclasses import replace

import numpy as np
import pytest

from plumetools.config import CaseConfig, MeshConfig
from plumetools.foamio.blockmesh import (
    BLOCKS,
    SPHERE_EDGES,
    arc_point,
    build_vertices,
    classify_boundary_faces,
    orient_block,
    render_block_mesh_dict,
    write_block_mesh_dict,
)

R, H, L = 0.5, 2.5, 5.0


@pytest.fixture
def cfg():
    return CaseConfig(mesh=MeshConfig(sphere_radius_m=R, box_half_width_m=H,
                                      box_length_m=L, n_tangential=20, n_radial=24))


@pytest.fixture
def verts():
    return build_vertices(R, H, L)


# --------------------------------------------------------------------------- #
# vertices
# --------------------------------------------------------------------------- #

def test_sixteen_vertices(verts):
    assert verts.shape == (16, 3)


def test_cube_and_equator_constants(verts):
    """a = R/sqrt(3) for the cap corners, b = R/sqrt(2) for the equator ring."""
    assert verts[0][0] == pytest.approx(R / math.sqrt(3))
    assert abs(verts[4][1]) == pytest.approx(R / math.sqrt(2))


def test_all_eight_inner_vertices_lie_on_the_sphere(verts):
    np.testing.assert_allclose(np.linalg.norm(verts[:8], axis=1), R, atol=1e-12)


def test_equator_ring_is_on_the_symmetry_plane(verts):
    np.testing.assert_allclose(verts[4:8, 0], 0.0, atol=1e-15)


def test_outer_vertices_are_on_the_box(verts):
    np.testing.assert_allclose(verts[8:12, 0], L)
    np.testing.assert_allclose(verts[12:16, 0], 0.0)
    assert np.all(np.abs(verts[8:16, 1:]) == H)


def test_geometry_is_parametric():
    other = build_vertices(1.0, 4.0, 9.0)
    np.testing.assert_allclose(np.linalg.norm(other[:8], axis=1), 1.0, atol=1e-12)
    np.testing.assert_allclose(other[8:12, 0], 9.0)
    assert np.all(np.abs(other[12:16, 1:]) == 4.0)


# --------------------------------------------------------------------------- #
# arcs
# --------------------------------------------------------------------------- #

def test_twelve_sphere_edges():
    assert len(SPHERE_EDGES) == 12
    assert len(set(map(frozenset, SPHERE_EDGES))) == 12


def test_every_arc_point_lies_on_the_sphere(verts):
    for v0, v1 in SPHERE_EDGES:
        p = arc_point(verts[v0], verts[v1], R)
        assert np.linalg.norm(p) == pytest.approx(R, abs=1e-12)


def test_arc_point_bisects_the_chord(verts):
    """R * (p1 + p2) / |p1 + p2| is equidistant from both endpoints."""
    for v0, v1 in SPHERE_EDGES:
        p = arc_point(verts[v0], verts[v1], R)
        assert np.linalg.norm(p - verts[v0]) == pytest.approx(
            np.linalg.norm(p - verts[v1]), rel=1e-12)


def test_equator_arc_midpoints_are_the_axis_crossings(verts):
    """The equator ring's arcs pass through (0, +-R, 0) and (0, 0, +-R)."""
    got = sorted(tuple(np.round(arc_point(verts[a], verts[b], R), 12))
                 for a, b in SPHERE_EDGES[8:])
    expected = sorted([(0.0, R, 0.0), (0.0, -R, 0.0), (0.0, 0.0, R), (0.0, 0.0, -R)])
    np.testing.assert_allclose(got, expected, atol=1e-12)


def test_antipodal_arc_is_rejected():
    with pytest.raises(ValueError, match="antipodal"):
        arc_point([R, 0, 0], [-R, 0, 0], R)


# --------------------------------------------------------------------------- #
# blocks
# --------------------------------------------------------------------------- #

def test_five_blocks():
    assert len(BLOCKS) == 5


def test_every_block_is_positively_oriented(verts):
    """OpenFOAM needs the first four vertices wound toward the second four.

    All five blocks require flipping relative to the naive vertex order, which is
    exactly why orientation is derived rather than hand-written -- a
    negative-volume block is not visible by inspection.
    """
    for name, inner_raw, outer_raw in BLOCKS:
        inner, outer = orient_block(inner_raw, outer_raw, verts)
        p = verts[list(inner)]
        n = np.cross(p[1] - p[0], p[2] - p[1])
        towards = verts[list(outer)].mean(axis=0) - p.mean(axis=0)
        assert np.dot(n, towards) > 0, f"block {name} is inverted"


def test_only_the_cap_block_avoids_the_symmetry_plane(verts):
    touching = [name for name, inner, _ in BLOCKS
                if np.any(np.abs(verts[list(inner)][:, 0]) < 1e-12)]
    assert touching == ["+y", "-y", "+z", "-z"]


# --------------------------------------------------------------------------- #
# boundary faces
# --------------------------------------------------------------------------- #

def test_face_counts(verts):
    inflow, sym, outer = classify_boundary_faces(verts, R, H, L)
    assert (len(inflow), len(sym), len(outer)) == (5, 4, 5)


def test_inflow_faces_tile_the_hemisphere(verts):
    """The five inner faces use exactly the eight sphere vertices."""
    inflow, _, _ = classify_boundary_faces(verts, R, H, L)
    assert sorted({i for f in inflow for i in f}) == list(range(8))


def test_inflow_normals_point_out_of_the_domain(verts):
    """The fluid is outside the cavity, so OpenFOAM's outward normal points in."""
    inflow, _, _ = classify_boundary_faces(verts, R, H, L)
    for f in inflow:
        p = verts[list(f)]
        n = np.cross(p[1] - p[0], p[2] - p[1])
        assert np.dot(n, -p.mean(axis=0)) > 0


def test_symmetry_faces_are_on_the_plane_and_point_minus_x(verts):
    _, sym, _ = classify_boundary_faces(verts, R, H, L)
    for f in sym:
        p = verts[list(f)]
        np.testing.assert_allclose(p[:, 0], 0.0, atol=1e-12)
        assert np.cross(p[1] - p[0], p[2] - p[1])[0] < 0


def test_outer_normals_point_away_from_the_box_centre(verts):
    _, _, outer = classify_boundary_faces(verts, R, H, L)
    centre = np.array([L / 2, 0.0, 0.0])
    for f in outer:
        p = verts[list(f)]
        n = np.cross(p[1] - p[0], p[2] - p[1])
        assert np.dot(n, p.mean(axis=0) - centre) > 0


# --------------------------------------------------------------------------- #
# rendered dictionary
# --------------------------------------------------------------------------- #

def test_dict_has_the_expected_sections(cfg):
    text = render_block_mesh_dict(cfg)
    for section in ("FoamFile", "vertices", "blocks", "edges", "boundary",
                    "mergePatchPairs"):
        assert section in text


def test_dict_declares_the_three_patches_with_the_right_types(cfg):
    text = render_block_mesh_dict(cfg)
    assert re.search(r"inflow\s*\{\s*type patch;", text)
    assert re.search(r"vacuum\s*\{\s*type patch;", text)
    assert re.search(r"sym\s*\{\s*type symmetry;", text)


def test_dict_never_emits_unspecified(cfg):
    """A generated mesh sidesteps HA-07 permanently."""
    assert "Unspecified" not in render_block_mesh_dict(cfg)


def test_dict_has_five_blocks(cfg):
    assert render_block_mesh_dict(cfg).count("hex (") == 5


# --------------------------------------------------------------------------- #
# curvature: arc edges vs projection onto a searchableSphere
# --------------------------------------------------------------------------- #

def test_arc_mode_emits_twelve_arcs_and_no_projection(cfg):
    text = render_block_mesh_dict(replace(cfg, mesh=replace(cfg.mesh, projection="arc")))
    # Count directive lines, not substrings: the header comment mentions "arc".
    assert len(re.findall(r"^\s*arc \d+ \d+ \(", text, re.M)) == 12
    assert "searchableSphere" not in text
    assert not re.search(r"^\s*project ", text, re.M)


def test_none_is_accepted_as_an_alias_for_arc(cfg):
    """Earlier configs wrote `projection: none`; they must keep working."""
    a = render_block_mesh_dict(replace(cfg, mesh=replace(cfg.mesh, projection="arc")))
    n = render_block_mesh_dict(replace(cfg, mesh=replace(cfg.mesh, projection="none")))
    assert a == n


def test_projection_mode_declares_the_sphere_primitive(cfg):
    text = render_block_mesh_dict(
        replace(cfg, mesh=replace(cfg.mesh, projection="searchable_sphere")))
    assert "geometry" in text
    assert "type    searchableSphere;" in text
    assert f"radius  {R}" in text
    # v1706 spells the centre key 'centre'; 'origin' would be silently ignored.
    assert "centre  (0 0 0);" in text


def test_projection_mode_projects_twelve_edges_and_five_faces(cfg):
    """Projecting edges alone is not enough.

    With arc edges the twelve block edges are exact but the face interiors remain
    ruled surfaces, so the patch dips to 16.3% inside the sphere at the cap face
    centre. The `faces` section is what makes it a true hemisphere.
    """
    text = render_block_mesh_dict(
        replace(cfg, mesh=replace(cfg.mesh, projection="searchable_sphere")))
    assert text.count("project ") == 17
    assert len(re.findall(r"^\s*project \d+ \d+ \(inflowSphere\)$", text, re.M)) == 12
    assert len(re.findall(r"^\s*project \([\d ]+\) inflowSphere$", text, re.M)) == 5
    assert "arc " not in text


def test_projected_faces_are_exactly_the_inflow_faces(cfg, verts):
    """The projected quads must be the same five the inflow patch declares."""
    text = render_block_mesh_dict(
        replace(cfg, mesh=replace(cfg.mesh, projection="searchable_sphere")))
    projected = {frozenset(int(i) for i in m.split())
                 for m in re.findall(r"^\s*project \(([\d ]+)\) inflowSphere$", text, re.M)}
    inflow, _, _ = classify_boundary_faces(verts, R, H, L)
    assert projected == {frozenset(f) for f in inflow}


def test_arc_mode_face_interiors_miss_the_sphere(verts):
    """Quantifies why projection exists, and pins the number in the docs.

    blockMesh fills an unprojected face by transfinite interpolation of its four
    edges; at the face centre that reduces to
    sum(edge midpoints)/2 - sum(corners)/4. The result is 16.3% inside the sphere
    on the cap face -- and it does NOT improve with n_tangential, because the face
    interior is determined by the edges alone.
    """
    inflow, _, _ = classify_boundary_faces(verts, R, H, L)
    worst = 0.0
    for quad in inflow:
        p = verts[list(quad)]
        mids = [arc_point(p[i], p[(i + 1) % 4], R) for i in range(4)]
        centre = np.sum(mids, axis=0) / 2.0 - np.sum(p, axis=0) / 4.0
        worst = max(worst, (R - np.linalg.norm(centre)) / R)
    assert worst == pytest.approx(0.1631, abs=1e-3)


def test_dict_resolution_follows_config(cfg):
    text = render_block_mesh_dict(cfg)
    assert "(20 20 24)" in text
    assert "2000 inflow quads" in text
    assert "48000 hexahedral cells" in text


def test_twenty_points_lie_on_the_sphere(cfg):
    """8 vertices + 12 arc interpolation points."""
    text = render_block_mesh_dict(cfg)
    coords = [tuple(map(float, m)) for m in
              re.findall(r"\(\s*(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*\)", text)]
    on_sphere = [p for p in coords if abs(np.linalg.norm(p) - R) < 1e-9]
    assert len(on_sphere) == 20


def test_dict_is_lf_terminated(tmp_path, cfg):
    (tmp_path / "system").mkdir()
    out = write_block_mesh_dict(tmp_path, cfg)
    assert b"\r\n" not in out.read_bytes()


def test_write_creates_system_directory(tmp_path, cfg):
    out = write_block_mesh_dict(tmp_path, cfg)
    assert out == tmp_path / "system" / "blockMeshDict"
    assert out.is_file()


# --------------------------------------------------------------------------- #
# rejected configurations
# --------------------------------------------------------------------------- #

def test_unknown_projection_is_rejected():
    cfg = CaseConfig(mesh=MeshConfig(projection="magic"))
    with pytest.raises(NotImplementedError, match="projection"):
        render_block_mesh_dict(cfg)


def test_unknown_mesh_type_is_rejected():
    cfg = CaseConfig(mesh=MeshConfig(type="something_else"))
    with pytest.raises(NotImplementedError, match="mesh.type"):
        render_block_mesh_dict(cfg)


def test_snappy_type_is_not_rendered_by_this_module():
    """block_mesh_ogrid and snappy_hex_sphere are separate generators; dispatch
    goes through plumetools.foamio.write_mesh_setup."""
    cfg = CaseConfig(mesh=MeshConfig(type="snappy_hex_sphere"))
    with pytest.raises(NotImplementedError):
        render_block_mesh_dict(cfg)


@pytest.mark.parametrize("kwargs, message", [
    ({"sphere_radius_m": 3.0}, "box_half_width_m"),
    ({"box_length_m": 0.25}, "box_length_m"),
    ({"n_tangential": 0}, "n_tangential"),
])
def test_impossible_geometry_is_rejected(kwargs, message):
    cfg = CaseConfig(mesh=MeshConfig(**kwargs))
    with pytest.raises(ValueError, match=message):
        render_block_mesh_dict(cfg)


# --------------------------------------------------------------------------- #
# outer patch type -- the standard-dsmcFoam workaround
# --------------------------------------------------------------------------- #

def test_outer_patch_defaults_to_patch(cfg):
    """Physically right: particles leave a plume domain through the outer boundary."""
    text = render_block_mesh_dict(cfg)
    assert re.search(r"vacuum\s*\{\s*type patch;", text)


def test_outer_patch_can_be_a_wall(cfg):
    """The only lever against FreeStream injecting on every patch-type boundary.

    FreeStream.C:57-62 appends every isType<polyPatch> patch with no selection
    list, and isType<> is an exact match -- so a wall is excluded. The cost is a
    reflecting, non-absorbing outer boundary; see docs/solver-compatibility.md.
    """
    text = render_block_mesh_dict(
        replace(cfg, mesh=replace(cfg.mesh, outer_patch_type="wall")))
    assert re.search(r"vacuum\s*\{\s*type wall;", text)
    assert "WALL: reflects" in text, "the physics cost must be stated in the dict"


def test_outer_patch_type_does_not_affect_inflow_or_sym(cfg):
    text = render_block_mesh_dict(
        replace(cfg, mesh=replace(cfg.mesh, outer_patch_type="wall")))
    assert re.search(r"inflow\s*\{\s*type patch;", text)
    assert re.search(r"sym\s*\{\s*type symmetry;", text)
