"""The snappyHexMesh setup: background box + searchableSphere primitive.

Runs without OpenFOAM -- the dictionaries are checked as text and as geometry.
Whether snappyHexMesh accepts them is a needs_openfoam question.
"""

from __future__ import annotations

import re
from dataclasses import replace

import numpy as np
import pytest

from plumetools.config import CaseConfig, MeshConfig
from plumetools.foamio import mesh_pipeline, write_mesh_setup
from plumetools.foamio.snappy import (
    background_divisions,
    location_in_mesh,
    render_background_block_mesh_dict,
    render_mesh_quality_dict,
    render_snappy_hex_mesh_dict,
    surface_cell_size,
    write_snappy_setup,
)

R, H, L = 0.5, 2.5, 5.0


@pytest.fixture
def cfg():
    return CaseConfig(mesh=MeshConfig(
        type="snappy_hex_sphere", sphere_radius_m=R, box_half_width_m=H,
        box_length_m=L, background_cell_size_m=0.125, refinement_level=3))


# --------------------------------------------------------------------------- #
# background mesh sizing
# --------------------------------------------------------------------------- #

def test_background_divisions_give_near_cubic_cells(cfg):
    nx, ny, nz = background_divisions(cfg)
    assert (nx, ny, nz) == (40, 40, 40)
    assert L / nx == pytest.approx(2 * H / ny)


def test_background_divisions_scale_with_cell_size(cfg):
    coarse = replace(cfg, mesh=replace(cfg.mesh, background_cell_size_m=0.25))
    assert background_divisions(coarse) == (20, 20, 20)


def test_cell_size_larger_than_the_sphere_is_rejected(cfg):
    """A cavity smaller than one background cell would fall between cells."""
    bad = replace(cfg, mesh=replace(cfg.mesh, background_cell_size_m=0.75))
    with pytest.raises(ValueError, match="exceeds the sphere radius"):
        background_divisions(bad)


def test_surface_resolution_comes_from_the_octree_not_the_background(cfg):
    """The refinement level, not the background size, sets the surface cell.

    A coarser background with one more level must land on the same surface size --
    that equivalence is what lets the far field be cheap without coarsening the
    plume source.
    """
    assert surface_cell_size(cfg) == pytest.approx(0.125 / 2 ** 3)
    coarser = replace(cfg, mesh=replace(
        cfg.mesh, background_cell_size_m=0.25, refinement_level=4))
    assert surface_cell_size(coarser) == pytest.approx(surface_cell_size(cfg))
    assert background_divisions(coarser) == (20, 20, 20)


def test_nonpositive_cell_size_is_rejected(cfg):
    bad = replace(cfg, mesh=replace(cfg.mesh, background_cell_size_m=0.0))
    with pytest.raises(ValueError, match="must be positive"):
        background_divisions(bad)


# --------------------------------------------------------------------------- #
# locationInMesh -- the classic way to get a wrong region
# --------------------------------------------------------------------------- #

def test_location_in_mesh_is_in_the_fluid(cfg):
    p = np.array(location_in_mesh(cfg))
    assert np.linalg.norm(p) > R, "must be outside the cavity"
    assert 0 < p[0] < L
    assert abs(p[1]) < H and abs(p[2]) < H


def test_location_in_mesh_avoids_every_symmetry_plane(cfg):
    """On a plane of symmetry, a face, or an edge, snappyHexMesh can keep the
    wrong region or fail outright."""
    x, y, z = location_in_mesh(cfg)
    assert y != 0.0 and z != 0.0
    assert abs(y) != abs(z)
    for v, limit in ((x, L), (y, H), (z, H)):
        assert v not in (0.0, limit, -limit)


@pytest.mark.parametrize("size", [0.5, 0.25, 0.2, 0.125, 0.1, 0.0625])
def test_location_in_mesh_never_lands_on_a_background_cell_face(cfg, size):
    """The failure mode being guarded: a fixed fraction of the domain is on a cell
    face whenever it happens to be a multiple of the cell size. ``0.5 * L`` is,
    for every even ``nx`` -- which the default 40 is.

    Snapping to the containing cell's centre makes it interior for any size, so
    the seed tracks the grid rather than drifting onto it.
    """
    sized = replace(cfg, mesh=replace(cfg.mesh, background_cell_size_m=size))
    nx, ny, nz = background_divisions(sized)
    point = location_in_mesh(sized)

    for value, lo, span, n in ((point[0], 0.0, L, nx),
                               (point[1], -H, 2 * H, ny),
                               (point[2], -H, 2 * H, nz)):
        h = span / n
        offset = ((value - lo) / h) % 1.0
        assert offset == pytest.approx(0.5), (
            f"{value} sits {offset:.3f} of the way through its cell -- "
            f"0.0 would be a cell face")


def test_location_in_mesh_stays_outside_the_cavity_when_snapped(cfg):
    """Snapping to a cell centre moves the seed by up to half a cell; it must not
    move it into the region snappyHexMesh is about to discard."""
    for size in (0.5, 0.25, 0.125, 0.0625):
        sized = replace(cfg, mesh=replace(cfg.mesh, background_cell_size_m=size))
        assert np.linalg.norm(location_in_mesh(sized)) > R


# --------------------------------------------------------------------------- #
# background blockMeshDict
# --------------------------------------------------------------------------- #

def test_background_is_a_single_plain_box(cfg):
    text = render_background_block_mesh_dict(cfg)
    assert text.count("hex (") == 1
    assert "(40 40 40)" in text
    # The sphere is introduced by snappyHexMesh, not by curved edges here.
    assert "arc " not in text
    assert "searchableSphere" not in text


def test_background_declares_only_sym_and_vacuum(cfg):
    """`inflow` must NOT exist yet -- snappyHexMesh creates it when it carves."""
    text = render_background_block_mesh_dict(cfg)
    assert re.search(r"sym\s*\{\s*type symmetry;", text)
    assert re.search(r"vacuum\s*\{\s*type patch;", text)
    assert not re.search(r"^\s*inflow\s*$", text, re.M)


def test_background_outer_patch_honours_outer_patch_type(cfg):
    """Not cosmetic. Standard dsmcFoam's FreeStream injects on every patch-type
    boundary, so an outer `patch` becomes a second inflow and the run aborts with
    "Zero boundary temperature detected". The O-grid respects
    `outer_patch_type`; the background box has to as well, or switching mesh.type
    silently reintroduces the abort."""
    walled = replace(cfg, mesh=replace(cfg.mesh, outer_patch_type="wall"))
    text = render_background_block_mesh_dict(walled)
    assert re.search(r"vacuum\s*\{\s*type wall;", text)
    assert not re.search(r"vacuum\s*\{\s*type patch;", text)


def test_background_box_spans_the_configured_domain(cfg):
    text = render_background_block_mesh_dict(cfg)
    # Scope to the vertices section: the blocks section also holds triples,
    # for the cell divisions and the grading.
    block = text[text.index("vertices"):text.index("blocks")]
    pts = np.array([[float(v) for v in m] for m in re.findall(
        r"\(\s*(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*\)", block)])
    assert len(pts) == 8, "a box has eight corners"
    assert pts[:, 0].min() == 0.0 and pts[:, 0].max() == L
    assert pts[:, 1].min() == -H and pts[:, 1].max() == H
    assert pts[:, 2].min() == -H and pts[:, 2].max() == H


def test_background_has_six_boundary_faces(cfg):
    """One sym face plus five vacuum faces -- a closed box."""
    text = render_background_block_mesh_dict(cfg)
    body = text[text.index("boundary"):text.index("mergePatchPairs")]
    assert len(re.findall(r"^\s*\(\d+ \d+ \d+ \d+\)\s*(?://.*)?$", body, re.M)) == 6


# --------------------------------------------------------------------------- #
# snappyHexMeshDict
# --------------------------------------------------------------------------- #

def test_geometry_is_the_sphere_primitive_named_after_the_patch(cfg):
    """snappyHexMesh names the created patch after the geometry entry, so the
    entry must be called `inflow`."""
    text = render_snappy_hex_mesh_dict(cfg)
    assert re.search(r"geometry\s*\{\s*inflow\s*\{\s*type\s+searchableSphere;", text)
    assert f"radius  {R}" in text
    assert "centre  (0 0 0);" in text


def test_no_stl_is_referenced(cfg):
    """The point of a primitive: snapping targets the true sphere, not a
    triangulation, and nothing outside the repository is needed."""
    text = render_snappy_hex_mesh_dict(cfg)
    assert ".stl" not in text
    assert "triSurfaceMesh" not in text


def test_castellate_and_snap_are_on_layers_are_off(cfg):
    text = render_snappy_hex_mesh_dict(cfg)
    assert "castellatedMesh true;" in text
    assert "snap            true;" in text
    assert "addLayers       false;" in text   # DSMC needs no boundary-layer stack


def test_refinement_surface_uses_the_configured_level(cfg):
    text = render_snappy_hex_mesh_dict(cfg)
    assert "level (3 3);" in text
    finer = render_snappy_hex_mesh_dict(
        replace(cfg, mesh=replace(cfg.mesh, refinement_level=5)))
    assert "level (5 5);" in finer


def test_inflow_patch_type_is_patch_not_wall(cfg):
    """A `wall` patch would make the plume source a no-slip surface."""
    text = render_snappy_hex_mesh_dict(cfg)
    assert re.search(r"patchInfo\s*\{\s*type patch;", text)
    assert "type wall;" not in text


def test_location_in_mesh_is_written(cfg):
    text = render_snappy_hex_mesh_dict(cfg)
    m = re.search(r"locationInMesh \(([-\d.eE+ ]+)\);", text)
    assert m
    p = np.array([float(v) for v in m.group(1).split()])
    np.testing.assert_allclose(p, location_in_mesh(cfg), rtol=1e-6)
    assert np.linalg.norm(p) > R


def test_mesh_quality_dict_is_included_and_written(cfg, tmp_path):
    assert '#include "meshQualityDict"' in render_snappy_hex_mesh_dict(cfg)
    names = [p.name for p in write_snappy_setup(tmp_path, cfg)]
    assert "meshQualityDict" in names, "the #include would dangle"


def test_negative_refinement_level_is_rejected(cfg):
    bad = replace(cfg, mesh=replace(cfg.mesh, refinement_level=-1))
    with pytest.raises(ValueError, match="refinement_level"):
        render_snappy_hex_mesh_dict(bad)


def test_mesh_quality_dict_has_the_expected_keys(cfg):
    text = render_mesh_quality_dict(cfg)
    for key in ("maxNonOrtho", "minVol", "minTetQuality", "minDeterminant", "relaxed"):
        assert key in text


# --------------------------------------------------------------------------- #
# writing and dispatch
# --------------------------------------------------------------------------- #

def test_writes_three_dictionaries(cfg, tmp_path):
    paths = write_snappy_setup(tmp_path, cfg)
    assert [p.name for p in paths] == ["blockMeshDict", "snappyHexMeshDict", "meshQualityDict"]
    assert all(p.is_file() and p.stat().st_size > 0 for p in paths)
    assert all(p.parent.name == "system" for p in paths)


def test_written_dictionaries_are_lf_terminated(cfg, tmp_path):
    for p in write_snappy_setup(tmp_path, cfg):
        assert b"\r\n" not in p.read_bytes(), p.name


def test_write_snappy_setup_rejects_the_wrong_mesh_type(tmp_path):
    ogrid = CaseConfig(mesh=MeshConfig(type="block_mesh_ogrid"))
    with pytest.raises(ValueError, match="snappy_hex_sphere"):
        write_snappy_setup(tmp_path, ogrid)


def test_dispatch_selects_the_right_generator(cfg, tmp_path):
    assert [p.name for p in write_mesh_setup(tmp_path, cfg)] == [
        "blockMeshDict", "snappyHexMeshDict", "meshQualityDict"]

    ogrid = CaseConfig(mesh=MeshConfig(type="block_mesh_ogrid"))
    assert [p.name for p in write_mesh_setup(tmp_path / "o", ogrid)] == ["blockMeshDict"]


def test_pipeline_matches_the_mesh_type(cfg):
    assert mesh_pipeline(cfg) == ["blockMesh", "snappyHexMesh"]
    assert mesh_pipeline(CaseConfig(mesh=MeshConfig())) == ["blockMesh"]


def test_dispatch_rejects_an_unknown_type(tmp_path):
    bad = CaseConfig(mesh=MeshConfig(type="nope"))
    with pytest.raises(NotImplementedError, match="mesh.type"):
        write_mesh_setup(tmp_path, bad)
    with pytest.raises(NotImplementedError, match="mesh.type"):
        mesh_pipeline(bad)
