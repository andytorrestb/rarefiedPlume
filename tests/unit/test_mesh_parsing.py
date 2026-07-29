"""polyMesh parsing, against inputs that broke the pre-refactor readers."""

from __future__ import annotations

import numpy as np
import pytest

from _helpers import HEADER, SQUARE_POINTS, TRI_FACES, make_case, patch, write
from plumetools.mesh.boundary import read_boundary, require_patch
from plumetools.mesh.polymesh import read_faces, read_points
from plumetools.mesh.sets import read_face_set


def test_reads_points(tri_case):
    pts = read_points(tri_case)
    assert pts.shape == (6, 3)
    np.testing.assert_allclose(pts[2], [1, 1, 0])


def test_reads_triangular_faces(tri_case):
    assert read_faces(tri_case) == [[0, 1, 2], [0, 2, 3]]


def test_reads_quadrilateral_faces(quad_case):
    """blockMesh produces these; the legacy code assumed triangles throughout."""
    assert read_faces(quad_case) == [[0, 1, 2, 3], [1, 4, 5, 2]]


def test_reads_mixed_polygon_faces(mixed_case):
    """The archived 2d-wedge mesh mixes 3- and 4-vertex faces."""
    assert [len(f) for f in read_faces(mixed_case)] == [3, 3, 4]


def test_reads_boundary(tri_case):
    patches = read_boundary(tri_case)
    assert list(patches) == ["inflow"]
    assert patches["inflow"].n_faces == 2
    assert patches["inflow"].start_face == 0
    assert patches["inflow"].type == "patch"


def test_tolerates_ingroups_before_nfaces(ingroups_case):
    """blockMesh emits `inGroups`. The legacy reader skipped exactly three lines
    after the patch name, so any extra key shifted it onto the wrong line (RB-04)."""
    assert read_boundary(ingroups_case)["inflow"].n_faces == 2


def test_tolerates_unspecified_patch_type(unspecified_case):
    """Pointwise V18.5R2 writes `type Unspecified;` -- finding HA-07."""
    assert read_boundary(unspecified_case)["inflow"].type == "Unspecified"


def test_handles_crlf_line_endings(crlf_case):
    """The legacy nFaces parse was `line.split(' ', 14)[-1][0:-2]`, correct only
    for exactly one trailing LF. A CRLF checkout silently shifted it (RB-02)."""
    assert read_boundary(crlf_case)["inflow"].n_faces == 2
    assert read_points(crlf_case).shape == (6, 3)
    assert read_face_set(crlf_case, "inflow") == [0, 1]


def test_patch_name_that_prefixes_another(prefix_case):
    """`while patch not in line` matched `inflow` against `inflowOuter` (RB-03)."""
    patches = read_boundary(prefix_case)
    assert patches["inflow"].n_faces == 2
    assert patches["inflowOuter"].n_faces == 1


def test_extra_patch_does_not_shift_the_others(panel_case):
    patches = read_boundary(panel_case)
    assert list(patches) == ["inflow", "panel", "vacuum"]
    assert [p.start_face for p in patches.values()] == [0, 1, 2]


def test_face_set_order_is_preserved(tmp_path):
    """Face order is the contract with dsmcFoam+; it must not be sorted."""
    case = make_case(tmp_path / "order", points=SQUARE_POINTS, faces=TRI_FACES,
                     patches=[patch("inflow", 2, 0)], sets={"inflow": [1, 0]})
    assert read_face_set(case, "inflow") == [1, 0]


# --------------------------------------------------------------------------- #
# failure modes -- the legacy code hung or guessed instead
# --------------------------------------------------------------------------- #

def test_missing_patch_raises_promptly(tri_case):
    """The legacy `while patch not in line: readline()` spun forever at EOF."""
    patches = read_boundary(tri_case)
    with pytest.raises(KeyError, match="nope"):
        require_patch(patches, "nope")


def test_missing_patch_error_lists_available_patches(panel_case):
    patches = read_boundary(panel_case)
    with pytest.raises(KeyError, match="inflow"):
        require_patch(patches, "typo")


def test_missing_set_raises_instead_of_hanging(tri_case):
    """The error names the fix, since a missing set usually means topoSet was skipped."""
    with pytest.raises(FileNotFoundError, match="topoSet"):
        read_face_set(tri_case, "absent")


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_points(tmp_path / "nothing-here")


def test_truncated_list_raises_not_hangs(tmp_path):
    """A file whose declared count exceeds its contents must fail, not loop."""
    mesh = tmp_path / "bad" / "constant" / "polyMesh"
    write(mesh / "points",
          HEADER.format(cls="vectorField", obj="points") + "\n5\n(\n(0 0 0)\n(1 0 0)\n)\n")
    with pytest.raises((ValueError, IndexError)):
        read_points(tmp_path / "bad")


def test_parsers_terminate_on_garbage(tmp_path):
    """Bounded parsing: no input should make a reader spin (RB-01)."""
    mesh = tmp_path / "junk" / "constant" / "polyMesh"
    write(mesh / "points", "not an OpenFOAM file at all\n")
    with pytest.raises(ValueError, match="list block"):
        read_points(tmp_path / "junk")
