"""Streamed reads give byte-for-byte what the whole-file versions gave.

One path stopped holding a whole file in memory on 2026-10-27:

* face ranges and face-set labels of ``constant/polyMesh/faces``
  (:mod:`plumetools.mesh.polymesh`), used for every patch read.

It is checked here against the code it replaced.
"""

from __future__ import annotations

import random

import pytest

from plumetools.mesh import polymesh


def write_faces(case, faces, *, header_comment=True):
    poly = case / "constant" / "polyMesh"
    poly.mkdir(parents=True, exist_ok=True)
    lines = ["/*--------------------------------*- C++ -*----------------------------------*\\",
             "  note: a comment ( with a paren",
             "\\*---------------------------------------------------------------------------*/"
             ] if header_comment else []
    lines += ["FoamFile", "{", "    format      ascii;", "    class       faceList;",
              "    object      faces;", "}", "", "", str(len(faces)), "("]
    lines += [f"{len(f)}({' '.join(map(str, f))})" for f in faces]
    lines += [")", "", "// ****** //"]
    (poly / "faces").write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def mesh(tmp_path):
    rng = random.Random(4)
    faces = [[rng.randrange(1000) for _ in range(rng.choice((3, 4, 4, 5, 6)))]
             for _ in range(500)]
    write_faces(tmp_path, faces)
    return tmp_path, faces


@pytest.mark.parametrize("start,n", [(0, 1), (0, 500), (17, 40), (499, 1), (250, 0)])
def test_a_face_range_is_what_the_whole_file_read_gave(mesh, start, n):
    case, faces = mesh
    assert polymesh.read_face_range(case, start, n) == faces[start:start + n]
    assert polymesh.read_face_range(case, start, n) == \
        polymesh.read_faces(case)[start:start + n]


def test_several_ranges_come_from_one_pass_in_the_order_asked(mesh):
    case, faces = mesh
    ranges = [(400, 50), (10, 5), (12, 20), (0, 0), (450, 50)]   # overlapping too
    got = polymesh.read_face_ranges(case, ranges)
    assert got == [faces[s:s + n] for s, n in ranges]


def test_face_set_labels_keep_their_order(mesh):
    case, faces = mesh
    labels = [499, 3, 3, 250, 0, 17]
    assert polymesh.read_faces_at(case, labels) == [faces[i] for i in labels]


@pytest.mark.parametrize("call,message", [
    (lambda c: polymesh.read_face_range(c, 490, 20), "faces 490 .. 510 requested"),
    (lambda c: polymesh.read_faces_at(c, [500]), "labels 500 .. 500 requested"),
])
def test_out_of_range_is_refused(mesh, call, message):
    with pytest.raises(ValueError, match=message):
        call(mesh[0])


def test_a_malformed_or_miscounted_face_is_refused(tmp_path):
    write_faces(tmp_path, [[0, 1, 2], [1, 2, 3]])
    path = tmp_path / "constant" / "polyMesh" / "faces"
    text = path.read_text().replace("3(1 2 3)", "4(1 2 3)")
    path.write_text(text)
    with pytest.raises(ValueError, match="face 1: declared 4 vertices, found 3"):
        polymesh.read_face_range(tmp_path, 0, 2)
    path.write_text(text.replace("4(1 2 3)", "garbage"))
    with pytest.raises(ValueError, match="face 1: expected 'N"):
        polymesh.read_face_range(tmp_path, 1, 1)


def test_a_list_shorter_than_its_count_is_refused(tmp_path):
    write_faces(tmp_path, [[0, 1, 2], [1, 2, 3]])
    path = tmp_path / "constant" / "polyMesh" / "faces"
    text = path.read_text().replace("\n2\n(", "\n3\n(")
    path.write_text(text)
    with pytest.raises(ValueError, match=r"face 2: expected 'N\(a b c \.\.\.\)', got '\)'"):
        polymesh.read_face_range(tmp_path, 0, 3)
    path.write_text(text[:text.index("3(1 2 3)") + len("3(1 2 3)\n")])   # cut off
    with pytest.raises(ValueError, match="ends before face 2"):
        polymesh.read_face_range(tmp_path, 0, 3)
