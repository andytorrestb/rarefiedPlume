"""Streamed reads and writes give byte-for-byte what the whole-file versions gave.

Three paths stopped holding a whole file in memory on 2026-10-27:

* face ranges and face-set labels of ``constant/polyMesh/faces``
  (:mod:`plumetools.mesh.polymesh`), used for every patch read;
* the per-cell values of a rarefaction diagnostic field
  (:func:`plumetools.rarefaction.foam._write_streamed`);
* the log summaries ``interface_case collect`` brings back from a cluster.

Each is checked here against the code it replaced.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from plumetools.cai2012ext import interface_case
from plumetools.foamio.primitives import write_lines
from plumetools.mesh import polymesh
from plumetools.rarefaction import foam


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


# --------------------------------------------------------------------------- #
# the rarefaction field writer
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("n", [0, 1, 6, 7, 8, 20])
def test_streamed_field_is_byte_identical_to_the_list_of_lines(tmp_path, monkeypatch, n):
    monkeypatch.setattr(foam, "WRITE_CHUNK", 7)
    values = np.random.default_rng(n).lognormal(size=n) * 1e19
    head = ["// head", "internalField   nonuniform List<scalar>", str(n), "("]
    tail = [")", ";", "", "// tail"]
    foam._write_streamed(tmp_path / "new", head, values, tail)
    write_lines(tmp_path / "old", head + [f"{v:.10g}" for v in values.tolist()] + tail)
    assert (tmp_path / "new").read_bytes() == (tmp_path / "old").read_bytes()


# --------------------------------------------------------------------------- #
# the log summary
# --------------------------------------------------------------------------- #

def legacy_summary(path):
    """interface_case.collect's per-log summary as it was at the tag."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return {"lines": len(lines), "tail": lines[-15:],
            "memory_lines": [x for x in lines if "MB" in x or "memory" in x.lower()][-5:],
            "execution_time": [x for x in lines if x.startswith("ExecutionTime")][-1:]}


@pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"])
def test_log_summary_matches_reading_the_whole_file(tmp_path, ending):
    body = []
    for step in range(60):
        body += [f"Time = {step}", f"    Memory {step} MB", "",
                 f"ExecutionTime = {step}.5 s  ClockTime = {step} s"]
    text = ending.join(body) + ending + "form\x0cfeed" + ending + "MeMoRy" + ending + "no newline"
    path = tmp_path / "log.dsmcFoam"
    path.write_bytes(text.encode("utf-8") + b"\xff bad byte")
    assert interface_case.log_summary(path) == legacy_summary(path)


def test_an_empty_log_summarises_as_empty(tmp_path):
    path = tmp_path / "log.dsmcFoam"
    path.write_text("")
    assert interface_case.log_summary(path) == legacy_summary(path)
