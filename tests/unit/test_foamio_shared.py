"""The shared OpenFOAM readers: plumetools.foamio.{volfield,times,patchfields}.

One field parser replaced three -- ``nozzleplume.post``, ``cai2012.post`` and
``markelov1999.foamfields`` -- so the corpus below is the union of what those
read, and each legacy parser is kept here, verbatim from the tag
``pre-refactor-2026-10-27``, as an oracle: on every input it accepted, the
shared reader must return the identical array.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from plumetools.foamio import patchfields, times, volfield
from plumetools.foamio.errors import FoamDataError


# --------------------------------------------------------------------------- #
# the legacy parsers, verbatim, as oracles
# --------------------------------------------------------------------------- #

def _legacy_nozzleplume_values(text, start, what):
    head = re.compile(r"\s*(uniform|nonuniform)\s*").match(text, start)
    if not head:
        raise ValueError(f"{what}: expected 'uniform' or 'nonuniform'")
    semicolon = text.find(";", head.end())
    if semicolon < 0:
        raise ValueError(f"{what}: unterminated entry")
    if head.group(1) == "uniform":
        token = text[head.end():semicolon].strip()
        numbers = np.asarray(
            [float(v) for v in token.strip("()").split()], dtype=np.float64)
        return numbers.reshape(1, 3) if token.startswith("(") else numbers
    kind = re.compile(r"List<(\w+)>").match(text, head.end())
    if not kind:
        raise ValueError(f"{what}: nonuniform value with no List<> type")
    opening = text.find("(", kind.end(), semicolon)
    if opening < 0:
        return np.empty((0, 3) if kind.group(1) == "vector" else (0,))
    closing = text.rfind(")", opening, semicolon)
    payload = text[opening + 1:closing].replace("(", " ").replace(")", " ")
    numbers = np.fromstring(payload, dtype=np.float64, sep=" ")
    return numbers.reshape(-1, 3) if kind.group(1) == "vector" else numbers


def legacy_nozzleplume_internal(text):
    match = re.search(r"^internalField", text, re.M)
    if not match:
        raise ValueError("no internalField entry")
    return _legacy_nozzleplume_values(text, match.end(), "internalField")


def legacy_nozzleplume_patch(text, patch):
    start = text.find("boundaryField")
    entry = re.compile(rf"^\s*{re.escape(patch)}\s*$", re.M).search(text, start)
    if not entry:
        raise ValueError("not in the boundaryField")
    body_start = text.index("{", entry.end())
    body_end = text.index("}", body_start)
    value = re.compile(r"^\s*value\s", re.M).search(text, body_start, body_end)
    if not value:
        raise ValueError("no 'value' entry")
    return _legacy_nozzleplume_values(text, value.end() - 1, "patch")


def legacy_cai2012_internal(text):
    head = re.compile(r"internalField\s+(uniform|nonuniform)", re.S).search(text)
    if not head:
        raise ValueError("no internalField entry")
    if head.group(1) == "uniform":
        token = re.search(r"internalField\s+uniform\s+(\([^)]*\)|[-+0-9.eE]+)\s*;", text)
        if not token:
            raise ValueError("malformed uniform internalField")
        value = token.group(1)
        if value.startswith("("):
            return np.asarray([[float(v) for v in re.findall(r"[-+0-9.eE]+", value)]],
                              dtype=np.float64)
        return np.asarray([float(value)], dtype=np.float64)
    kind = re.search(r"internalField\s+nonuniform\s+List<(\w+)>", text)
    if not kind:
        raise ValueError("nonuniform internalField with no List<> type")
    start = text.index("(", head.start())
    depth, i = 0, start
    while i < len(text):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                break
        i += 1
    payload = text[start + 1:i]
    values = np.fromstring(payload.replace("(", " ").replace(")", " "), sep=" ")
    return values.reshape(-1, 3) if kind.group(1) == "vector" else values


def legacy_markelov_patch(text, patch):
    start = text.find("boundaryField")
    match = re.search(rf"^\s*{re.escape(patch)}\s*$", text[start:], re.M)
    if not match:
        raise ValueError("not in the boundaryField")
    body_start = text.index("{", start + match.end())
    depth, i = 0, body_start
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                break
        i += 1
    body = text[body_start + 1:i]
    nonuniform = re.search(
        r"value\s+nonuniform\s+List<(\w+)>\s*(\d+)?\s*\((.*?)\)\s*;", body, re.S)
    if nonuniform:
        kind, _count, payload = nonuniform.groups()
        numbers = [float(v) for v in re.findall(r"[-+0-9.eE]+", payload)]
        arr = np.asarray(numbers, dtype=np.float64)
        return arr.reshape(-1, 3) if kind == "vector" else arr
    uniform = re.search(r"value\s+uniform\s+(\([^)]*\)|[-+0-9.eE]+)\s*;", body)
    if uniform:
        token = uniform.group(1)
        if token.startswith("("):
            return np.asarray(
                [[float(v) for v in re.findall(r"[-+0-9.eE]+", token)]],
                dtype=np.float64)
        return np.asarray([float(token)], dtype=np.float64)
    raise ValueError("no 'value' entry")


# --------------------------------------------------------------------------- #
# the corpus
# --------------------------------------------------------------------------- #

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       %s;
    location    "0.5";
    object      %s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 -3 0 0 0 0 0];

"""

FOOTER = "\n// ************************************************************************* //\n"


def field(internal: str, boundary: str = "", cls="volScalarField") -> str:
    return (HEADER % (cls, "f") + "internalField   " + internal + "\n\nboundaryField\n{\n"
            + boundary + "}\n" + FOOTER)


SCALARS = [1.5e20, 0.0, -2.25, 3.125e-7, 1.0e-300]
VECTORS = [(1.0, 2.0, 3.0), (-4.5, 0.0, 6.25e3), (7.0, 8.0, -9.0e-5)]

INTERNAL_CASES = {
    "uniform scalar": (field("uniform 2.5;"), np.array([2.5])),
    "uniform scalar exponent": (field("uniform 1.25e+20;"), np.array([1.25e20])),
    "uniform vector": (field("uniform (1 -2 3.5);", cls="volVectorField"),
                       np.array([[1.0, -2.0, 3.5]])),
    "nonuniform scalar": (field(
        "nonuniform List<scalar>\n5\n(\n" + "\n".join(f"{v:.10g}" for v in SCALARS)
        + "\n)\n;"), np.array(SCALARS)),
    "nonuniform vector": (field(
        "nonuniform List<vector>\n3\n(\n" + "\n".join(f"({x:.10g} {y:.10g} {z:.10g})"
                                                     for x, y, z in VECTORS)
        + "\n)\n;", cls="volVectorField"), np.array(VECTORS)),
    "inline scalar": (field("nonuniform List<scalar> 3(1 2 3);"), np.array([1.0, 2.0, 3.0])),
    "inline vector": (field("nonuniform List<vector> 2((1 2 3) (4 5 6));",
                            cls="volVectorField"),
                      np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])),
    "empty list": (field("nonuniform List<scalar> 0();"), np.empty(0)),
    "boundary after": (field(
        "nonuniform List<scalar>\n2\n(\n1\n2\n)\n;",
        "    outer\n    {\n        type            zeroGradient;\n    }\n"),
        np.array([1.0, 2.0])),
}


@pytest.mark.parametrize("name", sorted(INTERNAL_CASES))
def test_internal_field_matches_every_legacy_parser(tmp_path, name):
    text, expected = INTERNAL_CASES[name]
    path = tmp_path / "f"
    path.write_text(text, encoding="utf-8", newline="\n")
    got = volfield.read_internal_field(path)
    np.testing.assert_array_equal(got, expected)
    assert got.dtype == np.float64
    for legacy in (legacy_nozzleplume_internal, legacy_cai2012_internal):
        np.testing.assert_array_equal(got, legacy(text))


def test_crlf_files_read_the_same(tmp_path):
    text, expected = INTERNAL_CASES["nonuniform vector"]
    path = tmp_path / "f"
    path.write_bytes(text.replace("\n", "\r\n").encode())
    np.testing.assert_array_equal(volfield.read_internal_field(path), expected)


def test_a_header_comment_naming_the_keywords_is_not_the_entry(tmp_path):
    text = field("uniform 4;").replace(
        "FoamFile", "// note: internalField and boundaryField are written below\nFoamFile")
    path = tmp_path / "f"
    path.write_text(text, encoding="utf-8")
    np.testing.assert_array_equal(volfield.read_internal_field(path), [4.0])


def test_metadata_says_what_was_written(tmp_path):
    path = tmp_path / "f"
    path.write_text(INTERNAL_CASES["uniform vector"][0], encoding="utf-8")
    uniform = volfield.read_internal(path)
    assert uniform.uniform and uniform.value_type == "vector"
    assert uniform.declared_count is None and uniform.n_components == 3
    path.write_text(INTERNAL_CASES["nonuniform scalar"][0], encoding="utf-8")
    listed = volfield.read_internal(path)
    assert not listed.uniform and listed.value_type == "scalar"
    assert listed.declared_count == 5 and listed.n_values == 5


def test_tensor_types_keep_their_components(tmp_path):
    path = tmp_path / "f"
    path.write_text(field("nonuniform List<symmTensor> 2((1 2 3 4 5 6) (7 8 9 10 11 12));",
                          cls="volSymmTensorField"), encoding="utf-8")
    values = volfield.read_internal(path)
    assert values.value_type == "symmTensor" and values.values.shape == (2, 6)


# --------------------------------------------------------------------------- #
# patches
# --------------------------------------------------------------------------- #

def patch_entry(name, body):
    return f"    {name}\n    {{\n{body}    }}\n"


PATCH_FILE = field("uniform 0;", "".join([
    patch_entry("inflowOuter", "        type            zeroGradient;\n"),
    patch_entry("inflow",
                "        type            fixedValue;\n"
                "        value           nonuniform List<scalar>\n        3\n        (\n"
                "            1e+20\n            2e+20\n            3e+20\n        );\n"),
    patch_entry("plate", "        type            calculated;\n"
                         "        value           uniform 7;\n"),
    patch_entry("wall", "        type            calculated;\n"
                        "        value           nonuniform List<vector> 2((1 2 3) (4 5 6));\n"),
    patch_entry("outer", "        type            zeroGradient;\n"),
]))


@pytest.mark.parametrize("patch,expected", [
    ("inflow", np.array([1e20, 2e20, 3e20])),
    ("plate", np.array([7.0])),
    ("wall", np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])),
])
def test_patch_values_match_both_legacy_parsers(tmp_path, patch, expected):
    path = tmp_path / "f"
    path.write_text(PATCH_FILE, encoding="utf-8")
    got = volfield.read_patch_field(path, patch)
    np.testing.assert_array_equal(got, expected)
    for legacy in (legacy_nozzleplume_patch, legacy_markelov_patch):
        np.testing.assert_array_equal(got, legacy(PATCH_FILE, patch))


def test_a_nested_dictionary_does_not_end_the_patch(tmp_path):
    """markelov1999's parser balanced braces; nozzleplume's stopped at the first."""
    text = field("uniform 0;", patch_entry(
        "inflow", "        type            fixedValue;\n"
                  "        options { note x; }\n"
                  "        value           uniform 5;\n"))
    path = tmp_path / "f"
    path.write_text(text, encoding="utf-8")
    np.testing.assert_array_equal(volfield.read_patch_field(path, "inflow"), [5.0])
    np.testing.assert_array_equal(legacy_markelov_patch(text, "inflow"), [5.0])


def test_a_value_inside_a_nested_dictionary_is_not_the_patch_value(tmp_path):
    text = field("uniform 0;", patch_entry(
        "inflow", "        type            fixedValue;\n"
                  "        sub\n        {\n        value           uniform 9;\n        }\n"
                  "        value           uniform 5;\n"))
    path = tmp_path / "f"
    path.write_text(text, encoding="utf-8")
    np.testing.assert_array_equal(volfield.read_patch_field(path, "inflow"), [5.0])


@pytest.mark.parametrize("patch,message", [
    ("outer", "no 'value' entry"),
    ("panel", "patch 'panel' is not in the boundaryField"),
    ("inflo", "is not in the boundaryField"),
])
def test_patch_errors_say_why(tmp_path, patch, message):
    path = tmp_path / "f"
    path.write_text(PATCH_FILE, encoding="utf-8")
    with pytest.raises(FoamDataError, match=message):
        volfield.read_patch_field(path, patch)


@pytest.mark.parametrize("text,message", [
    (HEADER % ("volScalarField", "f"), "no internalField"),
    (field("nonuniform List<label> 2(1 2);"), "not a field type"),
    (field("uniform (1 2);"), "malformed uniform value"),
    (field("nonuniform List<vector> 2(1 2 3 4);"), "whole number"),
])
def test_unreadable_internal_fields_raise(tmp_path, text, message):
    path = tmp_path / "f"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(FoamDataError, match=message):
        volfield.read_internal(path)


def test_a_missing_file_is_named(tmp_path):
    with pytest.raises(FoamDataError, match="no field file at"):
        volfield.read_internal(tmp_path / "absent")
    assert issubclass(FoamDataError, ValueError)


# --------------------------------------------------------------------------- #
# what is refused: the malformed-field defect
# --------------------------------------------------------------------------- #
#
# Before 2026-10-27 a list declaring four values but holding one came back as
# one value, and every caller expanded any one-row result across the mesh.

def write(tmp_path, text, name="f"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_a_list_shorter_than_it_declares_is_refused(tmp_path):
    path = write(tmp_path, field("nonuniform List<scalar>\n4\n(\n7.5\n)\n;"))
    with pytest.raises(FoamDataError, match="declares 4 values but holds 1"):
        volfield.read_internal(path)


def test_a_vector_list_is_counted_in_vectors(tmp_path):
    path = write(tmp_path, field("nonuniform List<vector> 3((1 2 3) (4 5 6));"))
    with pytest.raises(FoamDataError, match="declares 3 values but holds 2"):
        volfield.read_internal(path)


def test_only_a_uniform_entry_is_expanded(tmp_path):
    uniform = volfield.read_internal(write(tmp_path, field("uniform 2;")))
    np.testing.assert_array_equal(uniform.over(4, "cells"), [2.0] * 4)
    vector = volfield.read_internal(write(tmp_path, field("uniform (1 2 3);"), "g"))
    assert vector.over(3, "cells").shape == (3, 3)
    one = volfield.read_internal(write(tmp_path, field("nonuniform List<scalar> 1(2);"), "h"))
    with pytest.raises(FoamDataError, match="holds 1 values, but there are 4 cells. "
                                            "These are not the same mesh"):
        one.over(4, "cells")
    np.testing.assert_array_equal(one.over(1, "cells"), [2.0])


def test_the_rarefaction_reader_no_longer_broadcasts_a_short_list(tmp_path):
    """The trigger as first observed: a malformed field copied across the mesh."""
    from plumetools.rarefaction import foam

    path = write(tmp_path, field("nonuniform List<scalar>\n4\n(\n7.5\n)\n;"))
    with pytest.raises(FoamDataError, match="declares 4"):
        foam._cell_values(path, 4)
    listed = write(tmp_path, field("nonuniform List<scalar> 1(7.5);"), "g")
    with pytest.raises(FoamDataError, match="not the same mesh"):
        foam._cell_values(listed, 4)
    np.testing.assert_array_equal(
        foam._cell_values(write(tmp_path, field("uniform 7.5;"), "h"), 4), [7.5] * 4)


@pytest.mark.parametrize("text,message", [
    (field("nonuniform List<scalar> 3(1 x 3);"), "not a number"),
    (field("nonuniform List<scalar>;"), "no list"),
    (field("uniform 0;").replace("format      ascii;", "format      binary;"),
     "binary format"),
])
def test_unreadable_payloads_are_refused(tmp_path, text, message):
    with pytest.raises(FoamDataError, match=message):
        volfield.read_internal(write(tmp_path, text))


def test_a_comment_mentioning_binary_is_not_the_format(tmp_path):
    text = field("uniform 3;").replace("FoamFile", "// format binary; would be refused\nFoamFile")
    np.testing.assert_array_equal(volfield.read_internal_field(write(tmp_path, text)), [3.0])


# --------------------------------------------------------------------------- #
# time directories
# --------------------------------------------------------------------------- #

def test_time_directories_are_numeric_and_after_zero(tmp_path):
    for name in ("0", "0.0001", "2e-05", "0.05", "constant", "system", "0.orig",
                 "processor0"):
        (tmp_path / name).mkdir()
    (tmp_path / "0.07").write_text("a file, not a time")
    assert times.time_names(tmp_path) == ["2e-05", "0.0001", "0.05"]
    assert times.time_names(tmp_path, after=0.0001) == ["0.05"]
    assert times.latest_time(tmp_path) == "0.05"
    assert times.time_value("1e-3") == 1e-3 and times.time_value("constant") is None


def test_no_solver_output_is_an_error_that_says_so(tmp_path):
    (tmp_path / "0").mkdir()
    with pytest.raises(FoamDataError, match="solver wrote nothing"):
        times.latest_time(tmp_path)


def test_a_missing_directory_raises_unless_allowed(tmp_path):
    with pytest.raises(FileNotFoundError):
        times.time_directories(tmp_path / "absent")
    assert times.time_directories(tmp_path / "absent", missing_ok=True) == []


# --------------------------------------------------------------------------- #
# patch entries
# --------------------------------------------------------------------------- #

def test_constraint_patches_keep_their_type():
    assert patchfields.patch_field_type("symmetry") == "symmetry"
    assert patchfields.patch_field_type("empty") == "empty"
    assert patchfields.patch_field_type("patch") == "zeroGradient"
    assert patchfields.patch_field_type("wall") == "zeroGradient"


def test_the_nonuniform_entry_is_written_as_before():
    """Verbatim what markelov1999.foamfields._nonuniform_entry wrote at the tag."""
    rows = ["            1e+20", "            2.5"]
    assert patchfields.nonuniform_value_entry("nozzle", "scalar", rows) == [
        "    nozzle",
        "    {",
        "        type            fixedValue;",
        "        value           nonuniform List<scalar>",
        "        2",
        "        (",
        "            1e+20",
        "            2.5",
        "        );",
        "    }",
    ]


# --------------------------------------------------------------------------- #
# the old names still work, and are the shared objects
# --------------------------------------------------------------------------- #

def test_old_import_paths_are_the_shared_objects():
    from plumetools.cai2012 import foamfields as cai_fields
    from plumetools.foamio import fields as legacy_fields
    from plumetools.markelov1999 import foamfields as markelov_fields
    from plumetools.mesh import patches, polymesh
    from plumetools.nozzleplume import foamfields as nozzle_fields
    from plumetools.nozzleplume import meshread, post
    from plumetools.rarefaction import foam

    for module in (cai_fields, legacy_fields, markelov_fields, nozzle_fields):
        assert module.patch_field_type is patchfields.patch_field_type
        assert module.CONSTRAINT_PATCH_TYPES is patchfields.CONSTRAINT_PATCH_TYPES
    assert meshread.PatchGeometry is patches.PatchGeometry
    assert meshread.face_geometry is patches.face_geometry
    assert meshread.read_patches is patches.read_patches
    assert meshread.read_cell_count is polymesh.read_cell_count
    assert post.read_internal_field is volfield.read_internal_field
    assert post.read_patch_field is volfield.read_patch_field
    assert post.latest_time is times.latest_time and foam.latest_time is times.latest_time


def test_error_classes_keep_what_each_handler_caught():
    from plumetools.cai2012 import post as cai_post
    from plumetools.nozzleplume import post as nozzle_post
    from plumetools.rarefaction import foam

    assert nozzle_post.PostError is FoamDataError is foam.PostError
    assert issubclass(cai_post.PostError, FoamDataError)


def test_cai2012_readers_still_raise_their_own_error(tmp_path):
    from plumetools.cai2012 import post as cai_post

    with pytest.raises(cai_post.PostError, match="no field file"):
        cai_post.read_internal_field(tmp_path / "absent")
    (tmp_path / "0").mkdir()
    with pytest.raises(cai_post.PostError, match="wrote nothing"):
        cai_post.latest_time(tmp_path)


def test_markelov_reader_keeps_file_not_found(tmp_path):
    from plumetools.markelov1999 import foamfields as markelov_fields

    with pytest.raises(FileNotFoundError, match="no field file"):
        markelov_fields.read_patch_field(tmp_path / "absent", "inflow")
