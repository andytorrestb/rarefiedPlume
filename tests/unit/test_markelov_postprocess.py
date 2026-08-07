"""Post-processing: log parsing, pressure windows, the study table.

The pressure tests build a synthetic cylinder patch with known face values, so
the area weighting and the window selection can be checked against arithmetic
done here rather than against the implementation.
"""

from __future__ import annotations

import csv
import json
import math
import textwrap
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from plumetools.markelov1999 import postprocess as pp
from plumetools.markelov1999.foamfields import read_patch_field
from plumetools.markelov1999.geometry import MarkelovGeometry

CHECK_MESH_LOG = """\
Mesh stats
    points:           45095
    faces:            113105
    cells:            34346
    faces per cell:   6.36
    boundary patches: 6

Overall number of cells of each type:
    hexahedra:     28584
    prisms:        1166
    wedges:        0
    pyramids:      0
    tet wedges:    0
    tetrahedra:    0
    polyhedra:     4596

Checking geometry...
    Max cell openness = 3.223373561e-16 OK.
    Max aspect ratio = 3.419460717 OK.
    Minimum face area = 1.611609618e-05. Maximum face area = 0.0006571493943.
    Min volume = 4.621948e-08. Max volume = 1.647404671e-05.  Total volume = 0.18.
    Mesh non-orthogonality Max: 37.5775334 average: 10.22883372
    Max skewness = 0.6879351379 OK.

Mesh OK.
End
"""


# --------------------------------------------------------------------------- #
# 1. mesh statistics
# --------------------------------------------------------------------------- #

def test_check_mesh_log_is_parsed(tmp_path):
    path = tmp_path / "log.checkMesh"
    path.write_text(CHECK_MESH_LOG, encoding="utf-8")
    parsed = pp.parse_check_mesh_log(path)

    assert parsed["n_cells"] == 34346
    assert parsed["n_points"] == 45095
    assert parsed["n_faces"] == 113105
    assert parsed["max_non_orthogonality"] == pytest.approx(37.5775334)
    assert parsed["average_non_orthogonality"] == pytest.approx(10.22883372)
    assert parsed["max_skewness"] == pytest.approx(0.6879351379)
    assert parsed["max_aspect_ratio"] == pytest.approx(3.419460717)
    assert parsed["check_mesh_ok"] is True


def test_the_trailing_period_after_a_volume_is_not_swallowed(tmp_path):
    """checkMesh writes "Max volume = 1.647404671e-05.  Cell volumes OK." -- a
    loose character class captures the sentence-ending period and float() then
    fails on '1.647404671e-05.'."""
    path = tmp_path / "log.checkMesh"
    path.write_text(CHECK_MESH_LOG, encoding="utf-8")
    parsed = pp.parse_check_mesh_log(path)

    assert parsed["min_cell_volume_m3"] == pytest.approx(4.621948e-08)
    assert parsed["max_cell_volume_m3"] == pytest.approx(1.647404671e-05)


def test_cell_type_mix_omits_the_types_with_none(tmp_path):
    path = tmp_path / "log.checkMesh"
    path.write_text(CHECK_MESH_LOG, encoding="utf-8")
    types = pp.parse_check_mesh_log(path)["cell_type_counts"]

    assert types == {"hexahedra": 28584, "prisms": 1166, "polyhedra": 4596}
    assert sum(types.values()) == 34346


def test_a_failed_check_mesh_is_recorded_as_such(tmp_path):
    path = tmp_path / "log.checkMesh"
    path.write_text(CHECK_MESH_LOG.replace("Mesh OK.", "***Failed 1 mesh checks."),
                    encoding="utf-8")
    assert pp.parse_check_mesh_log(path)["check_mesh_ok"] is False


def test_a_missing_log_is_not_an_error(tmp_path):
    assert pp.parse_check_mesh_log(tmp_path / "nope") == {}


# --------------------------------------------------------------------------- #
# 2. particle statistics
# --------------------------------------------------------------------------- #

SOLVER_LOG = """\
Time = 2e-07
    Particles inserted              = 724
    Number of dsmc particles        = 1226
ExecutionTime = 0.1 s

Time = 4e-07
    Particles inserted              = 785
    Collisions                      = 2
    Number of dsmc particles        = 2010
ExecutionTime = 0.2 s

Time = 6e-07
    Particles inserted              = 758
    Collisions                      = 4
    Number of dsmc particles        = 2768
ExecutionTime = 9.54 s
End
"""


def test_solver_log_is_parsed(tmp_path):
    path = tmp_path / "log.dsmcFoam"
    path.write_text(SOLVER_LOG, encoding="utf-8")
    stats = pp.parse_solver_log(path)

    assert stats.n_time_steps == 3
    assert stats.initial_particles == 1226
    assert stats.final_particles == 2768
    assert stats.particles_inserted_total == 724 + 785 + 758
    assert stats.collisions_total == 6
    assert stats.runtime_s == pytest.approx(9.54)


def test_mean_particles_uses_the_trailing_steady_fraction(tmp_path):
    path = tmp_path / "log.dsmcFoam"
    path.write_text(SOLVER_LOG, encoding="utf-8")
    # The last half of three samples is the last two: (2010 + 2768)/2.
    assert pp.parse_solver_log(path, steady_fraction=0.5).mean_particles == \
        pytest.approx((2010 + 2768) / 2)


def test_unparsed_fields_are_none_not_zero(tmp_path):
    """A study table showing 0 collisions for an unparsed log would read as a
    result; a blank reads as "not measured", which is the truth."""
    stats = pp.parse_solver_log(tmp_path / "nope")
    assert stats.initial_particles is None
    assert stats.collisions_total is None
    assert stats.runtime_s is None
    assert stats.n_time_steps == 0


def test_a_log_without_collisions_reports_none(tmp_path):
    path = tmp_path / "log.dsmcFoam"
    path.write_text("Time = 1\n    Number of dsmc particles        = 5\n",
                    encoding="utf-8")
    stats = pp.parse_solver_log(path)
    assert stats.collisions_total is None
    assert stats.final_particles == 5


# --------------------------------------------------------------------------- #
# angular windows
# --------------------------------------------------------------------------- #

def test_angular_difference_wraps_at_the_windward_generator():
    """The windward window is centred on 180 deg, exactly where the azimuth
    wraps: a face at -179 deg is 2 deg away, not 359."""
    assert float(pp.angular_difference_deg(-179.0, 180.0)) == pytest.approx(1.0)
    assert float(pp.angular_difference_deg(179.0, 180.0)) == pytest.approx(1.0)
    assert float(pp.angular_difference_deg(0.0, 180.0)) == pytest.approx(180.0)
    assert float(pp.angular_difference_deg(90.0, 0.0)) == pytest.approx(90.0)


def test_angular_difference_is_vectorised():
    out = pp.angular_difference_deg(np.array([-179.0, 179.0, 0.0]), 180.0)
    np.testing.assert_allclose(out, [1.0, 1.0, 180.0])


# --------------------------------------------------------------------------- #
# wall pressure, on a synthetic cylinder patch
# --------------------------------------------------------------------------- #

def build_cylinder_case(tmp_path, values_by_azimuth, *, n_azimuth=72, n_axial=9):
    """A one-patch mesh whose faces tile a cylinder, with a chosen fD per face.

    Each face is a small planar quad on the cylinder surface, wound so its normal
    points radially outward -- out of the fluid and into the wall, which is how
    OpenFOAM winds a boundary face.
    """
    geom = MarkelovGeometry()
    r, xc = geom.cylinder_radius_m, geom.cylinder_centre_x_m

    points, faces, centres, fd = [], [], [], []
    dpsi = 2 * math.pi / n_azimuth
    dz = 0.02
    for j in range(n_axial):
        z0 = (j - n_axial // 2) * dz
        for i in range(n_azimuth):
            psi0, psi1 = i * dpsi - math.pi, (i + 1) * dpsi - math.pi
            base = len(points)
            for psi in (psi0, psi1):
                for z in (z0, z0 + dz):
                    points.append((xc + r * math.cos(psi), r * math.sin(psi), z))
            # order: (psi0,z0) (psi0,z1) (psi1,z1) (psi1,z0) -> outward normal
            faces.append([base + 0, base + 1, base + 3, base + 2])
            psi_mid = 0.5 * (psi0 + psi1)
            centres.append((xc + r * math.cos(psi_mid), r * math.sin(psi_mid),
                            z0 + 0.5 * dz))
            outward = np.array([math.cos(psi_mid), math.sin(psi_mid), 0.0])
            fd.append(values_by_azimuth(math.degrees(psi_mid), z0 + 0.5 * dz) * outward)

    mesh = tmp_path / "constant" / "polyMesh"
    mesh.mkdir(parents=True, exist_ok=True)
    header = ('FoamFile{{version 2.0;format ascii;class {cls};object {obj};}}\n')

    mesh.joinpath("points").write_text(
        header.format(cls="vectorField", obj="points")
        + f"{len(points)}\n(\n"
        + "\n".join(f"({x} {y} {z})" for x, y, z in points) + "\n)\n",
        encoding="utf-8", newline="\n")
    mesh.joinpath("faces").write_text(
        header.format(cls="faceList", obj="faces")
        + f"{len(faces)}\n(\n"
        + "\n".join("4(" + " ".join(str(i) for i in f) + ")" for f in faces) + "\n)\n",
        encoding="utf-8", newline="\n")
    mesh.joinpath("boundary").write_text(
        header.format(cls="polyBoundaryMesh", obj="boundary")
        + "1\n(\n    cylinder\n    {\n        type wall;\n"
        f"        nFaces {len(faces)};\n        startFace 0;\n" + "    }\n)\n",
        encoding="utf-8", newline="\n")

    time_dir = tmp_path / "0.004"
    time_dir.mkdir(exist_ok=True)
    rows = "\n".join(f"({v[0]} {v[1]} {v[2]})" for v in fd)
    time_dir.joinpath("fDMean").write_text(
        header.format(cls="volVectorField", obj="fDMean")
        + "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform (0 0 0);\n"
        "boundaryField\n{\n    cylinder\n    {\n        type calculated;\n"
        f"        value nonuniform List<vector>\n{len(fd)}\n(\n{rows}\n)\n;\n"
        "    }\n}\n", encoding="utf-8", newline="\n")
    return geom


@pytest.fixture
def pressure_cfg():
    from test_markelov_mesh import make_cfg
    return make_cfg()


def test_windward_and_leeward_pressures_are_measured_where_expected(
        tmp_path, pressure_cfg):
    """Windward is the source-facing generator at azimuth 180; leeward is the
    base at 0. Swapping them silently swaps the whole result."""
    def values(psi_deg, z):
        return 10.0 if abs(pp.angular_difference_deg(psi_deg, 180.0)) < 30 else 2.0

    geom = build_cylinder_case(tmp_path, values)
    windows, extras = pp.wall_pressure(tmp_path, pressure_cfg, geom)
    by_name = {w.name: w for w in windows}

    assert by_name["windward"].pressure_pa == pytest.approx(10.0, rel=1e-6)
    assert by_name["leeward"].pressure_pa == pytest.approx(2.0, rel=1e-6)
    assert extras["pressure_ratio"] == pytest.approx(5.0, rel=1e-6)


def test_pressure_is_the_wall_normal_component_not_the_magnitude(
        tmp_path, pressure_cfg):
    """fD carries shear as well as pressure. The normal component is the
    pressure; the tangential part is reported separately, not folded in."""
    def values(psi_deg, z):
        return 7.0

    geom = build_cylinder_case(tmp_path, values)
    windows, _ = pp.wall_pressure(tmp_path, pressure_cfg, geom)
    for window in windows:
        assert window.pressure_pa == pytest.approx(7.0, rel=1e-6)
        assert window.shear_pa == pytest.approx(0.0, abs=1e-9)


def test_a_window_averages_many_faces_not_one(tmp_path, pressure_cfg):
    """Requirement 13: a single-face reading would depend on where snappyHexMesh
    happened to put that face."""
    geom = build_cylinder_case(tmp_path, lambda psi, z: 1.0)
    windows, _ = pp.wall_pressure(tmp_path, pressure_cfg, geom)
    for window in windows:
        assert window.n_faces > 1
        assert window.area_m2 > 0.0


def test_the_window_respects_its_axial_limit(tmp_path, pressure_cfg):
    """Faces outside the axial half-height must not contribute, so a value that
    only exists there cannot move the average."""
    def values(psi_deg, z):
        return 100.0 if abs(z) > 0.1 else 1.0

    geom = build_cylinder_case(tmp_path, values)
    windows, _ = pp.wall_pressure(tmp_path, pressure_cfg, geom)
    for window in windows:
        assert window.pressure_pa == pytest.approx(1.0, rel=1e-6)


def test_the_window_respects_its_angular_limit(tmp_path, pressure_cfg):
    def values(psi_deg, z):
        return 100.0 if abs(pp.angular_difference_deg(psi_deg, 90.0)) < 5 else 1.0

    geom = build_cylinder_case(tmp_path, values)
    windows, _ = pp.wall_pressure(tmp_path, pressure_cfg, geom)
    by_name = {w.name: w for w in windows}
    assert by_name["windward"].pressure_pa == pytest.approx(1.0, rel=1e-6)
    assert by_name["midspan_side"].pressure_pa > 1.0   # this window contains it


def test_a_window_capturing_no_faces_is_an_error(tmp_path, pressure_cfg):
    """Rather than a NaN that would reach the study table."""
    narrow = replace(pressure_cfg, post=replace(pressure_cfg.post, windows=(
        {"name": "windward", "azimuth_deg": 180.0, "half_angle_deg": 0.001,
         "axial_half_height_m": 1e-9},)))
    geom = build_cylinder_case(tmp_path, lambda psi, z: 1.0)
    with pytest.raises(ValueError, match="captured no faces"):
        pp.wall_pressure(tmp_path, narrow, geom)


def test_a_nonpositive_leeward_pressure_suppresses_the_ratio(tmp_path, pressure_cfg):
    """An under-sampled wall averages negative before enough hits accumulate.
    Dividing by it would put a meaningless -- possibly large -- number in the
    table."""
    def values(psi_deg, z):
        return 10.0 if abs(pp.angular_difference_deg(psi_deg, 180.0)) < 30 else -0.001

    geom = build_cylinder_case(tmp_path, values)
    _windows, extras = pp.wall_pressure(tmp_path, pressure_cfg, geom)
    assert extras["pressure_ratio"] is None
    assert "too few particles" in extras["pressure_ratio_note"]


def test_missing_results_are_reported_not_guessed(tmp_path, pressure_cfg):
    geom = MarkelovGeometry()
    with pytest.raises(FileNotFoundError, match="no time directory"):
        pp.wall_pressure(tmp_path, pressure_cfg, geom)


def test_latest_time_directory_excludes_zero(tmp_path):
    for name in ("0", "0.002", "0.004", "constant", "system"):
        (tmp_path / name).mkdir()
    assert pp.latest_time_dir(tmp_path).name == "0.004"


# --------------------------------------------------------------------------- #
# field reading
# --------------------------------------------------------------------------- #

def test_a_zero_gradient_patch_field_is_an_explicit_error(tmp_path):
    """zeroGradient writes no values, so the solver's wall data was discarded at
    write time. That is a different fact from "the pressure was zero"."""
    path = tmp_path / "fD"
    path.write_text(
        "FoamFile{version 2.0;format ascii;class volVectorField;object fD;}\n"
        "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform (0 0 0);\n"
        "boundaryField\n{\n    cylinder\n    {\n        type zeroGradient;\n"
        "    }\n}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no 'value' entry"):
        read_patch_field(path, "cylinder")


def test_a_uniform_patch_value_is_read(tmp_path):
    path = tmp_path / "q"
    path.write_text(
        "FoamFile{version 2.0;format ascii;class volScalarField;object q;}\n"
        "dimensions [1 0 -3 0 0 0 0];\ninternalField uniform 0;\n"
        "boundaryField\n{\n    plate\n    {\n        type calculated;\n"
        "        value uniform 3.5;\n    }\n}\n", encoding="utf-8")
    np.testing.assert_allclose(read_patch_field(path, "plate"), [3.5])


def test_a_missing_patch_is_an_error(tmp_path):
    path = tmp_path / "q"
    path.write_text(
        "FoamFile{version 2.0;class volScalarField;object q;}\n"
        "boundaryField\n{\n    plate\n    {\n        type calculated;\n"
        "        value uniform 1;\n    }\n}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not in the boundaryField"):
        read_patch_field(path, "cylinder")


# --------------------------------------------------------------------------- #
# the study table
# --------------------------------------------------------------------------- #

def test_study_row_pulls_the_expected_columns():
    summary = {
        "case": {"name": "p025psi", "pressure_psi": 25.0,
                 "pressure_pa": 172368.93, "gap_in": 6.0},
        "pressure": {
            "windows": [
                {"name": "windward", "pressure_pa": 1.5e-3},
                {"name": "leeward", "pressure_pa": 3.0e-4},
            ],
            "pressure_ratio": 5.0,
        },
        "mesh": {"n_cells": 34346},
        "particles": {"mean_particles": 3.4e6, "runtime_s": 1234.5},
        "resolution_audit": {"region_occupancy": {"cylinder": 19.8, "plate": 3.4}},
    }
    row = pp.summary_row(summary)

    assert row["case_name"] == "p025psi"
    assert row["pressure_psi"] == 25.0
    assert row["windward_pressure_pa"] == 1.5e-3
    assert row["leeward_pressure_pa"] == 3.0e-4
    assert row["pressure_ratio"] == 5.0
    assert row["n_cells"] == 34346
    assert row["particles_per_cell_cylinder"] == 19.8
    assert row["particles_per_cell_plate"] == 3.4
    assert row["runtime_s"] == 1234.5


def test_missing_values_become_blank_not_zero():
    """A zero would read as a result; a blank reads as "not measured"."""
    row = pp.summary_row({"case": {"name": "p005psi"}})
    assert row["windward_pressure_pa"] == ""
    assert row["pressure_ratio"] == ""
    assert row["n_cells"] == ""


def test_the_table_has_the_required_columns():
    for column in ("pressure_psi", "pressure_pa", "gap_in",
                   "windward_pressure_pa", "leeward_pressure_pa", "pressure_ratio",
                   "n_cells", "n_particles_mean", "particles_per_cell_cylinder",
                   "particles_per_cell_plate", "runtime_s"):
        assert column in pp.STUDY_COLUMNS


def test_the_table_is_written_in_the_given_order(tmp_path):
    rows = [pp.summary_row({"case": {"name": n, "pressure_psi": p}})
            for n, p in (("p005psi", 5), ("p025psi", 25), ("p100psi", 100))]
    path = pp.write_study_table(tmp_path / "results" / "study-table.csv", rows)

    with open(path, encoding="utf-8") as f:
        written = list(csv.DictReader(f))
    assert [r["case_name"] for r in written] == ["p005psi", "p025psi", "p100psi"]
    assert list(written[0]) == list(pp.STUDY_COLUMNS)


def test_the_table_is_byte_reproducible(tmp_path):
    rows = [pp.summary_row({"case": {"name": "p005psi", "pressure_psi": 5}})]
    a = pp.write_study_table(tmp_path / "a.csv", rows).read_bytes()
    b = pp.write_study_table(tmp_path / "b.csv", rows).read_bytes()
    assert a == b


def test_plots_are_skipped_gracefully_without_matplotlib(tmp_path, monkeypatch):
    """matplotlib is an optional extra; its absence must not fail the whole
    post-processing run."""
    rows = [pp.summary_row({"case": {"name": "p005psi", "pressure_psi": 5}})]
    table = pp.write_study_table(tmp_path / "t.csv", rows)

    import builtins
    real_import = builtins.__import__

    def no_matplotlib(name, *args, **kwargs):
        if name.startswith("matplotlib"):
            raise ImportError("no matplotlib")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_matplotlib)
    assert pp.plot_study(table, tmp_path) == []


def test_summary_is_json_serialisable():
    summary = {"case": {"name": "p005psi"}, "mesh": {"n_cells": 1}}
    json.dumps(pp.summary_row(summary))
