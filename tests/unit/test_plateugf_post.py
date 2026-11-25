"""Plate post-processing and the Cai comparison on a synthetic plate (phase B, B3)."""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest
import yaml

from plumetools.plateugf import post

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "cases" / "plate-unigas" / "baseCase"
X = 5.0
# Upstream face: a 3 (y) x 2 (z) grid of 1 m quads at x = 5, then two edge faces.
YS, ZS = [0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0]


def _mesh(case: Path) -> int:
    pts, faces = [], []
    index = {}
    for j, y in enumerate(YS):
        for k, z in enumerate(ZS):
            index[(j, k)] = len(pts)
            pts.append((X, y, z))
    for k in range(2):
        for j in range(3):
            a, b, c, d = index[(j, k)], index[(j + 1, k)], index[(j + 1, k + 1)], index[(j, k + 1)]
            faces.append((a, b, c, d))                      # normal +x
    base = len(pts)
    pts += [(X, 3.0, 0.0), (X + 0.25, 3.0, 0.0), (X + 0.25, 3.0, 1.0), (X, 3.0, 1.0),
            (X, 3.0, 1.0), (X + 0.25, 3.0, 1.0), (X + 0.25, 3.0, 2.0), (X, 3.0, 2.0)]
    faces += [(base, base + 1, base + 2, base + 3), (base + 4, base + 5, base + 6, base + 7)]
    mesh = case / "constant" / "polyMesh"
    mesh.mkdir(parents=True)
    head = "FoamFile\n{{\n    version 2.0; format ascii; class {cls}; object {obj};\n}}\n"
    (mesh / "points").write_text(head.format(cls="vectorField", obj="points")
                                 + f"{len(pts)}\n(\n" + "\n".join(f"({x} {y} {z})" for x, y, z in pts)
                                 + "\n)\n", encoding="utf-8")
    (mesh / "faces").write_text(head.format(cls="faceList", obj="faces")
                                + f"{len(faces)}\n(\n" + "\n".join("4(%d %d %d %d)" % f for f in faces)
                                + "\n)\n", encoding="utf-8")
    (mesh / "boundary").write_text(head.format(cls="polyBoundaryMesh", obj="boundary")
                                   + "1\n(\n    plate\n    {\n        type wall;\n"
                                   f"        nFaces {len(faces)};\n        startFace 0;\n    }}\n)\n",
                                   encoding="utf-8")
    return len(faces)


def _field(path: Path, name: str, kind: str, values: list[str]) -> None:
    cls = "volVectorField" if kind == "vector" else "volScalarField"
    path.write_text(
        "FoamFile\n{\n    version 2.0; format ascii; class %s; object %s;\n}\n" % (cls, name)
        + "dimensions [1 -1 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n"
        + "    plate\n    {\n        type            calculated;\n"
        + f"        value           nonuniform List<{kind}> {len(values)}({' '.join(values)});\n"
        + "    }\n}\n", encoding="utf-8")


@pytest.fixture
def case(tmp_path):
    shutil.copy(TEMPLATE / "unigas.yaml", tmp_path / "unigas.yaml")
    cfg = yaml.safe_load((TEMPLATE / "case.yaml").read_text(encoding="utf-8"))
    cfg["plate"]["standoff_m"] = X
    cfg["mesh"]["refinement_regions"] = [
        {"name": "exit", "level": 5, "shape": "cone", "x_max_m": 0.5,
         "radius1_m": 0.4, "radius2_m": 0.4}]
    (tmp_path / "case.yaml").write_text(yaml.safe_dump(cfg), encoding="utf-8")
    n = _mesh(tmp_path)
    t = tmp_path / "0.05"
    t.mkdir()
    fd = [f"({1e-6 * (i + 1)} {2e-8} {3e-8})" for i in range(n - 2)] + ["(0 9e-7 0)", "(0 9e-7 0)"]
    _field(t / "fD_Ar", "fD_Ar", "vector", fd)
    _field(t / "surfaceHeatTransfer_Ar", "surfaceHeatTransfer_Ar", "scalar",
           [str(1e-4 * (i + 1)) for i in range(n)])
    return tmp_path


def test_read_plate_uses_unigas_fields(case):
    loads = post.read_plate(case)
    assert loads.time == "0.05"
    assert loads.force_density[0] == pytest.approx([1e-6, 2e-8, 3e-8])
    assert loads.heat_flux[2] == pytest.approx(3e-4)
    assert loads.multiplicity == 4
    assert loads.pressure[:6] == pytest.approx([1e-6 * (i + 1) for i in range(6)])


def test_edge_faces_are_left_out_of_the_comparison(case):
    result = post.compare(case)
    assert result["faces"] == 6
    assert len(result["table"]["y_m"]) == 6
    assert all(c > 0 for c in result["table"]["p_cai_pa"])


def test_missing_mean_fields_are_refused(case):
    (case / "0.05" / "fD_Ar").unlink()
    with pytest.raises(post.PostError, match="fD_Ar"):
        post.read_plate(case)


def test_profiles_take_the_symmetry_rows(case):
    loads = post.read_plate(case)
    prof = post.profiles(loads, {"p": loads.pressure})
    assert prof["along_y"]["z_m"] == pytest.approx([0.5, 0.5, 0.5])
    assert prof["along_y"]["y_m"] == pytest.approx([0.5, 1.5, 2.5])
    assert prof["along_z"]["y_m"] == pytest.approx([0.5, 0.5])
    assert prof["along_z"]["p"] == pytest.approx([1e-6, 4e-6])


def test_run_post_writes_results_and_figures(case):
    written = {p.name for p in post.run_post(case)}
    for name in ("plate_loads.csv", "plate_summary.yaml", "plate_profiles.csv",
                 "cai_comparison.csv", "cai_comparison.yaml", "plate-pressure.png",
                 "plate-shear.png", "cai-profiles.png", "cai-ratio-maps.png"):
        assert name in written, name
    summary = yaml.safe_load((case / "results" / "cai_comparison.yaml").read_text(encoding="utf-8"))
    assert summary["whole_plate"]["normal_force_sim_n"] == pytest.approx(4 * 21e-6, rel=1e-9)


def test_profiles_survive_snapping_noise(case):
    # snappyHexMesh leaves face centres micrometres apart: nudge one point of
    # the z = 0 row so two of its faces sit a little lower than the third.
    points = case / "constant" / "polyMesh" / "points"
    text = points.read_text(encoding="utf-8")
    points.write_text(text.replace(f"({X} 2.0 0.0)", f"({X} 2.0 -1e-05)", 1), encoding="utf-8")
    loads = post.read_plate(case)
    prof = post.profiles(loads, {"p": loads.pressure})
    assert len(prof["along_y"]["y_m"]) == 3
    assert len(prof["along_z"]["z_m"]) == 2


def _controls(case: Path, startup_end: float, sampling_end: float) -> None:
    system = case / "system"
    system.mkdir(exist_ok=True)
    for leg, end in (("startup", startup_end), ("sampling", sampling_end)):
        (system / f"controlDict.{leg}").write_text(f"endTime         {end};\n", encoding="utf-8")


def test_a_write_before_averaging_is_refused(case):
    _controls(case, startup_end=0.1, sampling_end=0.2)
    with pytest.raises(post.PostError, match="before averaging starts"):
        post.run_post(case)
    assert (case / "results" / "cai_comparison.yaml").exists() is False


def test_a_startup_write_can_be_drawn_when_asked_and_is_labelled(case):
    _controls(case, startup_end=0.1, sampling_end=0.2)
    post.run_post(case, allow_startup=True)
    summary = yaml.safe_load((case / "results" / "cai_comparison.yaml").read_text(encoding="utf-8"))
    assert summary["leg"] == "startup"


def test_a_sampling_write_records_its_window(case):
    _controls(case, startup_end=0.04, sampling_end=0.06)
    post.run_post(case)
    summary = yaml.safe_load((case / "results" / "cai_comparison.yaml").read_text(encoding="utf-8"))
    assert summary["leg"] == "sampling"
    assert summary["averaging_window_s"] == pytest.approx(0.01)


def test_centre_metric_is_a_disk_average_not_one_face(case):
    result = post.compare(case)
    centre = result["summary"]["centre_disk"]
    t = result["table"]
    r = np.hypot(t["y_m"], t["z_m"])
    inside = r <= centre["radius_m"]
    assert inside.sum() >= 2                     # never one face alone
    expected = (t["p_sim_pa"][inside] * t["area_m2"][inside]).sum() / t["area_m2"][inside].sum()
    assert centre["p_sim_pa"] == pytest.approx(expected)
    assert centre["faces"] == int(inside.sum())
    shear = result["summary"]["shear_force"]
    assert shear["sim_n"] == pytest.approx(4 * (t["tau_sim_pa"] * t["area_m2"]).sum())
