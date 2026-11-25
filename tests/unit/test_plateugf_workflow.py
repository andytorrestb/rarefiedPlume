"""The plate family's case workflow: restart plan, build check, smoke check."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from plumetools.plateugf import workflow

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "cases" / "plate-unigas" / "baseCase"

HEADER = """FoamFile
{{
    version 2.0; format ascii; class {cls}; object {name};
}}
dimensions [0 0 0 0 0 0 0];
internalField uniform 1;
boundaryField
{{
    plate
    {{
        type            calculated;
        value           nonuniform List<{kind}> 3({values});
    }}
    vacuum
    {{
        type            calculated;
        value           uniform 0;
    }}
}}
"""


def _control(path: Path, end: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"application uniGasFoam;\nendTime         {end};\nwriteInterval 0.01;\n",
                    encoding="utf-8")


@pytest.fixture
def case_dir(tmp_path):
    shutil.copy(TEMPLATE / "case.yaml", tmp_path / "case.yaml")
    shutil.copy(TEMPLATE / "unigas.yaml", tmp_path / "unigas.yaml")
    _control(tmp_path / "system" / "controlDict.startup", 0.06)
    _control(tmp_path / "system" / "controlDict.sampling", 0.09)
    return tmp_path


def _plate_fields(time_dir: Path, p_values: str = "1 2 0") -> None:
    time_dir.mkdir(parents=True, exist_ok=True)
    for name in ("p_Ar", "surfaceShearStress_Ar", "surfaceHeatTransfer_Ar", "rhoN_Ar"):
        values = p_values if name == "p_Ar" else "1 1 1"
        (time_dir / name).write_text(HEADER.format(cls="volScalarField", name=name,
                                                   kind="scalar", values=values), encoding="utf-8")
    (time_dir / "fD_Ar").write_text(HEADER.format(
        cls="volVectorField", name="fD_Ar", kind="vector",
        values="(1 0 0) (2 0 0) (0 0 0)"), encoding="utf-8")


def test_run_window_from_control_dicts(case_dir):
    assert workflow.run_window(case_dir) == (0.06, 0.09)


def test_resume_plan_fresh_case_is_startup(case_dir):
    plan = workflow.resume_plan(case_dir, 1)
    assert plan.leg == "startup"
    assert plan.refusal is None
    assert "LEG=startup" in plan.shell()


def test_check_build_mismatch_reports(case_dir, tmp_path_factory, monkeypatch):
    appbin = tmp_path_factory.mktemp("appbin")
    (appbin / "uniGasFoam.commit").write_text("0123456789abcdef\n", encoding="utf-8")
    monkeypatch.setenv("FOAM_USER_APPBIN", str(appbin))
    problem = workflow.check_build(case_dir)
    assert problem is not None and "bd8ba5c" in problem
    (appbin / "uniGasFoam.commit").write_text("bd8ba5c3d9ffbbf95a30629eae7012b913188f4a\n",
                                              encoding="utf-8")
    assert workflow.check_build(case_dir) is None


def test_smoke_check_fails_without_plate_fields(case_dir):
    (case_dir / "0.01").mkdir()
    results = workflow.smoke_check(case_dir, "0.01")
    assert any(not ok for _, ok, _ in results)


def test_smoke_check_passes_with_loaded_plate(case_dir):
    _plate_fields(case_dir / "0.01")
    results = workflow.smoke_check(case_dir, "0.01")
    assert all(ok for _, ok, _ in results), results


def test_smoke_check_fails_on_unloaded_plate(case_dir):
    _plate_fields(case_dir / "0.01", p_values="0 0 0")
    results = dict((n, ok) for n, ok, _ in workflow.smoke_check(case_dir, "0.01"))
    assert results["plate_loaded"] is False


BOUNDARY = """FoamFile
{
    version 2.0; format ascii; class polyBoundaryMesh; object boundary;
}
2
(
    plate
    {
        type            wall;
        nFaces          1;
        startFace       4;
    }
    symmetry
    {
        type            symmetry;
        nFaces          2;
        startFace       5;
    }
)
"""

OWNER = """FoamFile
{
    version 2.0; format ascii; class labelList; object owner;
}
7
(
0 0 1 1 2 0 1
)
"""


def _symmetry_mesh(case_dir: Path, rho_values: str) -> None:
    mesh = case_dir / "constant" / "polyMesh"
    mesh.mkdir(parents=True, exist_ok=True)
    (mesh / "boundary").write_text(BOUNDARY, encoding="utf-8")
    (mesh / "owner").write_text(OWNER, encoding="utf-8")
    _plate_fields(case_dir / "0.01")
    text = (case_dir / "0.01" / "rhoN_Ar").read_text(encoding="utf-8")
    text = text.replace("internalField uniform 1;",
                        f"internalField nonuniform List<scalar> 3({rho_values});")
    (case_dir / "0.01" / "rhoN_Ar").write_text(text, encoding="utf-8")


def test_symmetry_check_reads_cells_beside_the_plane(case_dir):
    _symmetry_mesh(case_dir, "1e14 2e14 0")          # symmetry faces 5, 6 -> cells 0, 1
    results = {n: ok for n, ok, _ in workflow.smoke_check(case_dir, "0.01")}
    assert results["gas_at_symmetry"] is True


def test_symmetry_check_fails_when_plane_cells_are_empty(case_dir):
    _symmetry_mesh(case_dir, "0 0 5e14")
    results = {n: ok for n, ok, _ in workflow.smoke_check(case_dir, "0.01")}
    assert results["gas_at_symmetry"] is False


def test_leg_after_average_start_is_sampling(case_dir):
    t = case_dir / "0.07"
    (t / "lagrangian" / "uniGas").mkdir(parents=True)
    (t / "uniform").mkdir()
    (t / "uniform" / "time").write_text("deltaT 1e-06;\n", encoding="utf-8")
    plan = workflow.resume_plan(case_dir, 1)
    assert plan.base.resume == "0.07"
    assert plan.leg == "sampling"


POINTS = """FoamFile
{
    version 2.0; format ascii; class vectorField; object points;
}
6
(
(0 0 0)
(1 3e-05 0)
(1 -2e-05 1)
(0 0 1)
(2 5e-05 0)
(2 4e-05 1)
)
"""

FACES = """FoamFile
{
    version 2.0; format ascii; class faceList; object faces;
}
2
(
4(1 4 5 2)
4(0 1 2 3)
)
"""

SYM_BOUNDARY = """FoamFile
{
    version 2.0; format ascii; class polyBoundaryMesh; object boundary;
}
1
(
    symmetry
    {
        type            symmetry;
        nFaces          1;
        startFace       1;
    }
)
"""


def test_flatten_symmetry_puts_symmetry_points_on_the_plane(tmp_path):
    from plumetools.mesh.polymesh import read_points_array
    mesh = tmp_path / "constant" / "polyMesh"
    mesh.mkdir(parents=True)
    (mesh / "points").write_text(POINTS, encoding="utf-8")
    (mesh / "faces").write_text(FACES, encoding="utf-8")
    (mesh / "boundary").write_text(SYM_BOUNDARY, encoding="utf-8")
    moved = workflow.flatten_symmetry(tmp_path)
    points = read_points_array(tmp_path)
    assert moved == 2
    assert points[[0, 1, 2, 3], 1].tolist() == [0.0, 0.0, 0.0, 0.0]
    assert points[4, 1] == 5e-05 and points[5, 1] == 4e-05   # not on the symmetry face
    assert points[:, 0].tolist() == [0, 1, 1, 0, 2, 2]
