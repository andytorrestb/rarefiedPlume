"""The interface track end to end, without OpenFOAM.

Generation of ``Interface/`` from ``interfaces.yaml`` (and nothing in
``Cases/``), the policy that keeps validation runs at Kn = 0.01 out of the
production matrix, the mesh plan, the per-face ``0/`` fields on a synthetic
hemisphere patch in patch order, the identity that refuses old parcels on a new
mesh, and the PBS guard and label.
"""

from __future__ import annotations

import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from _helpers import make_case, patch
from plumetools.cai2012.study import load_study
from plumetools.cai2012ext import (
    generation,
    hpc,
    interface,
    interface_case,
    interface_mesh,
    interface_study,
    policy,
)
from plumetools.markelov1999.foamfields import read_patch_field

REPO = Path(__file__).resolve().parents[2]
FAMILY = REPO / "cases" / "cai2012-kn-extension"
SMOKE = "smoke_Kn0p001_B0p05"


def make_family(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    ignore = shutil.ignore_patterns("system", "constant", "results",
                                    "case-summary.json")
    shutil.copytree(FAMILY / "baseCase", root / "baseCase", ignore=ignore)
    shutil.copytree(FAMILY / "baseCaseInterface", root / "baseCaseInterface")
    for name in ("study.yaml", "interfaces.yaml", "viz.yaml"):
        shutil.copy2(FAMILY / name, root / name)
    return root


def run_generator(family: Path, *args) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(FAMILY / "generate_interfaces.py"),
         "--study", str(family / "study.yaml"), *args],
        capture_output=True, text=True, timeout=900)


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    family = make_family(tmp_path_factory.mktemp("iface") / "family")
    result = run_generator(family)
    assert result.returncode == 0, result.stdout + result.stderr
    return family


def manifest_of(family: Path) -> dict:
    return yaml.safe_load((family / "interface-manifest.yaml").read_text(encoding="utf-8"))


def clone(generated: Path, tmp_path: Path) -> Path:
    target = tmp_path / "family"
    shutil.copytree(generated, target)
    return target


def edit_interfaces(family: Path, change) -> None:
    path = family / "interfaces.yaml"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    change(data)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


# --------------------------------------------------------------------------- #
# generation
# --------------------------------------------------------------------------- #

def test_the_track_is_generated_beside_the_matrix_not_in_it(generated):
    names = sorted(p.name for p in (generated / "Interface").iterdir())
    assert names == sorted(["Kn0p001_np1x_sphere", "Kn0p0001_np1x_sphere",
                            "val_Kn0p01_B0p02", "val_Kn0p01_B0p03", "val_Kn0p01_B0p05",
                            SMOKE])
    assert not (generated / "Cases").exists()
    assert not (generated / "manifest.yaml").exists()
    for name in names:
        case = generated / "Interface" / name
        # The overlay replaced the direct-exit scripts; the rest are baseCase's.
        assert (case / "Allmesh").read_bytes() == (
            generated / "baseCaseInterface" / "Allmesh").read_bytes()
        assert (case / "runInflow.py").read_bytes() == (
            generated / "baseCaseInterface" / "runInflow.py").read_bytes()
        assert (case / "runCase.sh").read_bytes() == (
            generated / "baseCase" / "runCase.sh").read_bytes()
        assert interface.load_interface(case) is not None


def test_the_manifest_records_every_verdict(generated):
    manifest = manifest_of(generated)
    rows = {row["name"]: row for row in manifest["cases"]}
    for name in ("val_Kn0p01_B0p02", "val_Kn0p01_B0p03", "val_Kn0p01_B0p05"):
        assert rows[name]["Kn"] == pytest.approx(0.01)
        assert rows[name]["purpose"] == "interface-validation"
        assert rows[name]["selection_status"] == interface.STATUS_NO_INTERFACE
        assert rows[name]["status"] == interface_case.STATUS_INFEASIBLE
    # interfaces.yaml has stated both production cases' core angle and validity
    # radius since 2026-10-25, so a sphere is selected and planned for each.
    # Neither may run as production: both miss numerical criteria (the mesh
    # budget, cell / lambda) and no overlap evidence is named.
    for name in ("Kn0p001_np1x_sphere", "Kn0p0001_np1x_sphere"):
        row = rows[name]
        assert row["selection_status"] == interface.STATUS_ADMISSIBLE
        assert row["status"] == interface_case.STATUS_UNSUITABLE
        assert row["numerical_status"] == interface_case.STATUS_UNSUITABLE
        assert row["coupling_validation"] == "pending"
        assert not row["runnable"]
    assert rows[SMOKE]["status"] == interface_case.STATUS_SMOKE
    assert rows[SMOKE]["n_cells"] <= 50_000
    assert rows[SMOKE]["estimated_parcels"] <= 200_000
    validation = manifest["validation"]
    assert validation["overlap_study_possible"] is False
    assert validation["B_at_smallest_admissible_radius"]["axis"] > 0.05


def test_generation_is_idempotent(generated, tmp_path):
    family = clone(generated, tmp_path)
    result = run_generator(family)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count("keep (identical)") == 6


def test_a_physics_change_is_a_conflict_and_evidence_is_not(generated, tmp_path):
    family = clone(generated, tmp_path)

    def evidence(data):
        data["production"][0]["interface"]["validation_evidence"] = "x/overlap.yaml"
    edit_interfaces(family, evidence)
    result = run_generator(family)
    assert result.returncode == 0, result.stderr
    assert "STALE" in result.stdout

    def threshold(data):
        data["production"][0]["interface"]["bird_threshold"] = 0.03
    edit_interfaces(family, threshold)
    result = run_generator(family)
    assert result.returncode == 1
    assert "conflicts" in result.stderr


def _load(family: Path):
    return interface_study.load(family / "interfaces.yaml",
                                load_study(family / "study.yaml"))


@pytest.mark.parametrize("change,match", [
    (lambda d: d["production"][0].update(kn_case="Kn1"), "only at"),
    (lambda d: d["validation"].update(knudsen=0.1), "full-domain benchmark"),
    (lambda d: d["production"][0].update(overrides={"exit": {"T0_K": 310.0}}),
     "may not override"),
    (lambda d: d["smoke"][0].update(name="Kn10_np1x"), "collides"),
    (lambda d: d["production"][0]["interface"].update(purpose="development-smoke"),
     "purpose must be"),
    (lambda d: d["validation"]["interface"].update(bird_threshold=0.05),
     "bird_thresholds"),
])
def test_the_track_policy_refuses(tmp_path, change, match):
    family = make_family(tmp_path / "family")
    edit_interfaces(family, change)
    with pytest.raises(policy.ExtensionError, match=match):
        _load(family)


def test_the_direct_matrix_still_refuses_kn_001(tmp_path):
    """The exception is the interface track's alone."""
    study = load_study(FAMILY / "study.yaml")
    policy.check_study(study)
    manifest = {"cases": [{"name": "val", "Kn": 0.01, "path": "Cases/val"}]}
    assert policy.check_manifest(manifest)


# --------------------------------------------------------------------------- #
# the mesh plan
# --------------------------------------------------------------------------- #

def test_retained_cells_are_the_lattice_outside_the_sphere(generated):
    derived = interface_case.derive(generated / "Interface" / SMOKE)
    plan = derived.mesh_plan
    xs = np.concatenate([s.cell_centres() for s in plan.lattice.x_segments])
    ys = np.concatenate([s.cell_centres() for s in plan.lattice.y_segments])
    zs = np.concatenate([s.cell_centres() for s in plan.lattice.z_segments])
    x, y, z = np.meshgrid(xs, ys, zs, indexing="ij")
    outside = int(np.count_nonzero(x * x + y * y + z * z > plan.radius_m ** 2))
    assert plan.retained_cells == outside
    assert plan.lattice.n_cells > outside
    assert plan.core_x_m > plan.radius_m and plan.core_half_m > plan.radius_m


def test_mesh_dictionaries_carve_the_sphere(generated, tmp_path):
    case = clone(generated, tmp_path) / "Interface" / SMOKE
    derived = interface_case.derive(case)
    interface_case.write_mesh_dictionaries(derived)
    snappy = (case / "system" / "snappyHexMeshDict").read_text(encoding="utf-8")
    assert "searchableSphere" in snappy and "nozzle" in snappy
    import re
    radius = float(re.search(r"radius\s+([-+0-9.eE]+);", snappy).group(1))
    assert radius == pytest.approx(derived.selection.radius_m, rel=1e-6)
    assert not (case / "system" / "topoSetDict").exists()
    properties = (case / "constant" / "dsmcProperties").read_text(encoding="utf-8")
    assert "INTERFACE MODE" in properties
    assert f"{derived.run.n_equivalent_particles:.6e}" in properties


def test_an_infeasible_case_writes_no_mesh(generated, tmp_path):
    case = clone(generated, tmp_path) / "Interface" / "val_Kn0p01_B0p05"
    derived = interface_case.derive(case)
    with pytest.raises(policy.ExtensionError, match="NO_ADMISSIBLE_INTERFACE"):
        interface_case.write_mesh_dictionaries(derived)
    assert not (case / "system" / "blockMeshDict").exists()


# --------------------------------------------------------------------------- #
# the inflow on a synthetic hemisphere patch
# --------------------------------------------------------------------------- #

def hemisphere_mesh(case: Path, radius: float, *, n_theta=18, n_phi=36,
                    reverse=False, scale=1.0) -> int:
    """A boundary-only mesh: the sphere as 'nozzle', OpenFOAM's winding."""
    rings = [(0.0, 0.0)]
    points = [(radius * scale, 0.0, 0.0)]
    index = {}
    for i in range(1, n_theta + 1):
        theta = 0.5 * math.pi * i / n_theta
        for j in range(n_phi):
            phi = 2.0 * math.pi * j / n_phi
            index[(i, j)] = len(points)
            points.append((radius * scale * math.cos(theta),
                           radius * scale * math.sin(theta) * math.cos(phi),
                           radius * scale * math.sin(theta) * math.sin(phi)))
    del rings
    faces = []
    for i in range(0, n_theta):
        for j in range(n_phi):
            k = (j + 1) % n_phi
            face = ([0, index[(1, j)], index[(1, k)]] if i == 0 else
                    [index[(i, j)], index[(i + 1, j)], index[(i + 1, k)], index[(i, k)]])
            v = np.array([points[p] for p in face])
            fan = np.cross(v[1:-1] - v[0], v[2:] - v[0]).sum(axis=0)
            toward_source = fan @ v.mean(axis=0) < 0.0
            if toward_source == reverse:
                face = face[::-1]
            faces.append(face)
    n_sphere = len(faces)
    far = len(points)
    points += [(4 * radius, -1, -1), (4 * radius, 1, -1), (4 * radius, 1, 1),
               (4 * radius, -1, 1), (0, 2 * radius, 0), (0, 3 * radius, 0),
               (0, 3 * radius, 1), (0, 2 * radius, 1)]
    faces += [[far, far + 1, far + 2, far + 3], [far + 4, far + 7, far + 6, far + 5]]
    make_case(case, points=points, faces=faces, patches=[
        patch("nozzle", n_sphere, 0),
        patch("vacuum", 1, n_sphere),
        patch("upstreamVacuum", 1, n_sphere + 1)])
    return n_sphere


@pytest.fixture
def meshed_smoke(generated, tmp_path):
    case = clone(generated, tmp_path) / "Interface" / SMOKE
    derived = interface_case.derive(case)
    n_faces = hemisphere_mesh(case, derived.selection.radius_m)
    interface_case.write_identity(derived)
    return case, derived, n_faces


def test_inflow_is_written_per_face_in_patch_order(meshed_smoke):
    case, derived, n_faces = meshed_smoke
    written, document = interface_case.write_inflow(derived)
    n = read_patch_field(case / "0" / "boundaryNumberDensity_Ar", "nozzle")
    u = read_patch_field(case / "0" / "boundaryU", "nozzle")
    t = read_patch_field(case / "0" / "boundaryT", "nozzle")
    assert n.shape == (n_faces,) and u.shape == (n_faces, 3) and t.shape == (n_faces,)
    assert np.all(t > 0.0)
    from plumetools import geometry as facegeom
    from plumetools.inflow import read_patch_geometry
    _, faces, points = read_patch_geometry(case, "nozzle")
    centres = facegeom.centroids(faces, points)
    theta = np.degrees(np.arccos(np.clip(centres[:, 0] / np.linalg.norm(centres, axis=1),
                                         -1, 1)))
    core = theta <= derived.settings.core_half_angle_deg
    assert np.all(n[~core] == 0.0) and np.all(n[core] > 0.0)
    # Face i's value is face i's state: radial velocity along its own centroid.
    r_hat = centres / np.linalg.norm(centres, axis=1)[:, None]
    np.testing.assert_allclose(u / np.linalg.norm(u, axis=1)[:, None], r_hat, atol=1e-9)
    balance = document["flux_balance"]
    assert balance["within_tolerance"] is True or derived.settings.scientific is False
    assert balance["from_written_fields_gross_mass_kg_per_s"] == pytest.approx(
        balance["gross_injected_mass_kg_per_s"], rel=1e-8)
    # The other patches inject nothing.
    assert "nonuniform" not in (case / "0" / "boundaryNumberDensity_Ar").read_text(
        ).split("vacuum")[-1]


def test_a_patch_wound_the_wrong_way_is_refused(generated, tmp_path):
    case = clone(generated, tmp_path) / "Interface" / SMOKE
    derived = interface_case.derive(case)
    hemisphere_mesh(case, derived.selection.radius_m, reverse=True)
    interface_case.write_identity(derived)
    with pytest.raises(policy.ExtensionError, match="oriented back"):
        interface_case.write_inflow(derived)


def test_a_sphere_of_the_wrong_radius_is_refused(generated, tmp_path):
    case = clone(generated, tmp_path) / "Interface" / SMOKE
    derived = interface_case.derive(case)
    hemisphere_mesh(case, derived.selection.radius_m, scale=1.1)
    interface_case.write_identity(derived)
    with pytest.raises(policy.ExtensionError, match="not the sphere"):
        interface_case.write_inflow(derived)


def test_inflow_refuses_a_mesh_built_for_another_configuration(meshed_smoke):
    case, derived, _ = meshed_smoke
    path = case / interface_case.IDENTITY_FILE
    stored = yaml.safe_load(path.read_text())
    stored["R_interface_m"] *= 1.01
    path.write_text(yaml.safe_dump(stored))
    with pytest.raises(policy.ExtensionError, match="R_interface_m"):
        interface_case.write_inflow(derived)


# --------------------------------------------------------------------------- #
# the PBS steps
# --------------------------------------------------------------------------- #

def test_guard_allows_kn_001_only_for_validation(generated, tmp_path):
    family = clone(generated, tmp_path)
    case = family / "Interface" / "val_Kn0p01_B0p05"
    problem = hpc.guard(case)
    # Not the Kn rule: the interface itself is what is refused.
    assert problem is not None and "NO_ADMISSIBLE_INTERFACE" in problem
    raw = yaml.safe_load((case / "interface.yaml").read_text())
    raw["interface"]["purpose"] = "production"
    (case / "interface.yaml").write_text(yaml.safe_dump(raw))
    assert "Kn = 0.01 is one of Cai's" in hpc.guard(case)


def test_meshing_running_and_labelling_are_three_decisions(generated, tmp_path):
    """The committed production interface: the guard lets it be meshed (the
    sphere is admissible and planned), the label refuses to run it as
    production (numerically unsuitable, coupling pending), and as an
    exploratory run it is labelled NON-PRODUCTION."""
    family = clone(generated, tmp_path)
    case = family / "Interface" / "Kn0p001_np1x_sphere"
    assert hpc.guard(case) is None
    with pytest.raises(policy.ExtensionError, match="REFUSED"):
        hpc.label(case)
    assert not (case / hpc.MARKER).exists()
    raw = yaml.safe_load((case / "interface.yaml").read_text())
    raw["interface"]["run_mode"] = "exploratory"
    (case / "interface.yaml").write_text(yaml.safe_dump(raw))
    assert hpc.guard(case) is None
    assert "unsuitable-for-production" in hpc.label(case)
    assert (case / hpc.MARKER).is_file()


def test_guard_refuses_an_incomplete_production_case(generated, tmp_path):
    family = clone(generated, tmp_path)
    case = family / "Interface" / "Kn0p001_np1x_sphere"
    raw = yaml.safe_load((case / "interface.yaml").read_text())
    raw["interface"]["core_half_angle_deg"] = None
    (case / "interface.yaml").write_text(yaml.safe_dump(raw))
    problem = hpc.guard(case)
    assert problem and "CONFIGURATION_INCOMPLETE" in problem
    assert interface_case.derive(case).selection.status == interface.STATUS_INCOMPLETE


def test_guard_refuses_results_from_another_mesh(meshed_smoke):
    case, derived, _ = meshed_smoke
    assert hpc.guard(case) is None
    (case / "0.001").mkdir()
    assert hpc.guard(case) is None
    (case / interface_case.IDENTITY_FILE).unlink()
    assert "never resumed" in hpc.guard(case)


def test_label_marks_smoke_runs_and_refuses_pending_production(generated, tmp_path):
    family = clone(generated, tmp_path)
    smoke = family / "Interface" / SMOKE
    text = hpc.label(smoke)
    assert (smoke / hpc.MARKER).is_file() and "development-smoke" in text
    production = family / "Interface" / "Kn0p001_np1x_sphere"
    raw = yaml.safe_load((production / "interface.yaml").read_text())
    raw["interface"].update(core_half_angle_deg=30.0,
                            source_validity_min_radius_over_D=1.0)
    (production / "interface.yaml").write_text(yaml.safe_dump(raw))
    assessment = interface_case.assess(interface_case.derive(production))
    assert assessment["coupling_validation"]["status"] == "pending"
    assert not assessment["runnable"]
    with pytest.raises(policy.ExtensionError, match="REFUSED"):
        hpc.label(production)
    assert not (production / hpc.MARKER).exists()
    raw["interface"]["run_mode"] = "exploratory"
    (production / "interface.yaml").write_text(yaml.safe_dump(raw))
    assert "exploratory" in hpc.label(production)
    assert (production / hpc.MARKER).is_file()


def test_production_needs_passed_evidence_at_the_same_threshold(generated, tmp_path):
    family = clone(generated, tmp_path)
    case = family / "Interface" / "Kn0p001_np1x_sphere"
    raw = yaml.safe_load((case / "interface.yaml").read_text())
    raw["interface"].update(core_half_angle_deg=30.0,
                            source_validity_min_radius_over_D=1.0,
                            validation_evidence="overlap.yaml")
    (case / "interface.yaml").write_text(yaml.safe_dump(raw))
    for verdict, threshold, expected in (("pass", 0.05, "passed"),
                                         ("pass", 0.03, "threshold-mismatch"),
                                         ("inconclusive", 0.05, "inconclusive")):
        (case / "overlap.yaml").write_text(yaml.safe_dump(
            {"verdict": verdict, "interface": {"bird_threshold": threshold}}))
        coupling = interface_case.assess(interface_case.derive(case))["coupling_validation"]
        assert coupling["status"] == expected


def test_the_direct_exit_post_context_is_the_parent_derivation(generated, tmp_path):
    """A case with no interface.yaml is post-processed exactly as before."""
    case = tmp_path / "direct"
    shutil.copytree(FAMILY / "baseCase", case,
                    ignore=shutil.ignore_patterns("system", "constant", "results"))
    context = interface_case.post_context(case)
    derived = generation.derive(case)
    assert context.mode == "direct" and context.interface is None
    assert context.run == derived.run
    assert context.nozzle_area_error == derived.nozzle.area_error


@pytest.mark.skipif(sys.platform == "win32",
                    reason="bash on native Windows is WSL's and cannot open a "
                           "Windows path; the suite runs under WSL")
def test_overlay_scripts_parse():
    if shutil.which("bash") is None:
        pytest.skip("no bash")
    for name in ("Allmesh", "runPost.sh", "runViz.sh"):
        result = subprocess.run(["bash", "-n", str(FAMILY / "baseCaseInterface" / name)],
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr


def test_interface_mesh_patch_problems_are_specific():
    report = {"vertex_radius_over_R": [0.97, 1.0], "fitted_radius_over_R": 1.0,
              "fitted_centre_over_D": [0.0, 0.0, 0.0], "min_x_over_D": 0.0,
              "faces_facing_the_source": 0, "faces_in_core": 3, "area_error": 0.0}
    problems = interface_mesh.patch_problems(report, snap_tolerance=0.02,
                                             min_core_faces=12, area_tolerance=0.05)
    assert len(problems) == 2


# --------------------------------------------------------------------------- #
# the derivation's shape, and the module boundaries
# --------------------------------------------------------------------------- #

def test_planning_is_all_or_nothing(generated):
    """A case is meshable exactly when every planned quantity exists."""
    smoke = interface_case.derive(generated / "Interface" / SMOKE)
    assert smoke.meshable and isinstance(smoke.planned, interface_case.InterfacePlan)
    assert smoke.mesh_plan is smoke.planned.mesh_plan and smoke.run is smoke.planned.run
    infeasible = interface_case.derive(generated / "Interface" / "val_Kn0p01_B0p05")
    assert not infeasible.meshable and infeasible.planned is None
    assert (infeasible.mesh_plan, infeasible.run, infeasible.parcels) == (None, None, None)
