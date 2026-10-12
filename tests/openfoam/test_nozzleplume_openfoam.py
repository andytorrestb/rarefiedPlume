"""Tier 2: the nozzle-plume family against the real OpenFOAM toolchain.

Marked ``needs_openfoam`` and deselected by default::

    pytest -m "not needs_openfoam"                      # the fast suite
    pytest tests/openfoam -m needs_openfoam             # these

What the Tier-0 tests cannot answer
-----------------------------------
Everything here is a question about OpenFOAM, not about this repository's
Python:

* does ``snappyHexMesh`` carve the plate the plan describes, exactly, and leave
  the exit plane at the one refinement level the plan predicts?
* does ``topoSet`` select the nozzle faces the Python rule says it will?
* does ``plumeFieldInflow`` inject at the Maxwellian flux, on the nozzle only?
* does the plate feel the load it must?

The last is answered by a case with a known result. The jet is run
**collisionless** onto a diffuse plate large enough to catch all of it. Then
nothing but conservation fixes the totals:

.. code-block:: text

    force  =  jet momentum flux  +  mass flow * sqrt(pi R T_wall / 2)
    heat   =  jet energy flux    -  molecule flow * 2 k T_wall

-- what arrives, plus the recoil of diffuse re-emission; what arrives, less what
re-emission carries off. Getting both right needs the inflow flux, the wall
model's temperature, the ``fD`` and ``q`` accumulation, ``fieldAverage``, the
symmetry multiplicity and the post-processing normals all to be right at once,
which is why it is the test.

The case is the shipped ``cases/r-4d/baseCase`` -- its own ``Allmesh``,
``Allrun`` and ``Allpost`` -- with ``case.yaml`` shrunk to a quarter-symmetric
3 m box. The gas and the exit state are the R-4D's, untouched.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess

import numpy as np
import pytest
import yaml

from conftest import PLUME_LIB, REPO, SOLVER_COMMANDS, missing, plume_library_path
from plumetools.mesh.boundary import read_boundary
from plumetools.nozzleplume import exitstate, geometry, mesh, meshread, post
from plumetools.nozzleplume.config import load_case_config
from plumetools.nozzleplume.constants import BOLTZMANN_J_PER_K

pytestmark = pytest.mark.needs_openfoam

R4D_BASE_CASE = REPO / "cases" / "r-4d" / "baseCase"

#: The validation variant: what changes from the shipped case.yaml.
#:
#: A 10 x 10 m plate 2 m from the exit subtends 68 degrees from the axis, and a
#: speed-ratio-5 jet has nothing out there, so the plate catches all of it. The
#: slowest molecules that matter (U - 4 sigma) reach it in 1.4 ms, so averaging
#: from 1.5 ms sees a steady load. Collisionless, the gas that leaves the plate
#: never returns to it, so nothing else has to settle.
VALIDATION = {
    "meta": {"case_name": "nozzleplume-validation"},
    "plate": {"standoff_m": 2.0, "size_y_m": 10.0, "size_z_m": 10.0},
    "domain": {"x_max_m": 3.0, "y_half_m": 6.0, "z_half_m": 6.0,
               "symmetry_planes": ["y", "z"]},
    "mesh": {"refinement_regions": [
        {"name": "exitCore", "level": 4, "x_max_m": 1.0, "radius_m": 0.5},
        {"name": "nearField", "level": 3, "x_max_m": 2.0, "radius_m": 1.0},
    ]},
    "dsmc": {
        "collision_model": "NoBinaryCollision",
        "n_equivalent_particles": 1.0e+17,
        "end_time_s": 3.5e-3,
        "average_start_s": 1.5e-3,
    },
}


#: What the shipped case IS: its inputs and its scripts. Anything else in the
#: directory is something a run left there.
CASE_SOURCES = ("case.yaml", "viz.yaml", "open.foam", "All*")


def _write_case(destination):
    """Copy baseCase's own files and apply the validation variant to case.yaml.

    The files, not the directory. The shipped case is meant to be run where it
    is, and one that has been holds a mesh, time directories and processor
    directories -- gigabytes of them. A copy of that is a case ./Allmesh finds
    already built, so these tests would pass only on a checkout nobody had used.
    """
    destination.mkdir(parents=True)
    for pattern in CASE_SOURCES:
        for source in sorted(R4D_BASE_CASE.glob(pattern)):
            if source.is_file():
                shutil.copy2(source, destination / source.name)
    for script in destination.glob("All*"):
        script.chmod(0o755)

    path = destination / "case.yaml"
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    for section, values in VALIDATION.items():
        data[section].update(values)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return destination


def _env():
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
    return env


def _script(name, case, *args, timeout=7200):
    result = subprocess.run(
        ["bash", f"./{name}", *args], cwd=str(case), capture_output=True,
        text=True, timeout=timeout, env=_env(), check=False)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-40:])
        pytest.fail(f"./{name} failed:\n{tail}")
    return result


def _plan(case):
    cfg = load_case_config(case)
    geom = geometry.from_config(cfg)
    exit_state = exitstate.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    return cfg, geom, exit_state, plan, mesh.predicted_nozzle_faces(plan, geom)


@pytest.fixture(scope="module")
def meshed(openfoam, tmp_path_factory):
    """The validation case, meshed by its own ./Allmesh."""
    absent = missing(("topoSet", "createPatch"))
    if absent:
        pytest.skip(f"mesh utilities not on PATH: {absent}")
    case = _write_case(tmp_path_factory.mktemp("nozzleplume") / "case")
    _script("Allmesh", case)
    return case


@pytest.fixture(scope="module")
def solved(meshed, tmp_path_factory):
    """A run of the meshed case to a steady plate load, post-processed."""
    absent = missing(SOLVER_COMMANDS)
    if absent:
        pytest.skip(f"dsmcFoam utilities not on PATH: {absent}")
    if plume_library_path() is None:
        pytest.skip(
            f"{PLUME_LIB} is not built; run applications/dsmcBoundaryModels/Allwmake")

    case = tmp_path_factory.mktemp("nozzleplume-run") / "case"
    shutil.copytree(meshed, case)
    for script in case.glob("All*"):
        script.chmod(0o755)
    _script("Allrun", case, "--serial")
    _script("Allpost", case)
    return case


# --------------------------------------------------------------------------- #
# the mesh OpenFOAM built is the mesh that was planned
# --------------------------------------------------------------------------- #

def test_allmesh_ends_by_verifying_the_mesh(meshed):
    """./Allmesh ran `verify` and would have failed the fixture had it not
    passed; this asserts the pipeline got that far for the right reason."""
    log = (meshed / "log.checkMesh").read_text(encoding="utf-8")
    assert "Mesh OK." in log
    assert (meshed / "system" / "snappyHexMeshDict").is_file()


def test_the_patches_have_the_types_the_physics_needs(meshed):
    patches = read_boundary(meshed)
    assert {name: info.type for name, info in patches.items()} == {
        "upstreamVacuum": "patch", "symmetry": "symmetry", "vacuum": "patch",
        "plate": "wall", "nozzle": "patch"}


def test_toposet_selects_the_faces_the_plan_predicts(meshed):
    """The Python rule and cylinderToFace are the same rule, applied to the
    lattice snappyHexMesh leaves on the exit plane."""
    *_, nozzle = _plan(meshed)
    patch = meshread.read_patches(meshed, ["nozzle"])["nozzle"]
    assert patch.n_faces == nozzle.n_faces == 28           # a quarter of 112
    # To a part in a million, not to round-off: where the nozzle meets a
    # symmetry plane, snappyHexMesh's smoothing leaves those vertices a few
    # nanometres off it (measured: 4e-8 of the area).
    assert patch.area_m2 == pytest.approx(nozzle.area_m2, rel=1e-6)
    assert np.all(np.abs(patch.vertices[:, 0]) < 1e-9)
    assert np.all(patch.normals[:, 0] < -0.999999)


def test_snappy_carves_the_plate_exactly(meshed):
    """Feature snapping: no chamfered edge, and an upstream face whose area is
    the product of its two sides."""
    _, geom, *_ = _plan(meshed)
    plate = meshread.read_patches(meshed, ["plate"])["plate"]
    x0, x1, y0, y1, z0, z1 = geom.modelled_plate_bounds

    lo, hi = plate.vertices.min(axis=0), plate.vertices.max(axis=0)
    assert lo == pytest.approx([x0, y0, z0], abs=1e-6)
    assert hi == pytest.approx([x1, y1, z1], abs=1e-6)

    front = plate.normals[:, 0] > 0.999
    assert plate.areas[front].sum() == pytest.approx((y1 - y0) * (z1 - z0), rel=1e-6)
    assert not np.any(np.abs(plate.normals).max(axis=1) < 0.999)   # no oblique face


def test_the_cell_estimate_is_close_to_the_mesh(meshed):
    *_, plan, _ = _plan(meshed)
    assert meshread.read_cell_count(meshed) == pytest.approx(
        plan.estimated_cells, rel=0.10)


# --------------------------------------------------------------------------- #
# the solver does what the dictionaries ask
# --------------------------------------------------------------------------- #

def test_inflow_is_on_the_nozzle_and_nowhere_else(solved):
    log = (solved / "log.dsmcFoam").read_text(encoding="utf-8", errors="replace")
    assert "plumeFieldInflow: injecting across 1(nozzle)" in log
    assert "Selecting BinaryCollisionModel none" in log


def test_the_injection_rate_is_the_maxwellian_flux(solved):
    """Per step: flux * meshed area * deltaT / weight."""
    cfg, _, exit_state, _, nozzle = _plan(solved)
    log = (solved / "log.dsmcFoam").read_text(encoding="utf-8", errors="replace")
    inserted = np.array([int(m) for m in re.findall(
        r"Particles inserted\s*=\s*(\d+)", log)])
    delta_t = float(re.search(
        r"^deltaT\s+([^;]+);", (solved / "system" / "controlDict").read_text(),
        re.M).group(1))
    expected = (exit_state.number_flux_per_m2_s * nozzle.area_m2 * delta_t
                / cfg.dsmc.n_equivalent_particles)

    assert inserted.size > 2000
    assert inserted.mean() == pytest.approx(expected, rel=0.01)


def test_the_open_boundaries_delete_molecules(solved):
    """Were an open boundary a wall, the parcel count would be every parcel ever
    injected. It is not -- but only just, this early: the plate re-emits at the
    wall's thermal speed, a tenth of the jet's, and by 3.5 ms only the fast tail
    of that has crossed the 2 m back to the exit plane. About 2% have gone."""
    log = (solved / "log.dsmcFoam").read_text(encoding="utf-8", errors="replace")
    inserted = sum(int(m) for m in re.findall(r"Particles inserted\s*=\s*(\d+)", log))
    final = int(re.findall(r"Number of dsmc particles\s*=\s*(\d+)", log)[-1])
    assert final < 0.99 * inserted
    assert final > 0.90 * inserted


# --------------------------------------------------------------------------- #
# the plate feels what conservation says it must
# --------------------------------------------------------------------------- #

def _expected_loads(case):
    """Force and heat on a diffuse plate that catches a collisionless jet."""
    _, geom, exit_state, _, nozzle = _plan(case)
    area = nozzle.area_m2 * geom.symmetry_multiplicity      # the MESHED exit
    wall = geom.plate.temperature_K
    r = exit_state.species.gas_constant_j_per_kg_k

    force = (exit_state.thrust_n(area)
             + exit_state.mass_flow_kg_per_s(area) * math.sqrt(math.pi * r * wall / 2.0))
    heat = (exit_state.power_w(area)
            - exit_state.injection_rate_per_s(area) * 2.0 * BOLTZMANN_J_PER_K * wall)
    return force, heat


def test_the_plate_force_is_the_jet_momentum_plus_the_diffuse_recoil(solved):
    cfg, geom, *_ = _plan(solved)
    force, _ = _expected_loads(solved)
    loads = post.plate_loads(solved, cfg, geom)
    assert loads.force_n[0] == pytest.approx(force, rel=0.01)


def test_the_plate_heat_load_is_the_jet_power_less_the_re_emission(solved):
    cfg, geom, *_ = _plan(solved)
    _, heat = _expected_loads(solved)
    loads = post.plate_loads(solved, cfg, geom)
    assert loads.heat_load_w == pytest.approx(heat, rel=0.01)


def test_the_load_is_on_the_upstream_face_and_peaks_on_the_axis(solved):
    cfg, geom, *_ = _plan(solved)
    loads = post.plate_loads(solved, cfg, geom)
    front = loads.side == 1
    areas = loads.areas

    on_front = (loads.force_density[front, 0] * areas[front]).sum()
    assert on_front == pytest.approx(loads.modelled_force_n[0], rel=1e-3)

    peak = loads.centres[int(np.argmax(np.where(front, loads.pressure, -np.inf)))]
    assert math.hypot(peak[1], peak[2]) < 0.2            # within a cell of the axis


def test_the_two_lateral_forces_on_the_quarter_are_equal(solved):
    """The jet is axisymmetric, so the quarter plate is pushed outwards equally
    along y and z -- and the whole plate, by symmetry, along neither."""
    cfg, geom, *_ = _plan(solved)
    loads = post.plate_loads(solved, cfg, geom)
    _, fy, fz = loads.modelled_force_n
    assert fy > 0.0 and fz > 0.0
    assert fy == pytest.approx(fz, rel=0.05)
    assert loads.force_n[1] == loads.force_n[2] == 0.0


def test_allpost_writes_its_results(solved):
    results = solved / "results"
    for name in ("plate_loads.csv", "plate_summary.yaml", "centreline.csv"):
        assert (results / name).stat().st_size > 0, name
    summary = yaml.safe_load((results / "plate_summary.yaml").read_text())
    assert summary["symmetry_multiplicity"] == 4
    assert 1.05 < summary["axial_force_over_jet_thrust"] < 1.15


def test_the_centreline_leaves_the_exit_at_the_exit_state(solved):
    """Just off the nozzle the gas is the injected one. (Further out the plate's
    re-emitted molecules stream back through it, collisionless, and the moments
    are no longer those of a single beam.)"""
    cfg, geom, exit_state, *_ = _plan(solved)
    profile = post.centreline(solved, cfg, geom, exit_state)
    near = profile["x_m"] < 0.3
    assert near.sum() >= 3
    assert profile["U_over_U_exit"][near] == pytest.approx(1.0, abs=0.03)
    assert profile["n_over_n_exit"][near][1:] == pytest.approx(1.0, abs=0.08)


def test_rerunning_a_finished_case_is_a_no_op(solved):
    """./Allrun continues a case; one that has reached endTime has nothing to do
    and must not re-initialise it."""
    before = sorted(p.name for p in solved.iterdir())
    result = _script("Allrun", solved, "--serial")
    assert "already reached endTime" in result.stdout
    assert sorted(p.name for p in solved.iterdir()) == before
