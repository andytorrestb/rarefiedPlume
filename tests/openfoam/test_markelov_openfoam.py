"""The AIAA 99-3455 case against a real OpenFOAM installation.

Marked ``needs_openfoam`` and deselected by default. These answer the questions
the pure-Python tests cannot:

* do the generated dictionaries parse, and does the pipeline complete?
* is the mesh the geometry the config declares -- measured from the mesh?
* does the custom inflow library load, and does it select only the named patch?
* does it read a spatially varying number density, rather than collapsing it?
* do particles leave through an open ``patch`` boundary instead of reflecting?
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from conftest import PLUME_LIB, plume_library_path, run, run_or_fail

pytestmark = pytest.mark.needs_openfoam


# --------------------------------------------------------------------------- #
# dictionaries and meshing
# --------------------------------------------------------------------------- #

def test_generated_dictionaries_parse(meshed_case):
    """foamDictionary reads each one, which is a stricter test than "blockMesh
    did not crash" -- a dictionary can be malformed in a section nothing used."""
    if not run(["foamDictionary", "-help"], meshed_case).returncode == 0:
        pytest.skip("foamDictionary is not on PATH")

    for relative in ("system/blockMeshDict", "system/snappyHexMeshDict",
                     "system/meshQualityDict", "system/controlDict",
                     "system/dsmcInitialiseDict", "system/decomposeParDict",
                     "constant/dsmcProperties"):
        assert (meshed_case / relative).is_file(), f"{relative} was not generated"
        run_or_fail(["foamDictionary", relative], meshed_case)


def test_block_mesh_and_snappy_completed(meshed_case):
    assert (meshed_case / "log.blockMesh").is_file()
    assert (meshed_case / "log.snappyHexMesh").is_file()
    assert (meshed_case / "constant" / "polyMesh" / "points").is_file()
    assert "End" in (meshed_case / "log.snappyHexMesh").read_text(errors="replace")


def test_check_mesh_passes(meshed_case):
    log = (meshed_case / "log.checkMesh").read_text(errors="replace")
    assert "Mesh OK." in log, "checkMesh reported failures"
    assert "***" not in log or "Mesh OK." in log


def test_required_patch_names_exist_with_the_right_types(meshed_case):
    from plumetools.mesh.boundary import read_boundary

    patches = read_boundary(meshed_case)
    expected = {
        "inflow": "patch",
        "cylinder": "wall",
        "plate": "wall",
        "vacuum": "patch",
        "symmetry": "symmetry",
        "upstreamVacuum": "patch",
    }
    assert set(expected) <= set(patches), f"mesh has {sorted(patches)}"
    for name, geometric_type in expected.items():
        assert patches[name].type == geometric_type
        assert patches[name].n_faces > 0, f"{name} is empty"


def test_mesh_dimensions_agree_with_the_configuration(meshed_case):
    """plumetools.markelov1999.verify measures every dimension from the mesh's own
    vertices, so this is not a check that the generator is self-consistent."""
    from plumetools.markelov1999.verify import verify

    report = verify(meshed_case)
    assert any("mesh verification" not in line for line in report)
    assert not any(line.startswith("  FAIL") for line in report)


def test_verify_mesh_passes_as_a_command(meshed_case):
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-m", "plumetools.markelov1999.verify", str(meshed_case)],
        capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "mesh verification passed" in result.stdout


def test_the_six_inch_gap_is_measured_between_the_two_patches(meshed_case):
    """The paper's defining dimension, read out of the mesh rather than the
    config: the cylinder's downstream extent from the cylinder patch's vertices,
    the plate's upstream face from the plate patch's."""
    from plumetools.inflow import read_patch_geometry

    def patch_points(name):
        _labels, faces, points = read_patch_geometry(meshed_case, name)
        return points[sorted({p for f in faces for p in f})]

    cylinder_base = patch_points("cylinder")[:, 0].max()
    plate_face = patch_points("plate")[:, 0].min()

    gap = plate_face - cylinder_base
    assert gap == pytest.approx(0.1524, rel=1e-6)
    assert gap / 0.0254 == pytest.approx(6.0, rel=1e-6)


# --------------------------------------------------------------------------- #
# the custom inflow library
# --------------------------------------------------------------------------- #

def test_the_custom_library_is_built():
    path = plume_library_path()
    if path is None:
        pytest.skip(f"{PLUME_LIB} is not built; "
                    f"run applications/dsmcBoundaryModels/Allwmake")
    assert path.is_file()


def test_the_case_selects_the_custom_inflow_model(meshed_case):
    properties = (meshed_case / "constant" / "dsmcProperties").read_text()
    assert re.search(r"InflowBoundaryModel\s+plumeFieldInflow;", properties)
    assert "plumeFieldInflowCoeffs" in properties
    assert re.search(r"patches\s*\(\s*inflow\s*\)", properties)

    control = (meshed_case / "system" / "controlDict").read_text()
    assert PLUME_LIB in control, "controlDict must load the library"


def test_the_library_loads_and_injects_on_the_named_patch_only(solved_case):
    """The whole point of the model: FreeStream would inject on every patch-type
    boundary, which here would make the vacuum boundary a second source."""
    log = (solved_case / "log.dsmcFoam").read_text(errors="replace")

    assert "Selecting InflowBoundaryModel plumeFieldInflow" in log, (
        "the custom model was not selected; is the library on the library path?")
    assert re.search(r"plumeFieldInflow: injecting across 1\(inflow\)", log), (
        "the model must inject across exactly one patch, named inflow")
    assert "boundaryNumberDensity_N2" in log


def test_a_missing_library_fails_rather_than_falling_back(meshed_case, tmp_path):
    """Requirement 7: never silently fall back to uniform inflow. Allrun refuses
    to start when the library is absent, before spending a decomposePar on it."""
    import os
    import shutil
    import subprocess

    case = tmp_path / "nolib"
    shutil.copytree(meshed_case, case)
    for script in case.glob("All*"):
        script.chmod(0o755)

    env = dict(os.environ)
    # Point every library search variable somewhere empty.
    for variable in ("FOAM_USER_LIBBIN", "FOAM_SITE_LIBBIN", "FOAM_LIBBIN"):
        env[variable] = str(tmp_path / "empty")
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[2])

    result = subprocess.run(
        ["bash", "./Allrun", "--serial"], cwd=str(case), capture_output=True,
        text=True, timeout=600, env=env, check=False)

    assert result.returncode != 0, "a missing library must stop the run"
    assert "is not built" in result.stderr
    assert "no fallback" in result.stderr or "cannot express" in result.stderr


# --------------------------------------------------------------------------- #
# spatially varying inflow
# --------------------------------------------------------------------------- #

def test_the_inflow_density_field_is_spatially_varying(solved_case):
    """A uniform value would mean the angular structure the model exists to
    compute had been collapsed -- which is exactly what stock FreeStream forces."""
    from plumetools.markelov1999.foamfields import read_patch_field

    values = read_patch_field(solved_case / "0" / "boundaryNumberDensity_N2", "inflow")
    assert len(values) > 100, "the inflow patch should have many faces"
    assert values.min() > 0.0
    assert values.max() / values.min() > 5.0, (
        f"density spans only {values.max() / values.min():.2f}x; "
        f"the source-flow angular profile should vary far more than that")


def test_the_solver_accepts_the_varying_field_and_inserts_particles(solved_case):
    log = (solved_case / "log.dsmcFoam").read_text(errors="replace")
    inserted = [int(m) for m in re.findall(r"Particles inserted\s*=\s*(\d+)", log)]

    assert inserted, "no insertion was reported"
    assert all(n > 0 for n in inserted), "some step inserted nothing"
    # Steady injection: the boundary condition is time-independent, so the count
    # should fluctuate about a mean rather than drift or collapse.
    assert max(inserted) / min(inserted) < 2.0


def test_the_run_reaches_the_end_time(solved_case):
    log = (solved_case / "log.dsmcFoam").read_text(errors="replace")
    assert "End" in log
    counts = [int(m) for m in re.findall(
        r"Number of dsmc particles\s*=\s*(\d+)", log)]
    assert counts and counts[-1] > counts[0], "the domain should be filling"


# --------------------------------------------------------------------------- #
# vacuum outflow
# --------------------------------------------------------------------------- #

def test_particles_leave_through_the_open_boundaries(solved_case):
    """Requirement 7 asks whether an ordinary open patch already deletes outgoing
    particles once unwanted inflow is removed. It does: DSMCParcel::hitPatch
    returns false, and particle::hitBoundaryFace falls through its dispatch chain
    to keepParticle = false for anything that is not a wall, wedge, symmetry,
    cyclic or processor patch.

    Measured rather than assumed: with a reflecting boundary the particle count
    would equal everything ever put into the domain. It is lower, so particles
    are being removed. (The transit time is far longer than this short run, so
    most injected particles are still in flight -- the deficit is small but it is
    not zero, and with reflection it would be exactly zero.)
    """
    log = (solved_case / "log.dsmcFoam").read_text(errors="replace")

    inserted = sum(int(m) for m in re.findall(
        r"Particles inserted\s*=\s*(\d+)", log))
    counts = [int(m) for m in re.findall(
        r"Number of dsmc particles\s*=\s*(\d+)", log)]
    initial = (meshed_initial(solved_case) or 0)

    present = counts[-1]
    supplied = inserted + initial

    assert present < supplied, (
        f"{present} particles present but {supplied} were supplied "
        f"({inserted} injected + {initial} initialised). With a reflecting outer "
        f"boundary these would be equal.")


def meshed_initial(case: Path) -> int | None:
    """Particles created by dsmcInitialise, from its log."""
    log_path = case / "log.dsmcInitialise"
    if not log_path.is_file():
        return None
    match = re.search(r"Total number of molecules added:\s*(\d+)",
                      log_path.read_text(errors="replace"))
    return int(match.group(1)) if match else None


def test_no_particle_is_outside_the_domain(solved_case):
    """A reflecting or leaking boundary shows up here. Read from the written
    lagrangian positions."""
    positions = solved_case.glob("*/lagrangian/dsmc/positions")
    path = next((p for p in positions if p.stat().st_size > 0), None)
    if path is None:
        pytest.skip("no lagrangian positions were written")

    text = path.read_text(errors="replace")
    coordinates = np.array(
        [[float(v) for v in m] for m in re.findall(
            r"\(\s*(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*\)", text)])
    if len(coordinates) == 0:
        pytest.skip("could not parse the positions file")

    # The domain, with a tolerance of one background cell for barycentric storage.
    tolerance = 0.03
    assert coordinates[:, 0].min() > -tolerance
    assert coordinates[:, 0].max() < 0.9 + tolerance
    assert coordinates[:, 1].min() > -tolerance
    assert coordinates[:, 1].max() < 0.3 + tolerance
    assert abs(coordinates[:, 2]).max() < 0.35 + tolerance


# --------------------------------------------------------------------------- #
# post-processing
# --------------------------------------------------------------------------- #

def test_the_surface_pressure_field_is_written_with_values(solved_case):
    """fD and q must be 'calculated' on wall patches: hitWallPatch writes into
    their boundary values, and a zeroGradient patch field writes none at all --
    so the surface pressure would be discarded at write time."""
    from plumetools.markelov1999.foamfields import read_patch_field
    from plumetools.markelov1999.postprocess import latest_time_dir

    time_dir = latest_time_dir(solved_case)
    assert time_dir is not None

    for field in ("fD", "fDMean", "q"):
        path = time_dir / field
        if not path.is_file():
            continue
        values = read_patch_field(path, "cylinder")   # must not raise
        assert len(values) > 0


def test_post_processing_produces_a_case_summary(solved_case):
    import json
    import os
    import subprocess
    import sys

    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[2])
    result = subprocess.run(
        [sys.executable, "postProcess.py", "."], cwd=str(solved_case),
        capture_output=True, text=True, timeout=1800, env=env, check=False)
    assert result.returncode == 0, result.stdout + result.stderr

    summary = json.loads(
        (solved_case / "case-summary.json").read_text(encoding="utf-8"))
    assert summary["mesh"]["n_cells"] > 0
    assert summary["mesh"]["check_mesh_ok"] is True
    assert summary["particles"]["n_time_steps"] > 0
    assert "resolution_audit" in summary
    assert "flux" in summary, "runInflow.py's flux verification should survive"


def test_the_flux_verification_is_recorded_and_consistent(solved_case):
    """The analytical, quadrature and field-reconstructed fluxes must agree."""
    import json

    summary = json.loads(
        (solved_case / "case-summary.json").read_text(encoding="utf-8"))
    flux = summary["flux"]

    for name, error in flux["relative_errors"].items():
        assert error < 1e-6, f"{name} flux differs from the fields by {error:.2e}"
    for name, error in flux["quadrature_agreement"].items():
        assert error < 1e-6, f"{name} quadrature disagrees by {error:.2e}"

    assert flux["analytical"]["mass_kg_per_s"] > 0.0
    assert flux["analytical"]["energy_W"] > 0.0
    assert flux["analytical"]["momentum_x_N"] > 0.0


def test_the_dsmc_checks_all_passed(solved_case):
    import json

    summary = json.loads(
        (solved_case / "case-summary.json").read_text(encoding="utf-8"))
    assert summary["checks"]["passed"] is True
    assert summary["checks"]["n_failures"] == 0
