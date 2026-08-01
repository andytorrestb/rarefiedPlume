"""Bit-exact regression of the plumetools inflow path against the legacy code.

These goldens prove the refactor did not change behaviour. They do NOT prove the
behaviour is scientifically correct: findings SM-01 through SM-10 -- the ignored
per-case stagnation pressure, the split stagnation temperature, the N2/Ar species
mismatch, the theta-from-+z convention and the un-normalised angular function --
are all deliberately frozen into this snapshot. Agreement means "unchanged",
never "right". See tests/regression/README.md and docs/source-flow-model.md.

Runs entirely on the committed Pointwise mesh. No OpenFOAM, no subprocess.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "cases" / "3d-inflow"
GOLDEN = Path(__file__).resolve().parent / "golden"

RTOL = 1e-12

pytestmark = pytest.mark.skipif(
    not (GOLDEN / "3d_inflow_v0.npz").exists(),
    reason="golden not captured; run tests/regression/capture_golden.py",
)


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def golden():
    with np.load(GOLDEN / "3d_inflow_v0.npz") as data:
        return {
            "labels": [str(x) for x in data["labels"]],
            "rhoN": data["rhoN"],
            "U": data["U"],
            "T": data["T"],
        }


@pytest.fixture(scope="module")
def computed():
    """Inflow data from the extracted plumetools path.

    Computed with process spawning blocked: the extraction must read
    constant/polyMesh/boundary directly, where the legacy readMeshStats() shelled
    out to `checkMesh` and parsed its log (findings RB-01, RB-05).
    """
    import os
    import subprocess

    from plumetools.config import load_case_config
    from plumetools.inflow import compute_inflow

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "the inflow path must not spawn a subprocess; it reads "
            "constant/polyMesh/boundary directly"
        )

    with pytest.MonkeyPatch.context() as mp:
        for target, name in (
            (subprocess, "run"), (subprocess, "check_output"),
            (subprocess, "Popen"), (os, "system"),
        ):
            mp.setattr(target, name, forbidden)
        return compute_inflow(CASE, load_case_config(CASE))


# --------------------------------------------------------------------------- #
# 1. golden inflow data
# --------------------------------------------------------------------------- #

def test_face_order_matches_golden(computed, golden):
    """Face order is the contract with dsmcFoam+, not an implementation detail."""
    assert computed.labels == golden["labels"]


def test_face_count(computed):
    assert len(computed.labels) == 2044


def test_rhoN_matches_golden(computed, golden):
    np.testing.assert_allclose(computed.rhoN, golden["rhoN"], rtol=RTOL, atol=0)


def test_velocity_matches_golden(computed, golden):
    np.testing.assert_allclose(computed.U, golden["U"], rtol=RTOL, atol=0)


def test_temperature_matches_golden(computed, golden):
    np.testing.assert_allclose(computed.T, golden["T"], rtol=RTOL, atol=0)


# --------------------------------------------------------------------------- #
# 2. generated field files
# --------------------------------------------------------------------------- #

HEADER_LINES = 16  # banner + FoamFile block + separator; see README

NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _normalise(text: str) -> list[str]:
    return text.replace("\r\n", "\n").replace("\r", "\n").split("\n")


def _split_body(lines: list[str]) -> tuple[list[str], list[str]]:
    """Split a field file's body into (structural lines, numeric-payload lines).

    The payload is everything between the '(' that opens the inflow value list and
    its closing ');'.
    """
    start = next(i for i, l in enumerate(lines) if l.strip() == "(")
    end = next(i for i, l in enumerate(lines[start:], start) if l.strip() == ");")
    return lines[:start + 1] + lines[end:], lines[start + 1:end]


@pytest.fixture(scope="module")
def written_dir(tmp_path_factory, computed):
    """Write the three 0/ field files into a scratch dir, leaving the case clean."""
    from plumetools.config import load_case_config
    from plumetools.foamio.fields import write_inflow_fields

    out = tmp_path_factory.mktemp("zero")
    write_inflow_fields(out, computed, load_case_config(CASE))
    return out


@pytest.fixture(params=["boundaryU", "boundaryT", "boundaryNumberDensity_Ar"])
def field_file(request, written_dir):
    """(name, golden lines, written lines) for one generated 0/ field file."""
    name = request.param
    written = _normalise((written_dir / name).read_text())
    expected = _normalise((GOLDEN / f"{name}.txt").read_text())
    return name, expected, written


def test_field_file_structure_is_byte_identical(field_file):
    """Header, dimensions, patch names, BC types, counts and delimiters.

    Compared exactly (after newline normalisation) -- this is what catches a
    wrong face count, a dropped patch, or wrong dimensions.
    """
    name, expected, written = field_file
    exp_struct, _ = _split_body(expected)
    got_struct, _ = _split_body(written)
    assert got_struct == exp_struct, f"{name}: structural lines differ"


def test_field_file_values_match_numerically(field_file):
    """Numeric payload, parsed back out and compared at rtol=1e-12.

    Not a byte comparison: boundaryT emits `300` from a Python int where a config
    carrying T0_K: 300.0 yields `300.0`. See README for both formatting quirks.
    """
    name, expected, written = field_file
    _, exp_body = _split_body(expected)
    _, got_body = _split_body(written)

    assert len(got_body) == len(exp_body), f"{name}: value count differs"

    exp_vals = np.array([[float(v) for v in NUMBER.findall(l)] for l in exp_body])
    got_vals = np.array([[float(v) for v in NUMBER.findall(l)] for l in got_body])
    np.testing.assert_allclose(got_vals, exp_vals, rtol=RTOL, atol=0)


def test_field_files_use_lf_line_endings(written_dir):
    """plumetools writes newline='\\n' explicitly so output is platform-independent.

    The legacy code wrote via print() to a text-mode file, so it emitted CRLF on
    Windows and LF on Linux for identical inputs.
    """
    for name in ("boundaryU", "boundaryT", "boundaryNumberDensity_Ar"):
        raw = (written_dir / name).read_bytes()
        assert b"\r\n" not in raw, f"{name}: expected LF-only output"


# --------------------------------------------------------------------------- #
# 3. geometry invariants of the committed mesh
# --------------------------------------------------------------------------- #

def test_inflow_faces_are_triangles(computed):
    assert {len(f) for f in computed.faces} == {3}


def test_inflow_vertices_lie_on_the_sphere(computed):
    """The committed mesh's inflow surface is a sphere of radius exactly 0.5 m."""
    r = np.linalg.norm(computed.vertices, axis=1)
    np.testing.assert_allclose(r, 0.5, rtol=0, atol=1e-9)


def test_inflow_centroids_sit_just_inside_the_sphere(computed):
    """Triangle centroids fall inside the circumscribed sphere -- pure faceting.

    Measured span is 0.4989640 to 0.4997035 m, i.e. 0.06-0.21% inside the exact
    0.5 m surface. This is why the hard-coded r = 0.5 is the correct nominal
    radius and per-face |c| would be worse: using |c| would inject this faceting
    spread into f3 = (r_e/r)^2 as if it were physics (finding SM-09).
    """
    r = np.linalg.norm(computed.centroids, axis=1)
    assert 0.4989 <= r.min()
    assert r.max() <= 0.4998
    assert (r < 0.5).all()


def test_inflow_is_a_hemisphere_opening_toward_plus_x(computed):
    assert (computed.centroids[:, 0] >= 0.0).all()


# --------------------------------------------------------------------------- #
# 4. physical sanity -- reported, NOT scientific acceptance criteria
# --------------------------------------------------------------------------- #

def test_velocity_magnitude_is_the_limiting_velocity(computed):
    """|U| = v_l on every face (E7). Frozen value for gamma=1.4, T0=300 K, N2."""
    mag = np.linalg.norm(computed.U, axis=1)
    np.testing.assert_allclose(mag, 788.164111, rtol=1e-8)


def test_velocity_is_radially_outward(computed):
    """U points outward along the centroid direction (E7) -- to within faceting.

    Not exactly parallel, and the reason is worth recording. E7 builds U from the
    spherical angles, and those are computed as acos(z / 0.5) using the *nominal*
    sphere radius applied to the *centroid's* z. Because the centroid sits at
    |c| ~ 0.4994 rather than 0.5, the recovered polar angle is not the centroid's
    own polar angle, so U tilts off the radial direction by the faceting offset.

    Measured worst case is 3.1e-4 (0.014 deg) -- negligible physically, but it is
    a real consequence of mixing a fixed radius with per-face coordinates, and it
    would grow on a coarser inflow surface. Part of finding SM-09.
    """
    c = computed.centroids
    c_hat = c / np.linalg.norm(c, axis=1, keepdims=True)
    u_hat = computed.U / np.linalg.norm(computed.U, axis=1, keepdims=True)
    alignment = (c_hat * u_hat).sum(axis=1)
    assert alignment.min() > 0.999, "velocity is not outward-radial"
    assert (1.0 - alignment).max() < 1e-3, "radial misalignment larger than faceting"


def test_number_density_is_finite_and_positive(computed):
    """Guards the SM-06 NaN path: f_phi goes negative-base for |phi| > theta_l.

    Latent on this mesh (0/2044 faces exceed it) but live on the archived
    wake-cylinder meshes (99/99). Any new inflow surface could trip it.
    """
    assert np.isfinite(computed.rhoN).all()
    assert (computed.rhoN > 0).all()


def test_number_density_peaks_near_the_plume_axis(computed):
    """Density should be largest where the surface is closest to +x."""
    peak = computed.centroids[np.argmax(computed.rhoN)]
    assert peak[0] / np.linalg.norm(peak) > 0.9
