"""Reading the solver's output, extracting the four outputs, and scoring them.

The DSMC fields are synthesised from the analytical solution, so the extraction
can be checked against a known answer: if :func:`plumetools.cai2012.post.centerline`
is fed a field that *is* Cai's solution, it must report zero error.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import yaml

from _cai2012 import derived
from plumetools.cai2012 import analytical, post
from plumetools.cai2012.constants import BOLTZMANN_J_PER_K

MASS = 6.63e-26


# --------------------------------------------------------------------------- #
# a synthetic case
# --------------------------------------------------------------------------- #

def write_field(path, name, values, field_class):
    """Write a volume field the way OpenFOAM does, boundaryField and all."""
    values = np.asarray(values, dtype=np.float64)
    if field_class == "volVectorField":
        rows = [f"({v[0]:.10g} {v[1]:.10g} {v[2]:.10g})" for v in values]
        kind = "vector"
    else:
        rows = [f"{v:.10g}" for v in values]
        kind = "scalar"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
        f"    class {field_class};\n    object {name};\n}}\n\n"
        "dimensions      [0 0 0 0 0 0 0];\n\n"
        f"internalField   nonuniform List<{kind}>\n{len(rows)}\n(\n"
        + "\n".join(rows)
        + "\n)\n;\n\n"
        "boundaryField\n{\n    nozzle\n    {\n        type calculated;\n"
        "        value uniform 0;\n    }\n}\n",
        encoding="utf-8", newline="\n")
    return path


def make_case(directory, geom, exit_state, time="0.05", *, nx=40, nl=9):
    """A case directory holding the ANALYTICAL solution as if dsmcFoam wrote it.

    The **full three-dimensional** solution, not the centreline values smeared
    laterally: the extraction averages over a tube of finite radius, so a field
    that were uniform across the tube would let a broken tube correction pass.
    """
    d = geom.diameter_m
    x = np.linspace(0.0, 10.0 * d, nx) + 0.5 * (10.0 * d / nx)
    lateral = np.linspace(-0.5 * d, 0.5 * d, nl)
    xx, yy, zz = np.meshgrid(x, lateral, lateral, indexing="ij")
    centres = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    field = analytical.moments(centres[:, 0], centres[:, 1], centres[:, 2],
                               geom.radius_m, exit_state.speed_ratio,
                               n_radial=32, n_azimuthal=64)
    n = field["n_over_n0"] * exit_state.number_density_per_m3
    velocity = field["U_sqrt_beta0"] / math.sqrt(exit_state.beta0)
    temperature = field["T_over_T0"] * exit_state.T0_K
    u = velocity[:, 0]

    rho_m = n * MASS
    momentum = rho_m[:, None] * velocity
    gas_constant = BOLTZMANN_J_PER_K / MASS
    mean_square = (3.0 * gas_constant * temperature
                   + (velocity ** 2).sum(axis=1))
    kinetic = 0.5 * rho_m * mean_square

    out = directory / time
    write_field(out / "C", "C", centres, "volVectorField")
    write_field(out / "V", "V", np.full(len(centres), 1.0e-6), "volScalarField")
    write_field(out / "rhoNMean", "rhoNMean", n, "volScalarField")
    write_field(out / "rhoMMean", "rhoMMean", rho_m, "volScalarField")
    write_field(out / "momentumMean", "momentumMean", momentum, "volVectorField")
    write_field(out / "linearKEMean", "linearKEMean", kinetic, "volScalarField")
    return directory


@pytest.fixture
def case(tmp_path):
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    make_case(tmp_path / "case", geom, exit_state)
    return cfg, geom, exit_state, tmp_path / "case"


# --------------------------------------------------------------------------- #
# reading fields
# --------------------------------------------------------------------------- #

def test_nonuniform_scalar_list_is_read(tmp_path):
    path = write_field(tmp_path / "f", "f", [1.0, 2.0, 3.0], "volScalarField")
    assert post.read_internal_field(path) == pytest.approx([1.0, 2.0, 3.0])


def test_nonuniform_vector_list_is_read(tmp_path):
    values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    path = write_field(tmp_path / "f", "f", values, "volVectorField")
    read = post.read_internal_field(path)
    assert read.shape == (2, 3)
    assert np.allclose(read, values)


def test_the_boundary_field_is_not_mistaken_for_the_internal_one(tmp_path):
    """The legacy line-offset parsers in this repository went wrong here."""
    path = write_field(tmp_path / "f", "f", [1.0, 2.0], "volScalarField")
    assert post.read_internal_field(path).size == 2


def test_a_uniform_internal_field_is_read(tmp_path):
    path = tmp_path / "f"
    path.write_text("internalField   uniform 7;\nboundaryField\n{\n}\n",
                    encoding="utf-8")
    assert post.read_internal_field(path) == pytest.approx([7.0])


def test_a_uniform_vector_internal_field_is_read(tmp_path):
    path = tmp_path / "f"
    path.write_text("internalField   uniform (1 2 3);\nboundaryField\n{\n}\n",
                    encoding="utf-8")
    assert np.allclose(post.read_internal_field(path), [[1.0, 2.0, 3.0]])


def test_a_missing_file_is_an_error(tmp_path):
    with pytest.raises(post.PostError, match="no field file"):
        post.read_internal_field(tmp_path / "absent")


def test_a_file_without_an_internal_field_is_an_error(tmp_path):
    path = tmp_path / "f"
    path.write_text("boundaryField\n{\n}\n", encoding="utf-8")
    with pytest.raises(post.PostError, match="no internalField"):
        post.read_internal_field(path)


def test_latest_time_ignores_zero(tmp_path):
    for name in ("0", "0.01", "0.05", "constant", "system"):
        (tmp_path / name).mkdir()
    assert post.latest_time(tmp_path) == "0.05"


def test_a_case_with_no_results_is_an_error(tmp_path):
    (tmp_path / "0").mkdir()
    with pytest.raises(post.PostError, match="wrote nothing"):
        post.latest_time(tmp_path)


# --------------------------------------------------------------------------- #
# deriving U and T from the moments
# --------------------------------------------------------------------------- #

def test_read_case_recovers_the_velocity_and_temperature(case):
    """u = momentum/rhoM and 3RT = 2*linearKE/rhoM - |u|^2, inverted exactly."""
    _, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    expected = analytical.moments(
        sampled.centres[:, 0], sampled.centres[:, 1], sampled.centres[:, 2],
        geom.radius_m, exit_state.speed_ratio, n_radial=32, n_azimuthal=64)
    # 1e-7, not machine precision: the synthetic field file carries ten
    # significant figures, as OpenFOAM's writePrecision does.
    assert np.allclose(
        sampled.velocity,
        expected["U_sqrt_beta0"] / math.sqrt(exit_state.beta0), rtol=1e-6,
        atol=1e-6)
    assert np.allclose(sampled.temperature,
                       expected["T_over_T0"] * exit_state.T0_K, rtol=1e-6)


def test_empty_cells_get_zero_not_a_division_by_zero(tmp_path):
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    # Blank one cell out, as a vacuum cell would be.
    values = post.read_internal_field(directory / "0.05" / "rhoMMean")
    values[0] = 0.0
    write_field(directory / "0.05" / "rhoMMean", "rhoMMean", values,
                "volScalarField")
    sampled = post.read_case(directory, MASS)
    assert sampled.temperature[0] == 0.0
    assert np.all(np.isfinite(sampled.temperature))


def test_missing_cell_centres_names_the_command_that_writes_them(tmp_path):
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    (directory / "0.05" / "C").unlink()
    with pytest.raises(post.PostError, match="writeCellCentres"):
        post.read_case(directory, MASS)


def test_a_missing_averaged_field_names_field_average(tmp_path):
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    (directory / "0.05" / "rhoNMean").unlink()
    with pytest.raises(post.PostError, match="fieldAverage"):
        post.read_case(directory, MASS)


def test_missing_cell_volumes_fall_back_to_unweighted(tmp_path):
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    (directory / "0.05" / "V").unlink()
    sampled = post.read_case(directory, MASS)
    assert np.all(sampled.volumes == 1.0)


# --------------------------------------------------------------------------- #
# geometry from another time
#
# The mesh does not move, so C and V are identical in every time directory.
# A study that post-processes every written frame -- which is how an
# error-against-averaging-time curve is built -- would otherwise duplicate
# 60 MB of geometry per frame.
# --------------------------------------------------------------------------- #

def test_the_geometry_can_come_from_another_time(tmp_path):
    """Fields from one directory, C and V from another, same answer.

    If mesh_time were ignored, this would raise on the missing C rather than
    read it from where it was pointed -- and a study relying on it would write
    the geometry into every frame after all.
    """
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    reference = post.read_case(directory, MASS, time="0.05")

    (directory / "0.10").mkdir()
    for name in ("rhoNMean", "rhoMMean", "momentumMean", "linearKEMean"):
        (directory / "0.05" / name).rename(directory / "0.10" / name)

    sampled = post.read_case(directory, MASS, time="0.10", mesh_time="0.05")
    assert sampled.time == "0.10"
    assert np.allclose(sampled.centres, reference.centres)
    assert np.allclose(sampled.volumes, reference.volumes)
    assert np.allclose(sampled.number_density, reference.number_density)


def test_geometry_from_a_different_mesh_fails_on_the_cell_count(tmp_path):
    """Pairing one mesh's densities with another's centres would plot cleanly.

    Every sample would sit at the wrong place, and nothing about the resulting
    centreline would look unusual -- so the cell count is checked rather than
    trusted.
    """
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    coarse = make_case(tmp_path / "coarse", geom, exit_state, time="0.05", nx=8)
    (coarse / "0.05" / "C").replace(directory / "0.05" / "C")

    with pytest.raises(post.PostError, match="not the same mesh"):
        post.read_case(directory, MASS)


def test_an_absent_mesh_time_is_named_rather_than_reported_as_missing_centres(
        tmp_path):
    """A typo'd --mesh-time should say the directory is absent, not that
    postProcess -func writeCellCentres was never run: the fix is different."""
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case")
    directory = make_case(tmp_path / "case", geom, exit_state)
    with pytest.raises(post.PostError, match="no time directory"):
        post.read_case(directory, MASS, mesh_time="0.99")


# --------------------------------------------------------------------------- #
# A / B / C -- the centreline
# --------------------------------------------------------------------------- #

def test_the_centreline_of_the_analytical_field_has_no_error(case):
    """The extraction fed Cai's own solution must reproduce it.

    Not exactly zero: the values are binned in x, so a bin averages the solution
    over its width and the curvature shows up as a small bias.
    """
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    metrics = post.centerline_metrics(profile)
    assert metrics["density_max_rel_error"] < 0.05
    assert metrics["velocity_max_rel_error"] < 0.02
    assert metrics["temperature_max_rel_error"] < 0.05


def test_the_centreline_is_normalised_by_the_diameter(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    assert profile["x_over_D"].min() >= 0.0
    assert profile["x_over_D"].max() <= cfg.post.centerline_x_over_D_max


def test_the_centreline_density_is_normalised_by_n0(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    assert profile["n_over_n0"].max() <= 1.01
    assert profile["n_over_n0"].min() > 0.0


def test_the_centreline_velocity_is_a_speed_ratio(case):
    """U1 sqrt(beta0) is dimensionless and starts near S0 = 2."""
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    assert profile["U_sqrt_beta0"][0] == pytest.approx(2.0, rel=0.05)
    assert profile["U_sqrt_beta0"][-1] > profile["U_sqrt_beta0"][0]


def test_empty_bins_are_dropped_not_reported_as_zero(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    assert profile["x_over_D"].size <= cfg.post.centerline_points
    assert np.all(profile["n_cells"] > 0)


def test_the_tube_average_sits_below_the_on_axis_value(case):
    """The plume is peaked on the axis, so averaging over a tube of finite
    radius biases the density DOWN. That bias is what the ``*_analytical``
    columns exist to remove; ``*_analytical_axis`` is Cai's own curve."""
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    near = (profile["x_over_D"] > 0.5) & (profile["x_over_D"] < 3.0)
    assert np.all(profile["n_over_n0_analytical"][near]
                  <= profile["n_over_n0_analytical_axis"][near])


def test_a_wide_tube_still_scores_near_zero_against_the_matched_analytical(tmp_path):
    """The correction is what makes the tube radius a statistics choice rather
    than a physics one: at 0.25 D the raw bias reaches 6.8% at X/D = 1, and the
    scored error must not.
    """
    cfg, geom, exit_state, _, _ = derived(
        tmp_path / "case", post={"centerline_radius_over_D": 0.25})
    directory = make_case(tmp_path / "case", geom, exit_state)
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)

    near = (profile["x_over_D"] > 0.5) & (profile["x_over_D"] < 2.0)
    raw = np.abs(profile["n_over_n0"][near]
                 / profile["n_over_n0_analytical_axis"][near] - 1.0)
    corrected = np.abs(profile["n_over_n0"][near]
                       / profile["n_over_n0_analytical"][near] - 1.0)
    assert raw.max() > 0.02
    assert corrected.max() < 0.01


def test_every_centreline_column_has_the_same_length(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    assert len({v.size for v in profile.values()}) == 1


# --------------------------------------------------------------------------- #
# D -- the density plane
# --------------------------------------------------------------------------- #

def test_the_density_plane_takes_the_layers_straddling_the_plane(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    plane = post.density_plane(sampled, cfg, geom, exit_state)
    assert plane["plane"] == "xz"
    assert plane["lateral_axis"] == "Z"
    assert plane["n_over_n0"].size > 0
    assert plane["n_over_n0"].size < sampled.n_cells


def test_the_density_plane_is_normalised(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    plane = post.density_plane(sampled, cfg, geom, exit_state)
    assert plane["n_over_n0"].max() <= 1.01
    assert np.all(np.abs(plane["lateral_over_D"])
                  <= cfg.post.plane_half_over_D + 1e-9)


def test_the_xy_plane_is_selectable(tmp_path):
    cfg, geom, exit_state, _, _ = derived(tmp_path / "case", post={"plane": "xy"})
    directory = make_case(tmp_path / "case", geom, exit_state)
    sampled = post.read_case(directory, MASS)
    plane = post.density_plane(sampled, cfg, geom, exit_state)
    assert plane["lateral_axis"] == "Y"


def test_the_analytical_plane_is_a_regular_grid(case):
    """The lateral count is rounded up to an odd number so the axis is a grid
    line and the two halves mirror exactly."""
    cfg, geom, exit_state, _ = case
    plane = post.analytical_plane(cfg, geom, exit_state)
    n_half = (cfg.post.plane_points_lateral + 1) // 2
    assert plane["n_over_n0"].shape == (cfg.post.plane_points_x, 2 * n_half - 1)
    assert np.all(plane["n_over_n0"] > 0.0)


def test_the_analytical_plane_is_symmetric_about_the_axis(case):
    """Mirrored, not integrated twice -- so this is exact, not approximate."""
    cfg, geom, exit_state, _ = case
    plane = post.analytical_plane(cfg, geom, exit_state)
    assert np.array_equal(plane["n_over_n0"], plane["n_over_n0"][:, ::-1])
    assert np.allclose(plane["lateral_over_D"], -plane["lateral_over_D"][:, ::-1])


def test_the_analytical_plane_starts_off_the_singular_exit_plane(case):
    """The disk integral is singular in the exit plane itself."""
    cfg, geom, exit_state, _ = case
    plane = post.analytical_plane(cfg, geom, exit_state)
    assert plane["x_over_D"].min() > 0.0


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #

def test_metrics_of_an_exact_match_are_zero():
    exact = np.array([1.0, 0.5, 0.25])
    metrics = post.error_metrics(exact, exact, name="density")
    assert metrics["density_max_rel_error"] == 0.0
    assert metrics["density_rms_rel_error"] == 0.0
    assert metrics["density_n_points"] == 3


def test_metrics_are_relative_not_absolute():
    metrics = post.error_metrics([1.1, 0.55], [1.0, 0.5], name="density")
    assert metrics["density_max_rel_error"] == pytest.approx(0.1)
    assert metrics["density_mean_rel_error"] == pytest.approx(0.1)


def test_the_worst_point_is_located():
    metrics = post.error_metrics([1.0, 2.0, 1.0], [1.0, 1.0, 1.0], name="density")
    assert metrics["density_max_rel_error_at_index"] == 1
    assert metrics["density_max_rel_error"] == pytest.approx(1.0)


def test_rms_is_at_least_the_mean():
    metrics = post.error_metrics([1.0, 3.0], [1.0, 1.0], name="d")
    assert metrics["d_rms_rel_error"] >= metrics["d_mean_rel_error"]


def test_zero_valued_reference_points_are_dropped():
    metrics = post.error_metrics([1.0, 1.0], [0.0, 1.0], name="d")
    assert metrics["d_n_points"] == 1
    assert metrics["d_max_rel_error"] == 0.0


def test_mismatched_lengths_are_an_error():
    with pytest.raises(ValueError, match="measured values"):
        post.error_metrics([1.0], [1.0, 2.0], name="d")


def test_centreline_metrics_cover_all_three_quantities(case):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    metrics = post.centerline_metrics(post.centerline(sampled, cfg, geom, exit_state))
    for name in ("density", "velocity", "temperature"):
        for statistic in ("max", "mean", "rms"):
            assert f"{name}_{statistic}_rel_error" in metrics
    assert "density_max_rel_error_at_x_over_D" in metrics


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #

def test_results_are_written_as_csv_and_yaml(case, tmp_path):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    plane = post.density_plane(sampled, cfg, geom, exit_state)
    metrics = {"case": {"name": "Kn100"},
               "centerline": post.centerline_metrics(profile)}

    paths = post.write_results(tmp_path / "results", profile, plane, metrics)
    assert [p.name for p in paths] == [
        "centerline.csv", "density_plane.csv", "metrics.yaml"]

    header = paths[0].read_text(encoding="utf-8").splitlines()[0]
    assert header.split(",")[:3] == [
        "x_over_D", "n_over_n0", "n_over_n0_analytical"]
    assert "z_over_D" in paths[1].read_text(encoding="utf-8").splitlines()[0]

    document = yaml.safe_load(paths[2].read_text(encoding="utf-8"))
    assert document["case"]["name"] == "Kn100"


def test_csv_columns_of_different_lengths_are_an_error(tmp_path):
    with pytest.raises(ValueError, match="different lengths"):
        post.write_csv(tmp_path / "x.csv", {"a": [1, 2], "b": [1]})


def test_written_results_use_lf_line_endings(case, tmp_path):
    cfg, geom, exit_state, directory = case
    sampled = post.read_case(directory, MASS)
    profile = post.centerline(sampled, cfg, geom, exit_state)
    plane = post.density_plane(sampled, cfg, geom, exit_state)
    for path in post.write_results(tmp_path / "r", profile, plane, {}):
        assert b"\r\n" not in path.read_bytes()


# --------------------------------------------------------------------------- #
# the study table
# --------------------------------------------------------------------------- #

def test_study_row_converts_fractions_to_percent():
    row = post.study_row({
        "case": {"name": "Kn100", "Kn": 100.0, "cai_reported_max_percent": 0.07},
        "centerline": {"density_max_rel_error": 0.0123},
    })
    assert row["density_max_rel_error_percent"] == pytest.approx(1.23)
    assert row["cai_reported_max_percent"] == 0.07


def test_study_row_tolerates_a_missing_metric():
    row = post.study_row({"case": {"name": "x"}, "centerline": {}})
    assert row["density_max_rel_error_percent"] is None


def test_the_study_table_has_a_header_and_one_row_per_case(tmp_path):
    rows = [post.study_row({"case": {"name": n, "Kn": k},
                            "centerline": {"density_max_rel_error": 0.01}})
            for n, k in (("Kn100", 100.0), ("Kn0p1", 0.1))]
    path = post.write_study_table(tmp_path / "t.csv", rows)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split(",") == list(post.STUDY_COLUMNS)
    assert len(lines) == 3
