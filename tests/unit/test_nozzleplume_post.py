"""Mesh verification and post-processing, against surfaces whose answers are known.

Both read only the boundary of a mesh, so both are tested here on a synthetic
one built from the plan (see ``_nozzleplume.build_boundary_mesh``): the nozzle
faces the plan predicts and the six sides of the plate, with nothing behind
them. Every area, normal and load on it is known in closed form.

The same code is run against what OpenFOAM really builds in
``tests/openfoam/test_nozzleplume_openfoam.py``, where a collisionless jet on a
diffuse plate has to return the force and the heat that momentum and energy
conservation say it must.
"""

from __future__ import annotations

import numpy as np
import pytest
import yaml

from _nozzleplume import build_boundary_mesh, derived, write_field
from plumetools.nozzleplume import meshread, post, verify
from plumetools.nozzleplume.verify import MeshVerificationError

TIME = "0.004"

#: The minimal case's plate, upstream face first: a 4 x 2 m face in a 4 x 2 grid
#: of 1 m^2 quads, then the lee face and the four edges, one quad each.
FRONT, OTHERS = 8, 5


def case(tmp_path, **sections):
    """A planned minimal case with its synthetic boundary mesh on disk."""
    cfg, geom, exit_state, plan, nozzle, run = derived(tmp_path, **sections)
    build_boundary_mesh(tmp_path, geom, plan)
    return cfg, geom, exit_state, plan, nozzle, run


def plate_fields(tmp_path, force_density, heat_flux, plate="plate"):
    """Write fDMean and qMean with the given per-face plate values."""
    write_field(tmp_path / TIME / "fDMean", field_class="volVectorField",
                boundary={plate: force_density, "vacuum": None})
    write_field(tmp_path / TIME / "qMean", boundary={plate: heat_flux,
                                                     "vacuum": None})


# --------------------------------------------------------------------------- #
# reading patch geometry
# --------------------------------------------------------------------------- #

def test_face_geometry_of_a_quad_and_a_triangle():
    points = np.array([(0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0)], dtype=float)
    centres, areas = meshread.face_geometry([[0, 1, 2, 3], [0, 1, 2]], points)
    assert areas[0] == pytest.approx([0.0, 0.0, 2.0])          # 2 x 1, facing +z
    assert areas[1] == pytest.approx([0.0, 0.0, 1.0])
    assert centres[0] == pytest.approx([1.0, 0.5, 0.0])


def test_reversing_a_face_reverses_its_normal():
    points = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], dtype=float)
    _, areas = meshread.face_geometry([[0, 1, 2, 3], [3, 2, 1, 0]], points)
    assert areas[0] == pytest.approx(-areas[1])


def test_the_plate_patch_is_read_with_its_normals_into_the_plate(tmp_path):
    case(tmp_path)
    plate = meshread.read_patches(tmp_path, ["plate"])["plate"]
    assert plate.n_faces == FRONT + OTHERS
    assert plate.area_m2 == pytest.approx(2 * 4.0 * 2.0 + 2 * 0.25 * (4.0 + 2.0))
    # Upstream face: the fluid is at smaller x, so "out of the fluid" is +x.
    assert plate.normals[:FRONT] == pytest.approx(np.tile([1.0, 0.0, 0.0], (FRONT, 1)))
    assert plate.normals[FRONT] == pytest.approx([-1.0, 0.0, 0.0])
    # A closed surface: its area vectors sum to nothing.
    assert plate.area_vectors.sum(axis=0) == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)


def test_an_absent_patch_is_named_with_the_ones_present(tmp_path):
    case(tmp_path)
    with pytest.raises(KeyError, match="panel.*nozzle"):
        meshread.read_patches(tmp_path, ["panel"])


def test_the_cell_count_comes_from_the_owner_note(tmp_path):
    cfg, geom, _, plan, _, _ = derived(tmp_path)
    build_boundary_mesh(tmp_path, geom, plan, n_cells=529126)
    assert meshread.read_cell_count(tmp_path) == 529126
    assert meshread.read_cell_count(tmp_path / "nowhere") is None


# --------------------------------------------------------------------------- #
# verifying a mesh against case.yaml
# --------------------------------------------------------------------------- #

def test_the_planned_mesh_verifies(tmp_path):
    cfg, geom, _, plan, nozzle, _ = case(tmp_path)
    report = verify.verify(tmp_path, cfg, geom, plan, nozzle)
    assert not [line for line in report if "FAIL" in line]
    assert any("32 faces the plan predicted" in line for line in report)
    assert any("upstream face is 8 m^2" in line for line in report)


def test_a_symmetric_mesh_verifies_its_quarter(tmp_path):
    cfg, geom, _, plan, nozzle, _ = case(
        tmp_path, domain={"symmetry_planes": ["y", "z"]})
    report = verify.verify(tmp_path, cfg, geom, plan, nozzle)
    assert any("8 faces the plan predicted" in line for line in report)
    assert any("upstream face is 2 m^2" in line for line in report)
    assert any("symmetry plane y = 0 is flat" in line for line in report)


def test_a_free_plume_verifies_without_a_plate(tmp_path):
    cfg, geom, _, plan, nozzle, _ = case(tmp_path, plate=None)
    report = verify.verify(tmp_path, cfg, geom, plan, nozzle)
    assert not any("plate" in line for line in report)


def failure(tmp_path, **mesh_changes) -> str:
    cfg, geom, _, plan, nozzle, _ = derived(tmp_path)
    build_boundary_mesh(tmp_path, geom, plan, **mesh_changes)
    with pytest.raises(MeshVerificationError) as excinfo:
        verify.verify(tmp_path, cfg, geom, plan, nozzle)
    return str(excinfo.value).splitlines()[-1]


def test_a_nozzle_missing_faces_is_not_the_planned_mesh(tmp_path):
    assert "32 faces the plan predicted (meshed 30)" in failure(
        tmp_path, drop_nozzle_faces=2)


def test_a_plate_at_the_wrong_standoff_fails(tmp_path):
    line = failure(tmp_path, plate_bounds=(4.1, 4.35, -2.0, 2.0, -1.0, 1.0))
    assert "upstream face x is 4 m (measured 4.1)" in line


def test_a_plate_of_the_wrong_width_fails(tmp_path):
    assert "y max is 2 m (measured 1.9)" in failure(
        tmp_path, plate_bounds=(4.0, 4.25, -2.0, 1.9, -1.0, 1.0))


def test_a_plate_meshed_as_a_patch_fails(tmp_path):
    """It would absorb the plume."""
    assert "'plate' is geometric type 'wall' (got 'patch')" in failure(
        tmp_path, patch_types={"plate": "patch"})


def test_a_vacuum_boundary_meshed_as_a_wall_fails(tmp_path):
    assert "'vacuum' is geometric type 'patch' (got 'wall')" in failure(
        tmp_path, patch_types={"vacuum": "wall"})


def test_a_patch_case_yaml_does_not_describe_fails(tmp_path):
    """A stale mesh from an earlier configuration."""
    assert "unexpected: ['cylinder']" in failure(tmp_path, extra_patches=["cylinder"])


def test_a_mesh_over_the_cell_ceiling_fails(tmp_path):
    assert "within mesh.max_cells" in failure(tmp_path, n_cells=5_000_000)


def test_a_failure_carries_the_checks_that_passed_before_it(tmp_path):
    cfg, geom, _, plan, nozzle, _ = derived(tmp_path)
    build_boundary_mesh(tmp_path, geom, plan, drop_nozzle_faces=2)
    with pytest.raises(MeshVerificationError) as excinfo:
        verify.verify(tmp_path, cfg, geom, plan, nozzle)
    text = str(excinfo.value)
    assert text.count("  OK   ") >= 5 and text.count("  FAIL ") == 1


# --------------------------------------------------------------------------- #
# reading fields
# --------------------------------------------------------------------------- #

def test_internal_fields_are_read_in_every_form_openfoam_writes(tmp_path):
    scalar = write_field(tmp_path / "a", internal=[1.0, 2.5, -3.0])
    vector = write_field(tmp_path / "b", internal=[(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)],
                         field_class="volVectorField")
    uniform = write_field(tmp_path / "c")
    assert post.read_internal_field(scalar).tolist() == [1.0, 2.5, -3.0]
    assert post.read_internal_field(vector).tolist() == [[1.0, 2.0, 3.0],
                                                         [4.0, 5.0, 6.0]]
    assert post.read_internal_field(uniform).tolist() == [0.0]


def test_a_short_list_written_on_one_line_is_read(tmp_path):
    """OpenFOAM writes a list of a few entries inline."""
    path = tmp_path / "inline"
    path.write_text("internalField   nonuniform List<vector> 2((1 2 3) (4 5 6));\n"
                    "boundaryField\n{\n}\n", encoding="utf-8")
    assert post.read_internal_field(path).tolist() == [[1.0, 2.0, 3.0],
                                                       [4.0, 5.0, 6.0]]


def test_the_internal_field_is_not_confused_with_a_patch_of_the_same_shape(tmp_path):
    path = write_field(tmp_path / "f", internal=[1.0, 2.0],
                       boundary={"plate": [7.0, 8.0, 9.0]})
    assert post.read_internal_field(path).tolist() == [1.0, 2.0]
    assert post.read_patch_field(path, "plate").tolist() == [7.0, 8.0, 9.0]


def test_a_uniform_patch_value_comes_back_as_one_row(tmp_path):
    path = write_field(tmp_path / "f", boundary={"plate": (1.0, 2.0, 3.0)},
                       field_class="volVectorField")
    assert post.read_patch_field(path, "plate").tolist() == [[1.0, 2.0, 3.0]]


def test_a_patch_with_no_value_is_data_that_was_never_written(tmp_path):
    """That is what zeroGradient looks like on disk -- not a zero load."""
    path = write_field(tmp_path / "f", boundary={"plate": None})
    with pytest.raises(post.PostError, match="no 'value' entry"):
        post.read_patch_field(path, "plate")


def test_a_patch_the_field_does_not_have_is_reported(tmp_path):
    path = write_field(tmp_path / "f", boundary={"plate": [1.0]})
    with pytest.raises(post.PostError, match="'panel' is not in"):
        post.read_patch_field(path, "panel")


def test_the_latest_time_is_numeric_not_alphabetical(tmp_path):
    for name in ("0", "0.0005", "0.004", "9e-05", "constant", "processor0"):
        (tmp_path / name).mkdir()
    assert post.latest_time(tmp_path) == "0.004"


def test_a_case_with_no_results_says_so(tmp_path):
    (tmp_path / "0").mkdir()
    with pytest.raises(post.PostError, match="solver wrote nothing"):
        post.latest_time(tmp_path)


# --------------------------------------------------------------------------- #
# plate loads
# --------------------------------------------------------------------------- #

def test_a_uniform_pressure_gives_pressure_times_area(tmp_path):
    """100 Pa on the 8 m^2 upstream face, nothing elsewhere: 800 N along +x."""
    cfg, geom, exit_state, *_ = case(tmp_path)
    plate_fields(tmp_path,
                 [(100.0, 0.0, 0.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS,
                 [50.0] * FRONT + [0.0] * OTHERS)
    loads = post.plate_loads(tmp_path, cfg, geom)

    assert loads.time == TIME
    assert loads.force_n == pytest.approx([800.0, 0.0, 0.0])
    assert loads.pressure[:FRONT] == pytest.approx(100.0)
    assert loads.shear == pytest.approx(0.0, abs=1e-12)
    assert loads.heat_load_w == pytest.approx(50.0 * 8.0)
    assert loads.side.tolist() == [1] * FRONT + [-1] + [0] * 4


def test_pressure_is_the_normal_part_and_shear_the_rest(tmp_path):
    cfg, geom, *_ = case(tmp_path)
    plate_fields(tmp_path,
                 [(100.0, 30.0, 40.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS,
                 [0.0] * (FRONT + OTHERS))
    loads = post.plate_loads(tmp_path, cfg, geom)
    assert loads.pressure[:FRONT] == pytest.approx(100.0)
    assert loads.shear[:FRONT] == pytest.approx(50.0)
    assert loads.force_n == pytest.approx([800.0, 240.0, 320.0])


def test_pushing_on_the_lee_face_is_a_positive_pressure_there_too(tmp_path):
    """fD . n with n INTO the plate: gas behind the plate pushes it along -x."""
    cfg, geom, *_ = case(tmp_path)
    values = [(0.0, 0.0, 0.0)] * (FRONT + OTHERS)
    values[FRONT] = (-10.0, 0.0, 0.0)
    plate_fields(tmp_path, values, [0.0] * (FRONT + OTHERS))
    loads = post.plate_loads(tmp_path, cfg, geom)
    assert loads.pressure[FRONT] == pytest.approx(10.0)
    assert loads.force_n[0] == pytest.approx(-10.0 * 8.0)


def test_a_quarter_model_reports_the_whole_plate(tmp_path):
    """The mesh holds a 2 x 1 m quarter. The axial force is four times what it
    measures; the lateral force cancels against its mirror images."""
    cfg, geom, exit_state, *_ = case(tmp_path,
                                     domain={"symmetry_planes": ["y", "z"]})
    plate_fields(tmp_path,
                 [(100.0, 5.0, 7.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS,
                 [50.0] * FRONT + [0.0] * OTHERS)
    loads = post.plate_loads(tmp_path, cfg, geom)

    assert loads.multiplicity == 4
    assert loads.modelled_force_n == pytest.approx([200.0, 10.0, 14.0])
    assert loads.force_n == pytest.approx([800.0, 0.0, 0.0])
    assert loads.heat_load_w == pytest.approx(4 * 50.0 * 2.0)
    assert loads.summary(exit_state)["upstream_face"]["area_m2"] == pytest.approx(8.0)


def test_a_half_model_keeps_the_force_along_its_symmetry_plane(tmp_path):
    cfg, geom, *_ = case(tmp_path, domain={"symmetry_planes": ["y"]})
    plate_fields(tmp_path,
                 [(100.0, 5.0, 7.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS,
                 [0.0] * (FRONT + OTHERS))
    force = post.plate_loads(tmp_path, cfg, geom).force_n
    assert force == pytest.approx([800.0, 0.0, 2 * 7.0 * 4.0])


def test_the_summary_sets_the_loads_beside_the_jet(tmp_path):
    cfg, geom, exit_state, *_ = case(tmp_path)
    thrust, power = exit_state.thrust_n(), exit_state.power_w()
    plate_fields(tmp_path,
                 [(thrust / 8.0, 0.0, 0.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS,
                 [0.5 * power / 8.0] * FRONT + [0.0] * OTHERS)
    summary = post.plate_loads(tmp_path, cfg, geom).summary(exit_state)
    assert summary["axial_force_over_jet_thrust"] == pytest.approx(1.0)
    assert summary["heat_load_over_jet_power"] == pytest.approx(0.5)
    assert summary["upstream_face"]["mean_pressure_pa"] == pytest.approx(thrust / 8.0)


def test_the_peak_pressure_is_located(tmp_path):
    cfg, geom, exit_state, *_ = case(tmp_path)
    values = [(10.0, 0.0, 0.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS
    values[5] = (90.0, 0.0, 0.0)
    plate_fields(tmp_path, values, [0.0] * (FRONT + OTHERS))
    loads = post.plate_loads(tmp_path, cfg, geom)
    front = loads.summary(exit_state)["upstream_face"]
    assert front["peak_pressure_pa"] == pytest.approx(90.0)
    assert front["peak_pressure_at_m"] == pytest.approx(loads.centres[5].tolist())


def test_a_field_from_another_mesh_is_refused(tmp_path):
    """Pairing loads with the wrong faces would still draw a plausible map."""
    cfg, geom, *_ = case(tmp_path)
    plate_fields(tmp_path, [(1.0, 0.0, 0.0)] * 5, [0.0] * 5)
    with pytest.raises(post.PostError, match="not the same mesh"):
        post.plate_loads(tmp_path, cfg, geom)


def test_a_run_that_never_reached_the_averaging_window_says_so(tmp_path):
    cfg, geom, *_ = case(tmp_path)
    (tmp_path / TIME).mkdir()
    with pytest.raises(post.PostError, match="fieldAverage"):
        post.plate_loads(tmp_path, cfg, geom)


def test_a_free_plume_has_no_plate_loads(tmp_path):
    cfg, geom, *_ = case(tmp_path, plate=None)
    with pytest.raises(post.PostError, match="no plate"):
        post.plate_loads(tmp_path, cfg, geom)


def test_plate_results_are_written(tmp_path):
    cfg, geom, exit_state, *_ = case(tmp_path)
    plate_fields(tmp_path,
                 [(100.0, 0.0, 0.0)] * FRONT + [(0.0, 0.0, 0.0)] * OTHERS,
                 [50.0] * FRONT + [0.0] * OTHERS)
    loads = post.plate_loads(tmp_path, cfg, geom)
    table, summary = post.write_plate_results(tmp_path / "results", loads, exit_state)

    rows = table.read_text(encoding="utf-8").splitlines()
    assert rows[0].split(",")[:6] == ["x_m", "y_m", "z_m", "area_m2", "side",
                                      "pressure_pa"]
    assert len(rows) == 1 + FRONT + OTHERS
    document = yaml.safe_load(summary.read_text(encoding="utf-8"))
    assert document["force_on_plate_n"] == pytest.approx([800.0, 0.0, 0.0])


# --------------------------------------------------------------------------- #
# the centreline
# --------------------------------------------------------------------------- #

def centreline_case(tmp_path, cells):
    """A minimal case with ``cells`` as its whole sampled field.

    Each cell is ``(x, y, z, volume, n / n_exit, u / u_exit, T / T_exit)``.
    """
    cfg, geom, exit_state, *_ = derived(tmp_path)
    mass = exit_state.species.mass_kg
    r = exit_state.species.gas_constant_j_per_kg_k
    centres, volumes, n, rho, momentum, kinetic = [], [], [], [], [], []
    for x, y, z, volume, density, speed, temperature in cells:
        number = density * exit_state.number_density_per_m3
        u = speed * exit_state.speed_m_per_s
        centres.append((x, y, z))
        volumes.append(volume)
        n.append(number)
        rho.append(number * mass)
        momentum.append((number * mass * u, 0.0, 0.0))
        kinetic.append(number * mass * (0.5 * u * u
                                        + 1.5 * r * temperature * exit_state.temperature_K))
    directory = tmp_path / TIME
    write_field(directory / "C", internal=centres, field_class="volVectorField")
    write_field(directory / "V", internal=volumes)
    write_field(directory / "rhoNMean", internal=n)
    write_field(directory / "rhoMMean", internal=rho)
    write_field(directory / "momentumMean", internal=momentum,
                field_class="volVectorField")
    write_field(directory / "linearKEMean", internal=kinetic)
    return cfg, geom, exit_state


def test_the_centreline_recovers_density_velocity_and_temperature(tmp_path):
    """The moments dsmcFoam writes, inverted: u = momentum / rhoM and
    3 R T = <v^2> - u^2."""
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.5, 0.02, 0.02, 1e-4, 1.0, 1.0, 1.0),
        (1.5, 0.02, 0.02, 1e-4, 0.25, 0.9, 0.5),
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["n_over_n_exit"] == pytest.approx([1.0, 0.25])
    assert profile["U_over_U_exit"] == pytest.approx([1.0, 0.9])
    assert profile["T_over_T_exit"] == pytest.approx([1.0, 0.5])
    assert profile["x_m"] == pytest.approx([0.5, 1.5], abs=0.02)


def test_the_centreline_ignores_cells_outside_the_tube(tmp_path):
    """0.25 D of a 0.4 m exit is 0.1 m; a small cell 0.3 m out is not on axis."""
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.5, 0.02, 0.02, 1e-4, 1.0, 1.0, 1.0),
        (0.5, 0.30, 0.00, 1e-4, 9.0, 1.0, 1.0),
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["n_over_n_exit"] == pytest.approx([1.0])


def test_the_centreline_follows_the_axis_where_the_cells_outgrow_the_tube(tmp_path):
    """Downstream the octree's cells are bigger than the nozzle. A 0.25 m cell
    beside the axis has its centre 0.177 m out -- past a 0.1 m tube, which would
    sample nothing there -- but inside three quarters of its own size."""
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (3.0, 0.125, 0.125, 0.25 ** 3, 0.01, 1.0, 0.1),      # beside the axis
        (3.0, 0.375, 0.125, 0.25 ** 3, 9.0, 1.0, 0.1),       # the next ring out
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["n_over_n_exit"] == pytest.approx([0.01])


def test_the_centreline_stops_at_the_plate(tmp_path):
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.5, 0.02, 0.02, 1e-4, 1.0, 1.0, 1.0),
        (5.0, 0.02, 0.02, 1e-4, 1.0, 1.0, 1.0),      # behind the plate at 4 m
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["x_m"].size == 1


def test_the_centreline_density_is_volume_weighted(tmp_path):
    """Two cells in one bin: total molecules over total volume, not the mean of
    two densities -- which would let a sliver outvote a cell."""
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.500, 0.02, 0.02, 3e-4, 1.0, 1.0, 1.0),
        (0.505, 0.02, -0.02, 1e-4, 0.2, 1.0, 1.0),
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["n_over_n_exit"] == pytest.approx([(3 * 1.0 + 1 * 0.2) / 4])
    assert profile["n_cells"].tolist() == [2.0]


def test_a_shear_across_the_tube_is_not_read_as_temperature(tmp_path):
    """Two cold streams at different speeds in one bin. Thermal energy is taken
    cell by cell, so the spread between them stays out of T."""
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.500, 0.02, 0.02, 1e-4, 1.0, 1.0, 0.2),
        (0.505, 0.02, -0.02, 1e-4, 1.0, 0.5, 0.2),
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["T_over_T_exit"] == pytest.approx([0.2])
    assert profile["U_over_U_exit"] == pytest.approx([0.75])


def test_without_cell_centres_the_centreline_names_the_command(tmp_path):
    cfg, geom, exit_state, *_ = derived(tmp_path)
    (tmp_path / TIME).mkdir()
    with pytest.raises(post.PostError, match="writeCellCentres"):
        post.centreline(tmp_path, cfg, geom, exit_state)


def test_centreline_results_are_written(tmp_path):
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.5, 0.02, 0.02, 1e-4, 1.0, 1.0, 1.0)])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    path = post.write_centreline_results(tmp_path / "results", profile)
    header, row = path.read_text(encoding="utf-8").splitlines()
    assert header == ("x_m,x_over_D,n_over_n_exit,U_over_U_exit,T_over_T_exit,"
                      "n_cells,cell_size_m")
    assert float(row.split(",")[2]) == pytest.approx(1.0)


def test_each_sample_carries_the_size_of_the_cells_it_came_from(tmp_path):
    cfg, geom, exit_state = centreline_case(tmp_path, [
        (0.5, 0.02, 0.02, 0.0625 ** 3, 1.0, 1.0, 1.0),
        (3.0, 0.125, 0.125, 0.25 ** 3, 0.01, 1.0, 0.1),
    ])
    profile = post.centreline(tmp_path, cfg, geom, exit_state)
    assert profile["cell_size_m"] == pytest.approx([0.0625, 0.25])


def test_a_coarsening_of_the_sampling_cells_is_marked(tmp_path):
    """The plume is peaked on its axis, so cells twice the size average over
    more of the peak and read lower. The profile steps there; the flow does not.
    The figure has to say which steps are of that kind."""
    profile = {
        "x_m": np.array([0.5, 1.5, 2.5, 3.5, 4.5, 8.5]),
        # fine, fine, one bin straddling the interface, coarse, coarse, coarser
        "cell_size_m": np.array([0.0625, 0.0625, 0.09, 0.125, 0.125, 0.25]),
    }
    assert post.sampling_changes(profile) == [3.5, 8.5]


def test_a_refinement_towards_the_plate_is_marked_too(tmp_path):
    profile = {"x_m": np.array([17.0, 18.0, 19.0, 19.9]),
               "cell_size_m": np.array([0.25, 0.25, 0.125, 0.125])}
    assert post.sampling_changes(profile) == [19.0]


def test_a_uniform_mesh_has_nothing_to_mark():
    profile = {"x_m": np.linspace(0.0, 1.0, 5), "cell_size_m": np.full(5, 0.01)}
    assert post.sampling_changes(profile) == []
    assert post.sampling_changes({"x_m": np.array([])}) == []
