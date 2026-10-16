"""The rarefaction diagnostics, against answers known in closed form.

A few focused checks, not a census: the mean free path against the
hard-sphere formula and against the scalar gas module it restates; the
gradient against fields it must reproduce exactly; and the handling of cells
that cannot be trusted. The same pipeline on a real dsmcFoam run, and
ParaView reading it, is ``tests/openfoam/test_nozzleplume_openfoam.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import yaml

from _nozzleplume import write_case, write_field
from plumetools.nozzleplume import gas
from plumetools.rarefaction import core, foam, summary
from plumetools.rarefaction.settings import DiagnosticSettings, SettingsError

SPECIES = gas.VhsSpecies(name="exhaust", mass_kg=3.4e-26, diameter_m=3.5e-10)


def graph_of_boxes(boxes, mirror_y0=False):
    """A cell graph for axis-aligned boxes ``(lo, hi)``: face neighbours share a
    face of positive area. ``mirror_y0`` puts a symmetry plane at y = 0."""
    lo = np.array([b[0] for b in boxes], dtype=float)
    hi = np.array([b[1] for b in boxes], dtype=float)
    owner, neighbour = [], []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            for axis in range(3):
                others = [a for a in range(3) if a != axis]
                touch = (math.isclose(hi[i, axis], lo[j, axis])
                         or math.isclose(hi[j, axis], lo[i, axis]))
                overlap = all(min(hi[i, a], hi[j, a]) - max(lo[i, a], lo[j, a]) > 1e-12
                              for a in others)
                if touch and overlap:
                    owner.append(i)
                    neighbour.append(j)
    centres = 0.5 * (lo + hi)
    mirror = np.flatnonzero(np.isclose(lo[:, 1], 0.0)) if mirror_y0 else np.empty(0, int)
    offsets = np.zeros((mirror.size, 3))
    offsets[:, 1] = -2.0 * centres[mirror, 1]
    return core.CellGraph(centres=centres, owner=np.array(owner),
                          neighbour=np.array(neighbour), mirror_cells=mirror,
                          mirror_offsets=offsets)


def octree_boxes():
    """One coarse layer (1 m cells) against a fine one (0.5 m): a 2:1 face."""
    boxes = [((0, y, z), (1, y + 1, z + 1)) for y in (0, 1) for z in (0, 1)]
    boxes += [((x, y, z), (x + 0.5, y + 0.5, z + 0.5))
              for x in (1.0, 1.5) for y in np.arange(0, 2, 0.5)
              for z in np.arange(0, 2, 0.5)]
    return boxes


# --------------------------------------------------------------------------- #
# the mean free path
# --------------------------------------------------------------------------- #

def test_hard_sphere_mean_free_path_is_the_formula_and_the_gas_modules():
    n = np.array([1e14, 3e18, 1.306e22])
    t = np.array([50.0, 498.04, 9000.0])
    lam = core.mean_free_path(n, t, SPECIES)
    assert lam == pytest.approx(1.0 / (math.sqrt(2.0) * math.pi * 3.5e-10 ** 2 * n),
                                rel=1e-14)
    assert lam == pytest.approx([gas.mean_free_path(a, SPECIES, b)
                                 for a, b in zip(n, t)], rel=1e-14)
    # The R-4D exit, as its README quotes it.
    assert lam[2] == pytest.approx(0.141e-3, rel=0.01)


def test_vhs_mean_free_path_is_the_gas_modules_equation_g1():
    argon = gas.VhsSpecies(name="Ar", mass_kg=6.63e-26, diameter_m=4.17e-10,
                           omega=0.81, t_ref_K=273.0)
    n, t = np.array([1e19, 1e20]), np.array([100.0, 2000.0])
    for convention in ("vhs", "hard_sphere"):
        assert core.mean_free_path(n, t, argon, convention=convention) == \
            pytest.approx([gas.mean_free_path(a, argon, b, convention=convention)
                           for a, b in zip(n, t)], rel=1e-14)
    assert np.isnan(core.mean_free_path([0.0, -1.0], [300.0, 300.0], argon)).all()
    assert np.isnan(core.mean_free_path([1e19], [-5.0], argon)).all()


def test_translational_temperature_is_recovered_from_the_moments():
    r = SPECIES.gas_constant_j_per_kg_k
    n, u, t = 1e20, np.array([3000.0, -40.0, 5.0]), 750.0
    rho = n * SPECIES.mass_kg
    kinetic = rho * (0.5 * u @ u + 1.5 * r * t)
    assert core.translational_temperature([rho], [rho * u], [kinetic], r) == \
        pytest.approx([t], rel=1e-10)
    # One parcel, however often sampled, has no thermal spread: its moments
    # differ only by round-off, and that is not a temperature. Nor is rho = 0.
    single = rho * 0.5 * (u @ u) * (1.0 + 1e-12)
    assert np.isnan(core.translational_temperature(
        [rho, 0.0], [rho * u, (0, 0, 0)], [single, 0.0], r)).all()


# --------------------------------------------------------------------------- #
# the gradient
# --------------------------------------------------------------------------- #

def test_the_gradient_is_exact_for_a_linear_field_across_a_level_transition():
    graph = graph_of_boxes(octree_boxes())
    g = np.array([2.0, -3.0, 0.5])
    phi = 7.0 + graph.centres @ g
    grad = core.least_squares_gradient(phi, np.ones(phi.size, bool), graph)
    assert grad.valid.all()
    assert grad.vectors == pytest.approx(np.tile(g, (phi.size, 1)), abs=1e-12)

    # Garbage in an INVALID neighbour does not reach its neighbours' fit.
    valid = np.ones(phi.size, bool)
    valid[5] = False
    phi[5] = 1e30
    grad = core.least_squares_gradient(phi, valid, graph)
    assert not grad.valid[5]
    fitted = grad.valid
    assert grad.vectors[fitted] == pytest.approx(np.tile(g, (fitted.sum(), 1)), abs=1e-9)
    # A cell left with only coplanar valid neighbours gets no gradient at all
    # rather than a guess -- and only a neighbour of the lost cell can be one.
    lost = set(np.flatnonzero(valid & ~fitted))
    beside = set(graph.neighbour[graph.owner == 5]) | set(graph.owner[graph.neighbour == 5])
    assert lost <= beside and len(lost) < len(beside)


def test_a_cell_with_no_valid_neighbours_has_no_gradient():
    graph = graph_of_boxes(octree_boxes())
    valid = np.zeros(graph.n_cells, bool)
    valid[0] = True
    grad = core.least_squares_gradient(np.ones(graph.n_cells), valid, graph)
    assert not grad.valid.any() and np.isnan(grad.vectors).all()


def test_a_symmetry_plane_holds_the_normal_gradient_of_a_symmetric_field():
    """phi = y^2 is symmetric about y = 0. With the cell's mirror image the fit
    beside the plane is the central difference -- exact for a quadratic. One-sided
    it would be twice the true slope."""
    boxes = [((0, y, 0), (1, y + 0.5, 1)) for y in (0.0, 0.5, 1.0, 1.5)]
    boxes += [((1, y, 0), (2, y + 0.5, 1)) for y in (0.0, 0.5, 1.0, 1.5)]
    boxes += [((x, y, 1), (x + 1, y + 0.5, 2)) for x in (0, 1) for y in (0.0, 0.5, 1.0, 1.5)]
    for mirror, expected in ((True, 0.5), (False, 1.0)):
        graph = graph_of_boxes(boxes, mirror_y0=mirror)
        phi = graph.centres[:, 1] ** 2
        grad = core.least_squares_gradient(phi, np.ones(phi.size, bool), graph)
        beside = np.flatnonzero(np.isclose(graph.centres[:, 1], 0.25))
        assert grad.vectors[beside, 1] == pytest.approx(expected, abs=1e-12)


def context(n, t, graph, *, parcels=1e3, steps=100, settings=None, **inputs):
    size = n.size
    data = {"number_density": n, "temperature": t,
            "parcels_per_step": np.full(size, float(parcels)), "averaged_steps": steps,
            "cell_volume": np.full(size, 0.001), "cell_extent": np.full((size, 3), 0.1),
            "cell_graph": graph, **inputs}
    return core.Context({k: v for k, v in data.items() if v is not None},
                        species=SPECIES, settings=settings or DiagnosticSettings(),
                        extra={"KnD": core.knudsen_on("D", 0.4)})


def test_gradient_knudsen_numbers_of_known_fields():
    """n = n0 (1 + x/L): |grad n| / n = 1 / (L + x), so KnGradN = lambda / (L + x)
    cell by cell -- and a uniform field gives a valid zero."""
    boxes = [((x, y, z), (x + 0.1, y + 0.1, z + 0.1))
             for x in np.arange(0, 1.0, 0.1) for y in (0, 0.1) for z in (0, 0.1)]
    graph = graph_of_boxes(boxes)
    x = graph.centres[:, 0]
    n0, length = 1e20, 2.0
    n = n0 * (1.0 + x / length)
    t = 300.0 * (1.0 + x)                       # |grad T| / T = 1 / (1 + x)
    ctx = context(n, t, graph)
    outputs, skipped = core.evaluate(ctx, core.output_names(ctx))
    assert not skipped
    lam = core.mean_free_path(n, t, SPECIES)
    assert outputs["KnGradN"].values == pytest.approx(lam / (length + x), rel=1e-9)
    assert outputs["KnGradT"].values == pytest.approx(lam / (1.0 + x), rel=1e-9)
    assert outputs["KnD"].values == pytest.approx(lam / 0.4, rel=1e-14)
    assert outputs["cellSizeTarget"].values == pytest.approx(lam / 3.0, rel=1e-14)
    assert outputs["cellOverMfpMax"].values == pytest.approx(0.1 / lam, rel=1e-14)

    flat = context(np.full(x.size, n0), np.full(x.size, 300.0), graph)
    outputs, _ = core.evaluate(flat, ["KnGradN", "KnGradT"])
    for name in ("KnGradN", "KnGradT"):
        assert outputs[name].valid.all() and np.all(outputs[name].values == 0.0)


# --------------------------------------------------------------------------- #
# cells that cannot be trusted
# --------------------------------------------------------------------------- #

def test_invalid_cells_are_flagged_kept_out_of_the_numbers_and_summaries():
    """Empty, at the threshold (it takes MORE than 20), negative T, and fine --
    counting accumulated samples; then the default, parcels per step."""
    graph = graph_of_boxes([((i, 0, 0), (i + 1, 1, 1)) for i in range(4)])
    n = np.array([0.0, 1e20, 1e20, 1e20])
    t = np.array([np.nan, 300.0, -4.0, 300.0])
    per_step = np.array([0.0, 0.2, 5.0, 5.0])                    # x 100 steps
    accumulated = DiagnosticSettings(particle_count_basis="accumulated")
    ctx = context(n, t, graph, parcels=0, parcels_per_step=per_step,
                  settings=accumulated)
    status = ctx["status"]
    assert status.tolist() == [core.NO_PARCELS | core.BAD_DENSITY | core.BAD_TEMPERATURE,
                               core.FEW_PARCELS, core.BAD_TEMPERATURE, 0]
    outputs, _ = core.evaluate(ctx, core.output_names(ctx))
    lam = outputs["mfp"]
    assert lam.valid.tolist() == [False, False, False, True]
    # Nothing to compute from in the empty cell or the one with T < 0; the
    # under-sampled one keeps its own value, to be shown but not trusted.
    assert np.isnan(lam.values[[0, 2]]).all()
    assert lam.values[1] == pytest.approx(core.mean_free_path(1e20, 300.0, SPECIES))
    assert outputs["diagValid"].values.tolist() == [0, 0, 0, 1]
    # 0.1 m cells, lambda = 18 mm: not resolved where there is a value; no data
    # where there is none.
    assert outputs["mfpResolved"].values.tolist() == [-1, 0, -1, 0]
    assert outputs["sampledParcels"].values.tolist() == [0, 20, 500, 500]
    # Nothing is ever substituted: not an exit value, not a neighbour's.
    stats = summary.statistics(outputs["cellOverMfpV"].values,
                               outputs["cellOverMfpV"].valid, np.ones(4))
    assert stats["cells"] == 1

    # The default counts parcels per step, where 5 is not more than 20.
    ctx = context(n, t, graph, parcels=0, parcels_per_step=per_step)
    assert ctx["status"][3] == core.FEW_PARCELS


def test_a_missing_input_skips_only_what_needs_it_and_says_why():
    graph_free = context(np.full(3, 1e20), np.full(3, 300.0), None)
    graph_free.unavailable["cell_graph"] = "no C at t = 0.1; postProcess writes it"
    outputs, skipped = core.evaluate(graph_free, core.output_names(graph_free))
    assert {"KnGradN", "KnGradT", "diagValidGradN", "diagValidGradT"} == set(skipped)
    assert "postProcess writes it" in skipped["KnGradN"]
    assert {"mfp", "KnD", "cellOverMfpMax", "diagValid"} <= set(outputs)


def test_settings_refuse_an_invalid_value_that_could_be_a_diagnostic():
    with pytest.raises(SettingsError, match="negative"):
        DiagnosticSettings.from_mapping({"invalid_value": 0.0})
    with pytest.raises(SettingsError, match="unknown key"):
        DiagnosticSettings.from_mapping({"min_particles": 20})
    assert DiagnosticSettings.from_mapping(None) == DiagnosticSettings()


# --------------------------------------------------------------------------- #
# planning
# --------------------------------------------------------------------------- #

def test_a_ten_millimetre_floor_permits_level_5_of_a_half_metre_background():
    levels = summary.refinement_levels(0.5, 0.010)
    assert levels["finest_permitted_level"] == 5
    assert levels["finest_permitted_cell_size_m"] == 0.015625
    assert levels["first_excluded_cell_size_m"] == 0.0078125
    assert [row["permitted"] for row in levels["levels"]] == [True] * 6 + [False]
    # The HPC variant's 1 m background reaches the same sizes one level later.
    assert summary.refinement_levels(1.0, 0.010)["finest_permitted_level"] == 6
    # The coarsest level whose cell is no larger than the target.
    assert summary.required_level(np.array([0.02, 0.5, 0.6, 4.7e-5]), 0.5).tolist() \
        == [5, 0, 0, 14]


# --------------------------------------------------------------------------- #
# the pipeline, end to end, on a synthetic case
# --------------------------------------------------------------------------- #

TIME = "0.01"


def write_box_mesh(case, shape, lo, size):
    """A structured hex mesh in OpenFOAM's on-disk format; one patch per side."""
    nx, ny, nz = shape
    poly = case / "constant" / "polyMesh"
    poly.mkdir(parents=True)

    def p(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    def c(i, j, k):
        return i + nx * (j + ny * k)

    points = [(lo[0] + i * size[0], lo[1] + j * size[1], lo[2] + k * size[2])
              for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    internal, sides = [], {name: [] for name in ("xmin", "xmax", "ymin", "ymax",
                                                 "zmin", "zmax")}
    for k in range(nz):
        for j in range(ny):
            for i in range(nx + 1):
                face = [p(i, j, k), p(i, j + 1, k), p(i, j + 1, k + 1), p(i, j, k + 1)]
                if 0 < i < nx:
                    internal.append((face, c(i - 1, j, k), c(i, j, k)))
                else:
                    sides["xmin" if i == 0 else "xmax"].append(
                        (face[::-1] if i == 0 else face, c(min(i, nx - 1), j, k)))
    for k in range(nz):
        for j in range(ny + 1):
            for i in range(nx):
                face = [p(i, j, k), p(i, j, k + 1), p(i + 1, j, k + 1), p(i + 1, j, k)]
                if 0 < j < ny:
                    internal.append((face, c(i, j - 1, k), c(i, j, k)))
                else:
                    sides["ymin" if j == 0 else "ymax"].append(
                        (face[::-1] if j == 0 else face, c(i, min(j, ny - 1), k)))
    for k in range(nz + 1):
        for j in range(ny):
            for i in range(nx):
                face = [p(i, j, k), p(i + 1, j, k), p(i + 1, j + 1, k), p(i, j + 1, k)]
                if 0 < k < nz:
                    internal.append((face, c(i, j, k - 1), c(i, j, k)))
                else:
                    sides["zmin" if k == 0 else "zmax"].append(
                        (face[::-1] if k == 0 else face, c(i, j, min(k, nz - 1))))

    faces = [f for f, _, _ in internal]
    owner = [o for _, o, _ in internal]
    neighbour = [n for _, _, n in internal]
    boundary = []
    for name, entries in sides.items():
        boundary.append((name, len(faces), len(entries)))
        faces += [f for f, _ in entries]
        owner += [o for _, o in entries]

    def listing(name, cls, rows):
        (poly / name).write_text(
            f"FoamFile\n{{\n    format ascii;\n    class {cls};\n    object {name};\n}}\n\n"
            f"{len(rows)}\n(\n" + "\n".join(rows) + "\n)\n", encoding="utf-8")

    listing("points", "vectorField", [f"({x:.10g} {y:.10g} {z:.10g})" for x, y, z in points])
    listing("faces", "faceList", [f"4({' '.join(map(str, f))})" for f in faces])
    listing("owner", "labelList", [str(o) for o in owner])
    listing("neighbour", "labelList", [str(n) for n in neighbour])
    (poly / "boundary").write_text(
        "FoamFile\n{\n    format ascii;\n    class polyBoundaryMesh;\n    object boundary;\n}\n\n"
        f"{len(boundary)}\n(\n" + "".join(
            f"    {name}\n    {{\n        type {'symmetry' if name == 'ymin' else 'patch'};\n"
            f"        nFaces {count};\n        startFace {start};\n    }}\n"
            for name, start, count in boundary) + ")\n", encoding="utf-8")
    cells = [((lo[0] + (i + 0.5) * size[0]), (lo[1] + (j + 0.5) * size[1]),
              (lo[2] + (k + 0.5) * size[2]))
             for k in range(nz) for j in range(ny) for i in range(nx)]
    return np.array(cells)


def synthetic_case(tmp_path, *, steps=(400, 400)):
    """The minimal nozzle-plume case.yaml over a box of 1 x 0.5 x 0.5 m cells,
    symmetric about y = 0, with averaged fields from a known state: n falls
    linearly along x, and no parcel ever reached x > 3 m."""
    write_case(tmp_path)
    centres = write_box_mesh(tmp_path, (6, 4, 4), (0.0, 0.0, -1.0), (1.0, 0.5, 0.5))
    x = centres[:, 0]
    weight, volume = 1e15, 0.25                     # 250 parcels per step
    n = np.where(x < 3.0, 1e18 * (1.0 - x / 8.0), 0.0)
    r = SPECIES.gas_constant_j_per_kg_k
    u, t = 2000.0, 400.0
    rho = n * SPECIES.mass_kg
    directory = tmp_path / TIME
    write_field(directory / "C", internal=[tuple(c) for c in centres],
                field_class="volVectorField")
    write_field(directory / "V", internal=[volume] * x.size)
    write_field(directory / "rhoNMean", internal=n.tolist())
    write_field(directory / "rhoMMean", internal=rho.tolist())
    write_field(directory / "momentumMean", internal=[(m * u, 0.0, 0.0) for m in rho],
                field_class="volVectorField")
    write_field(directory / "linearKEMean",
                internal=(rho * (0.5 * u * u + 1.5 * r * t)).tolist())
    write_field(directory / "dsmcRhoNMean", internal=(n * volume / weight).tolist())
    names = ("rhoN", "rhoM", "momentum", "linearKE", "dsmcRhoN")
    state = directory / "uniform" / "functionObjects"
    state.mkdir(parents=True)
    (state / "functionObjectProperties").write_text(
        "fieldAverage1\n{\n" + "".join(
            f"    {name}\n    {{\n        totalIter {steps[name == 'dsmcRhoN']};\n"
            f"        totalTime {steps[name == 'dsmcRhoN'] * 1e-6};\n    }}\n"
            for name in names) + "}\n", encoding="utf-8")
    return centres, n



def test_fields_averaged_over_different_windows_are_refused(tmp_path):
    synthetic_case(tmp_path, steps=(400, 300))
    with pytest.raises(foam.PostError, match="different windows"):
        foam.load(tmp_path, DiagnosticSettings())
