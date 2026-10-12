"""Helpers for the nozzle-plume unit tests: a minimal, valid case on disk.

Every test that needs a config goes through :func:`load_case`, which writes a
``case.yaml`` and loads it with the real loader. Constructing
:class:`~plumetools.nozzleplume.config.CaseConfig` directly would skip
``validate()``, and the validation is a substantial part of what these tests are
checking.

The template is deliberately the *smallest* file that loads, not a copy of
``cases/r-4d/baseCase/case.yaml``: a fixture that mirrored the shipped case
would pass even if the shipped case and the schema drifted apart together.
``tests/unit/test_nozzleplume_config.py`` loads the real one separately.
"""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

from plumetools.nozzleplume import exitstate, geometry, mesh, run
from plumetools.nozzleplume.config import load_case_config

REPO = Path(__file__).resolve().parents[2]

#: The shipped R-4D case.
R4D_CASE = REPO / "cases" / "r-4d" / "baseCase"

#: The smallest case.yaml this family's loader accepts, as a data structure.
#:
#: Round numbers, so an expected value can be worked out by hand: a 0.4 m exit,
#: a 4 x 2 m plate 4 m downstream, a 6 x 6 x 4 m box of 0.5 m cells, and one
#: refinement cone over the exit that widens from 0.4 m to 0.6 m in its metre.
MINIMAL_CASE = {
    "model": "nozzle_exit_plume",
    "meta": {"case_name": "test"},
    "nozzle": {"exit_diameter_m": 0.4},
    "exit": {
        "temperature_K": 500.0,
        "speed_m_per_s": 3000.0,
        "number_density_per_m3": 1.0e+22,
    },
    "gas": {
        "species_name": "exhaust",
        "molecular_mass_kg": 3.4e-26,
        "molecular_diameter_m": 3.5e-10,
    },
    "plate": {
        "standoff_m": 4.0,
        "thickness_m": 0.25,
        "size_y_m": 4.0,
        "size_z_m": 2.0,
    },
    "domain": {"x_max_m": 6.0, "y_half_m": 3.0, "z_half_m": 2.0},
    "mesh": {
        "background_cell_size_m": 0.5,
        "plate_refinement_level": 2,
        "refinement_regions": [
            {"name": "exitCore", "level": 3, "x_max_m": 1.0,
             "radius1_m": 0.4, "radius2_m": 0.6},
        ],
    },
    "dsmc": {"n_equivalent_particles": 1.0e+17},
}


def _merge(target: dict, updates: dict) -> dict:
    """Merge ``updates`` into a deep copy of ``target``, one section deep.

    A section given as ``None`` is removed, which is how a test asks for a case
    with no plate.
    """
    result = copy.deepcopy(target)
    for section, values in (updates or {}).items():
        if values is None:
            result.pop(section, None)
        elif isinstance(values, dict) and isinstance(result.get(section), dict):
            result[section] = {**result[section], **values}
        else:
            result[section] = values
    return result


def write_case(directory: Path, **sections) -> Path:
    """Write a ``case.yaml`` under ``directory``, with per-section overrides.

    ``write_case(tmp_path, exit={"temperature_K": 300.0})`` gives the minimal
    case with that one key changed; ``plate=None`` gives it with no plate.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "case.yaml"
    path.write_text(
        yaml.safe_dump(_merge(MINIMAL_CASE, sections), sort_keys=False),
        encoding="utf-8")
    return path


def load_case(directory: Path, **sections):
    """Write and load a case config in one step."""
    write_case(directory, **sections)
    return load_case_config(directory)


def derive(cfg):
    """``(geom, exit_state, plan, nozzle, run)`` for a loaded config."""
    geom = geometry.from_config(cfg)
    exit_state = exitstate.from_config(cfg)
    plan = mesh.plan(cfg, geom, exit_state)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)
    settings = run.derive(cfg, exit_state, geom, plan, nozzle)
    return geom, exit_state, plan, nozzle, settings


def derived(directory: Path, **sections):
    """``(cfg, geom, exit_state, plan, nozzle, run)`` -- the whole chain."""
    cfg = load_case(directory, **sections)
    return (cfg, *derive(cfg))


# --------------------------------------------------------------------------- #
# a synthetic boundary mesh
# --------------------------------------------------------------------------- #
#
# verify and post read only the BOUNDARY of a mesh: the patch list, the patches'
# faces and the points they use. So a "mesh" here is exactly that -- the nozzle
# faces the plan predicts, the six sides of the plate, one face for each open
# boundary -- with no cells behind it. It is built from the plan, not from
# OpenFOAM, which is the point: these tests then check the measuring code
# against a surface whose every dimension is known, and tests/openfoam checks
# the same code against what OpenFOAM really builds.

def _quads(axis: int, at: float, lo: tuple, hi: tuple, n: tuple, sign: int) -> list:
    """An ``n[0] x n[1]`` grid of quads on the plane ``axis = at``.

    Each quad is four ``(x, y, z)`` corners, wound so its normal points along
    ``sign * axis`` -- OpenFOAM's convention being that a boundary face's normal
    points out of the fluid.
    """
    u_axis, v_axis = [a for a in range(3) if a != axis]
    quads = []
    for i in range(n[0]):
        for j in range(n[1]):
            u0 = lo[0] + (hi[0] - lo[0]) * i / n[0]
            u1 = lo[0] + (hi[0] - lo[0]) * (i + 1) / n[0]
            v0 = lo[1] + (hi[1] - lo[1]) * j / n[1]
            v1 = lo[1] + (hi[1] - lo[1]) * (j + 1) / n[1]
            corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
            # (u, v, axis) is right-handed for axis 0 and 2 and left-handed for
            # axis 1, so the winding that faces +axis differs between them.
            if (sign > 0) != (axis != 1):
                corners.reverse()
            quad = []
            for u, v in corners:
                point = [0.0, 0.0, 0.0]
                point[axis], point[u_axis], point[v_axis] = at, u, v
                quad.append(tuple(point))
            quads.append(quad)
    return quads


def plate_quads(bounds: tuple, n_front: tuple = (4, 2)) -> list:
    """The six sides of a plate box, normals pointing INTO it.

    The upstream face is first and is an ``n_front`` grid, so a test can put a
    load distribution on it; the other five sides are one quad each.
    """
    x0, x1, y0, y1, z0, z1 = bounds
    return (_quads(0, x0, (y0, z0), (y1, z1), n_front, +1)
            + _quads(0, x1, (y0, z0), (y1, z1), (1, 1), -1)
            + _quads(1, y0, (x0, z0), (x1, z1), (1, 1), +1)
            + _quads(1, y1, (x0, z0), (x1, z1), (1, 1), -1)
            + _quads(2, z0, (x0, y0), (x1, y1), (1, 1), +1)
            + _quads(2, z1, (x0, y0), (x1, y1), (1, 1), -1))


def nozzle_quads(geom, plan) -> list:
    """The ``x = 0`` lattice faces whose centres are inside the exit disk."""
    import math

    h = plan.exit_cell_size_m
    quads = []
    reach = int(math.ceil(geom.radius_m / h)) + 1
    for j in range(-reach, reach):
        for k in range(-reach, reach):
            y, z = (j + 0.5) * h, (k + 0.5) * h
            if bool(geom.is_on_nozzle(y, z)):
                quads += _quads(0, 0.0, (j * h, k * h), ((j + 1) * h, (k + 1) * h),
                                (1, 1), -1)
    return quads


def build_boundary_mesh(root: Path, geom, plan, *, plate_bounds=None,
                        patch_types=None, drop_nozzle_faces: int = 0,
                        extra_patches=(), n_cells: int = 1000,
                        n_front: tuple = (4, 2)) -> dict:
    """Write a boundary-only ``constant/polyMesh`` for a planned case.

    Args:
        root: the case directory.
        geom, plan: the case's geometry and mesh plan.
        plate_bounds: the plate as meshed; the configured plate by default.
        patch_types: patch name -> geometric type, overriding the correct ones.
        drop_nozzle_faces: leave this many nozzle faces out.
        extra_patches: names of patches ``case.yaml`` does not describe.
        n_cells: the cell count the mesh claims in its ``owner`` note.
        n_front: the grid on the plate's upstream face.

    Returns:
        ``{patch name: number of faces}``.
    """
    from _helpers import HEADER, make_case, patch, write

    patch_types = dict(patch_types or {})
    groups = [
        ("upstreamVacuum", "patch",
         _quads(0, 0.0, (2.0, 1.0), (2.5, 1.5), (1, 1), -1)),
        ("vacuum", "patch",
         _quads(0, geom.x_max_m, (0.0, 0.0), (0.5, 0.5), (1, 1), +1)),
    ]
    if geom.symmetry_planes:
        faces = []
        for axis in geom.symmetry_planes:
            index = "xyz".index(axis)
            faces += _quads(index, 0.0, (0.0, 0.0), (0.5, 0.5), (1, 1), -1)
        groups.append(("symmetry", "symmetry", faces))
    if geom.plate is not None:
        groups.append(("plate", "wall", plate_quads(
            plate_bounds or geom.modelled_plate_bounds, n_front)))
    nozzle = nozzle_quads(geom, plan)
    groups.append(("nozzle", "patch", nozzle[:len(nozzle) - drop_nozzle_faces]))
    groups += [(name, "patch", _quads(0, geom.x_max_m, (1.0, 1.0), (1.5, 1.5),
                                      (1, 1), +1)) for name in extra_patches]

    points, faces, patches, counts = [], [], [], {}
    index: dict = {}
    for name, kind, quads in groups:
        patches.append(patch(name, len(quads), len(faces),
                             ptype=patch_types.get(name, kind)))
        counts[name] = len(quads)
        for quad in quads:
            labels = []
            for corner in quad:
                if corner not in index:
                    index[corner] = len(points)
                    points.append(corner)
                labels.append(index[corner])
            faces.append(labels)

    make_case(root, points=points, faces=faces, patches=patches)
    write(Path(root) / "constant" / "polyMesh" / "owner",
          HEADER.format(cls="labelList", obj="owner").replace(
              "    object ",
              f'    note        "nPoints:{len(points)}  nCells:{n_cells}  '
              f'nFaces:{len(faces)}  nInternalFaces:0";\n    object ')
          + f"\n{len(faces)}\n(\n" + "\n".join("0" for _ in faces) + "\n)\n")
    return counts


def write_field(path: Path, *, internal=None, boundary=None,
                field_class: str = "volScalarField") -> Path:
    """Write an OpenFOAM field file the way ``dsmcFoam`` does.

    Args:
        path: ``<case>/<time>/<name>``.
        internal: the cell values -- a list of scalars or of 3-tuples -- or
            ``None`` for ``uniform 0``.
        boundary: patch name -> a list of values (written ``nonuniform``), a
            single scalar or 3-tuple (written ``uniform``), or ``None`` for a
            ``zeroGradient`` entry with no value at all.
    """
    def one(value) -> str:
        if isinstance(value, (tuple, list)):
            return "(" + " ".join(f"{v:.10g}" for v in value) + ")"
        return f"{value:.10g}"

    def many(values) -> str:
        kind = "vector" if isinstance(values[0], (tuple, list)) else "scalar"
        rows = "\n".join(one(v) for v in values)
        return f"nonuniform List<{kind}> \n{len(values)}\n(\n{rows}\n)\n"

    lines = [
        "FoamFile", "{", "    version     2.0;", "    format      ascii;",
        f"    class       {field_class};", f"    object      {path.name};", "}", "",
        "dimensions      [0 0 0 0 0 0 0];", "",
    ]
    lines.append("internalField   uniform 0;" if internal is None
                 else f"internalField   {many(list(internal))};")
    lines += ["", "boundaryField", "{"]
    for name, values in (boundary or {}).items():
        lines += [f"    {name}", "    {"]
        if values is None:
            lines.append("        type            zeroGradient;")
        elif isinstance(values, list):
            lines += ["        type            calculated;",
                      f"        value           {many(values)};"]
        else:
            lines += ["        type            calculated;",
                      f"        value           uniform {one(values)};"]
        lines.append("    }")
    lines += ["}", ""]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
