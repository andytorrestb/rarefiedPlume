"""Small independent physical/mesh fixtures for portable refinement tests."""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import yaml

from plumetools.nozzleplume.gas import VhsSpecies


def unit_inputs(*, edges, lambdas, parcels):
    """Use hard spheres of diameter 1 so requested MFPs are hand-checkable."""
    from plumetools.refinement.model import Geometry, Measurements

    n = len(edges)
    gas = VhsSpecies(name="X", mass_kg=1e-26, diameter_m=1.0, omega=0.5)
    geometry = Geometry(
        centres_m=np.column_stack((np.arange(n) + 0.5, np.zeros(n), np.zeros(n))),
        volumes_m3=np.ones(n), max_edge_m=np.asarray(edges), mesh_sha256="unit-mesh",
    )
    measured = Measurements(
        number_density=1.0 / (math.sqrt(2.0) * math.pi * np.asarray(lambdas)),
        temperature_K=np.full(n, 300.0), mean_parcels=np.asarray(parcels),
        window={"start_s": 1.0, "end_s": 2.0, "basis": "mean_per_cell", "source": "declared"},
        source_sha256="unit-measurements", provenance={},
    )
    return geometry, measured, gas


def tree_hash(root: Path) -> dict[str, str]:
    """Hash source files, using names that survive moving the fixture."""
    result = {}
    for path in sorted(Path(root).rglob("*")):
        if path.is_file():
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(65536), b""):
                    digest.update(chunk)
            result[path.relative_to(root).as_posix()] = digest.hexdigest()
    return result


def write_field(path, values, dimensions, *, uniform=False):
    """Write a real ASCII volume field on the fixture's one wall patch."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    a = np.asarray(values)
    vector = a.ndim == 2
    kind = "vector" if vector else "scalar"
    def fmt(value):
        if vector:
            return "(" + " ".join(format(float(v), ".17g") for v in value) + ")"
        return format(float(value), ".17g")
    if uniform:
        payload = "uniform " + fmt(a[0]) + ";"
    else:
        payload = f"nonuniform List<{kind}>\n{len(a)}\n(\n" + "\n".join(fmt(v) for v in a) + "\n);"
    text = (f'FoamFile\n{{\n version 2.0;\n format ascii;\n class vol{kind.title()}Field;\n'
            f' location "{path.parent.name}";\n object {path.name};\n}}\n'
            f'dimensions [{" ".join(str(v) for v in dimensions)}];\ninternalField {payload}\n'
            'boundaryField\n{\n walls { type zeroGradient; }\n}\n'
            '// ************************************************************************* //\n')
    path.write_text(text, encoding="utf-8")


def write_cartesian_case(root, *, profile="custom", lambdas=(0.75, 2.0), parcels=(40, 40)):
    """Two adjacent unit cubes, independently specified faces and owner labels."""
    from _helpers import HEADER, make_case, patch

    root = Path(root)
    points = [(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1, 2)]
    faces = [[1, 4, 10, 7], [0, 6, 9, 3], [2, 5, 11, 8],
             [0, 3, 4, 1], [1, 4, 5, 2], [6, 7, 10, 9], [7, 8, 11, 10],
             [0, 1, 7, 6], [1, 2, 8, 7], [3, 9, 10, 4], [4, 10, 11, 5]]
    make_case(root, points=points, faces=faces, patches=[patch("walls", 10, 1, "wall")])
    for name, labels in (("owner", [0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]), ("neighbour", [1])):
        (root / "constant/polyMesh" / name).write_text(
            HEADER.format(cls="labelList", obj=name) + f'\n{len(labels)}\n(\n'
            + "\n".join(str(v) for v in labels) + "\n)\n", encoding="utf-8")
    dimless, n_dim, t_dim = [0]*7, [0, -3, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0]
    write_field(root / "1/C", [[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]], [0, 1, 0, 0, 0, 0, 0])
    write_field(root / "1/V", [1, 1], [0, 3, 0, 0, 0, 0, 0])
    density = 1.0 / (math.sqrt(2) * math.pi * np.asarray(lambdas))
    fields = {"profile": profile}
    if profile == "custom":
        fields.update(number_density="nObserved", temperature="TObserved", mean_parcels="countObserved")
        names = ("nObserved", "TObserved", "countObserved")
    elif profile == "hystrath":
        names = ("rhoN_Xe", "Ttra_Xe", "dsmcNMean_Xe")
    else:
        names = ("rhoNMean", "unusedT", "dsmcRhoNMean")
    write_field(root / "1" / names[0], density, n_dim)
    write_field(root / "1" / names[2], parcels, dimless)
    if profile != "standard":
        write_field(root / "1" / names[1], [300, 300], t_dim)
    else:
        gas = VhsSpecies(name="Xe", mass_kg=1e-26, diameter_m=1.0)
        rho = density * gas.mass_kg
        write_field(root / "1/rhoMMean", rho, [1, -3, 0, 0, 0, 0, 0])
        write_field(root / "1/momentumMean", np.column_stack((rho*10, rho*0, rho*0)), [1, -2, -1, 0, 0, 0, 0])
        write_field(root / "1/linearKEMean", 0.5*rho*(100 + 3*gas.gas_constant_j_per_kg_k*300), [1, -1, -2, 0, 0, 0, 0])
        state = root / "1/uniform/functionObjects/functionObjectProperties"
        state.parent.mkdir(parents=True, exist_ok=True)
        state.write_text("\n".join(f"{name}\n{{ totalIter 10; totalTime 1; }}" for name in
                                  ("rhoN", "rhoM", "momentum", "linearKE", "dsmcRhoN")), encoding="utf-8")
    settings = {
        "gas": {"species_name": "Xe", "molecular_mass_kg": 1e-26,
                "molecular_diameter_m": 1.0, "omega": 0.5, "t_ref_K": 273.0},
        "fields": fields,
        "sampling": {} if profile == "standard" else {
            "start_s": 0, "end_s": 1, "basis": "mean_per_cell", "source": "declared"},
    }
    (root / "refinement.yaml").write_text(yaml.safe_dump(settings), encoding="utf-8")
    return root
