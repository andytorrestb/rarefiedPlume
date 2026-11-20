"""The uniGasFoam -> dsmcFoam field conversion (spec §5).

The converter must only re-express data: through the frozen
``plumetools.cai2012.post.read_case`` a converted case returns uniGasFoam's own
density, velocity and translational temperature, exactly.
"""

from __future__ import annotations

import hashlib
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from plumetools.cai2012.post import read_case
from plumetools.cai2012ugf.convert import ConversionError, convert_case
from plumetools.foamio.volfield import read_internal_field

REPO = Path(__file__).resolve().parents[2]
TEMPLATE = REPO / "cases" / "cai2012-kn-extension" / "baseCase" / "case.yaml"
MASS = 6.63e-26
K_B = 1.380649e-23

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {name};
}}
"""


def write_field(path: Path, name: str, values: np.ndarray, *, vector=False):
    values = np.asarray(values, dtype=float)
    kind = "vector" if vector else "scalar"
    rows = [f"({' '.join(repr(float(c)) for c in v)})" if vector else repr(float(v))
            for v in values]
    body = (HEADER.format(cls=f"vol{kind.capitalize()}Field", name=name)
            + "dimensions [0 0 0 0 0 0 0];\n"
            + f"internalField nonuniform List<{kind}>\n{len(rows)}\n(\n"
            + "\n".join(rows) + "\n)\n;\n"
            + "boundaryField\n{\n    nozzle { type calculated; value uniform 0; }\n"
              "    upstreamVacuum { type calculated; value uniform 0; }\n"
              "    vacuum { type calculated; value uniform 0; }\n}\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


@dataclass
class Synthetic:
    dir: Path
    rhoN: np.ndarray
    rhoM: np.ndarray
    U: np.ndarray
    T: np.ndarray
    startup_time: str = "0.03"
    sampled_time: str = "0.06"


@pytest.fixture
def synthetic_case(tmp_path):
    case = tmp_path / "case"
    case.mkdir()
    shutil.copy(TEMPLATE, case / "case.yaml")
    (case / "unigas.yaml").write_text("{}\n", encoding="utf-8")
    mesh = case / "constant" / "polyMesh"
    mesh.mkdir(parents=True)
    (mesh / "owner").write_text(
        'FoamFile\n{\n    version 2.0;\n    format ascii;\n    class labelList;\n'
        '    note "nPoints:8 nCells:5 nFaces:6 nInternalFaces:0";\n    object owner;\n}\n'
        "0\n()\n", encoding="utf-8")
    (mesh / "boundary").write_text(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n    class polyBoundaryMesh;\n"
        "    object boundary;\n}\n3\n(\n"
        "    nozzle\n    {\n        type patch;\n        nFaces 1;\n        startFace 0;\n    }\n"
        "    upstreamVacuum\n    {\n        type patch;\n        nFaces 1;\n        startFace 1;\n    }\n"
        "    vacuum\n    {\n        type patch;\n        nFaces 4;\n        startFace 2;\n    }\n"
        ")\n", encoding="utf-8")

    rhoN = np.array([6.6e21, 3.1e20, 4.5e18, 0.0, 1.2e17])
    rhoM = rhoN * MASS
    U = np.array([[700.0, 1.0, -2.0], [650.0, 30.0, 40.0], [600.0, -5.0, 3.0],
                  [0.0, 0.0, 0.0], [550.0, 100.0, -100.0]])
    T = np.array([300.0, 210.0, 90.5, 0.0, 30.25])
    for time in ("0.03", "0.06"):
        d = case / time
        write_field(d / "rhoN_Ar", "rhoN_Ar", rhoN)
        write_field(d / "rhoM_Ar", "rhoM_Ar", rhoM)
        write_field(d / "UMean_Ar", "UMean_Ar", U, vector=True)
        write_field(d / "translationalT_Ar", "translationalT_Ar", T)
        write_field(d / "uniGasRhoNMean_Ar", "uniGasRhoNMean_Ar",
                    np.where(rhoN > 0, np.array([40.0, 12.0, 3.0, 0.0, 1.5]), 0.001))
        write_field(d / "rhoN_ArWindow", "rhoN_ArWindow", rhoN * 1.01)
    sampled = case / "0.06"
    write_field(sampled / "C", "C", np.arange(15.0).reshape(5, 3) * 0.01, vector=True)
    write_field(sampled / "V", "V", np.full(5, 1e-6))
    return Synthetic(dir=case, rhoN=rhoN, rhoM=rhoM, U=U, T=T)


def digests(case: Path, pattern: str) -> dict:
    return {str(p.relative_to(case)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(case.glob(f"*/{pattern}"))}


def test_round_trip_is_exact(synthetic_case):
    convert_case(synthetic_case.dir)
    s = read_case(synthetic_case.dir, MASS, time=synthetic_case.sampled_time)
    occ = synthetic_case.rhoM > 0
    np.testing.assert_allclose(s.number_density, synthetic_case.rhoN, rtol=1e-12)
    np.testing.assert_allclose(s.velocity[occ], synthetic_case.U[occ], rtol=1e-12)
    np.testing.assert_allclose(s.temperature[occ], synthetic_case.T[occ], rtol=1e-12)


def test_empty_cells_stay_empty(synthetic_case):        # Review Focus 2
    convert_case(synthetic_case.dir)
    t = synthetic_case.dir / synthetic_case.sampled_time
    empty = synthetic_case.rhoM == 0
    assert (read_internal_field(t / "dsmcRhoNMean")[empty] == 0).all()   # source has 0.001 there
    assert (read_internal_field(t / "rhoMMean")[empty] == 0).all()
    assert (read_internal_field(t / "linearKEMean")[empty] == 0).all()


def test_no_mean_fields_before_average_start(synthetic_case):   # Review Focus 3
    convert_case(synthetic_case.dir)
    early = synthetic_case.dir / synthetic_case.startup_time
    assert (early / "rhoN").is_file() and not (early / "rhoNMean").exists()


def test_second_run_is_noop(synthetic_case):            # Review Focus 3
    convert_case(synthetic_case.dir)
    results = convert_case(synthetic_case.dir)
    assert all(not r.written for r in results)


def test_refuses_differing_target(synthetic_case):      # Review Focus 3
    convert_case(synthetic_case.dir)
    (synthetic_case.dir / synthetic_case.sampled_time / "rhoNMean").write_text("tampered")
    with pytest.raises(ConversionError, match="rhoNMean"):
        convert_case(synthetic_case.dir)


def test_never_touches_source_fields(synthetic_case):
    before = digests(synthetic_case.dir, "*_Ar")
    convert_case(synthetic_case.dir)
    assert digests(synthetic_case.dir, "*_Ar") == before


def test_header_names_converter_and_sources(synthetic_case):
    convert_case(synthetic_case.dir)
    text = (synthetic_case.dir / synthetic_case.sampled_time / "momentumMean").read_text()
    assert "plumetools.cai2012ugf.convert" in text and "rhoM_Ar" in text and "UMean_Ar" in text


def test_missing_window_field_is_reported_not_fatal(synthetic_case):
    (synthetic_case.dir / synthetic_case.sampled_time / "rhoN_ArWindow").unlink()
    r = {x.time: x for x in convert_case(synthetic_case.dir)}[synthetic_case.sampled_time]
    assert "rhoN" in r.skipped


def test_linear_ke_is_the_kinetic_energy_density(synthetic_case):
    convert_case(synthetic_case.dir)
    ke = read_internal_field(synthetic_case.dir / synthetic_case.sampled_time / "linearKEMean")
    expected = (0.5 * synthetic_case.rhoM * (synthetic_case.U ** 2).sum(axis=1)
                + 1.5 * synthetic_case.rhoN * K_B * synthetic_case.T)
    np.testing.assert_allclose(ke, expected, rtol=1e-15)
