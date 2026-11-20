"""Re-express uniGasFoam's output as the dsmcFoam fields plumetools reads.

About twenty modules -- several of them frozen -- read dsmcFoam's field names
(``rhoNMean``, ``momentumMean``, ``linearKEMean`` ...). uniGasFoam writes the
same physical moments under other names. This module writes the dsmcFoam names
from them, so the existing analysis runs unchanged (spec §5, Approach 1):

=================  ==========================================================
written            from uniGasFoam (``<f>`` mean field, ``<w>`` window field)
=================  ==========================================================
``rhoNMean``       ``rhoN_<f>``
``rhoMMean``       ``rhoM_<f>``
``momentumMean``   ``rhoM_<f> * UMean_<f>``
``linearKEMean``   ``0.5 rhoM |U|^2 + 1.5 rhoN k_B translationalT_<f>``
``dsmcRhoNMean``   ``uniGasRhoNMean_<f>``, 0 where ``rhoN_<f>`` is 0
``rhoN``           ``rhoN_<w>``
=================  ==========================================================

``plumetools.cai2012.post.read_case`` recovers U = momentum / rhoM and
T = (2 KE / rhoM - |U|^2) / 3R, so the round trip returns uniGasFoam's own
mass-weighted velocity and translational temperature exactly.

The ``*Mean`` fields are written only for writes after ``average_start_s``:
up to and including the write that ends the startup leg there, uniGasFoam's
fields are per-interval averages, and dsmcFoam had no ``*Mean`` either. ``rhoN`` is written wherever the window field exists.

Nothing is overwritten: a target that already exists is kept when it holds the
same values and refused otherwise.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from plumetools.foamio.errors import FoamDataError
from plumetools.foamio.primitives import BANNER, SEPARATOR, write_lines
from plumetools.foamio.volfield import read_internal

GENERATOR = "plumetools.cai2012ugf.convert"
K_BOLTZMANN = 1.380649e-23
#: Relative. The startup leg ends on a write at average_start_s, named to
#: timePrecision 8 and so rounded either way; that write holds the last startup
#: interval's average, not a sampled mean. Only writes clearly after it are sampled.
TIME_TOLERANCE = 1.0e-6

#: target -> (class, dimensions, needs the averaging window to have started)
MAPPING = {
    "rhoNMean": ("volScalarField", "[0 -3 0 0 0 0 0]", True),
    "rhoMMean": ("volScalarField", "[1 -3 0 0 0 0 0]", True),
    "momentumMean": ("volVectorField", "[1 -2 -1 0 0 0 0]", True),
    "linearKEMean": ("volScalarField", "[1 -1 -2 0 0 0 0]", True),
    "dsmcRhoNMean": ("volScalarField", "[0 -3 0 0 0 0 0]", True),
    "rhoN": ("volScalarField", "[0 -3 0 0 0 0 0]", False),
}


class ConversionError(RuntimeError):
    """A conversion that would change data already on disk."""


@dataclass
class ConversionResult:
    time: str
    written: list = field(default_factory=list)
    kept: list = field(default_factory=list)
    skipped: list = field(default_factory=list)


def _sources(fields) -> dict:
    f, w = fields.mean_field, fields.window_field
    return {
        "rhoNMean": (f"rhoN_{f}",),
        "rhoMMean": (f"rhoM_{f}",),
        "momentumMean": (f"rhoM_{f}", f"UMean_{f}"),
        "linearKEMean": (f"rhoN_{f}", f"rhoM_{f}", f"UMean_{f}", f"translationalT_{f}"),
        "dsmcRhoNMean": (f"uniGasRhoNMean_{f}", f"rhoN_{f}"),
        "rhoN": (f"rhoN_{w}",),
    }


def _values(time_dir: Path, name: str, n_cells: int) -> np.ndarray:
    return np.asarray(read_internal(time_dir / name).over(n_cells, f"cells in {time_dir}"),
                      dtype=np.float64)


def _compute(target: str, fields, read) -> np.ndarray:
    f = fields.mean_field
    if target == "rhoNMean":
        return read(f"rhoN_{f}")
    if target == "rhoMMean":
        return read(f"rhoM_{f}")
    if target == "momentumMean":
        return read(f"rhoM_{f}")[:, None] * read(f"UMean_{f}")
    if target == "linearKEMean":
        rho_m, u = read(f"rhoM_{f}"), read(f"UMean_{f}")
        return (0.5 * rho_m * (u ** 2).sum(axis=1)
                + 1.5 * read(f"rhoN_{f}") * K_BOLTZMANN * read(f"translationalT_{f}"))
    if target == "dsmcRhoNMean":
        return np.where(read(f"rhoN_{f}") > 0.0, read(f"uniGasRhoNMean_{f}"), 0.0)
    if target == "rhoN":
        return read(f"rhoN_{fields.window_field}")
    raise KeyError(target)


def _format(values: np.ndarray, vector: bool) -> list[str]:
    if vector:
        return [f"({v[0]:.17g} {v[1]:.17g} {v[2]:.17g})" for v in values]
    return [f"{v:.17g}" for v in values]


def _render(target: str, values: np.ndarray, sources: tuple, patches: dict) -> list[str]:
    cls, dims, _ = MAPPING[target]
    vector = cls == "volVectorField"
    kind = "vector" if vector else "scalar"
    zero = "(0 0 0)" if vector else "0"
    lines = [
        BANNER,
        f"// converted by {GENERATOR} from {', '.join(sources)}.",
        "// Do not edit: rerun the converter, which keeps identical values.",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        f"    class       {cls};",
        f"    object      {target};",
        "}",
        SEPARATOR,
        "",
        f"dimensions      {dims};",
        "",
        f"internalField   nonuniform List<{kind}>",
        str(len(values)),
        "(",
        *_format(values, vector),
        ")",
        ";",
        "",
        "boundaryField",
        "{",
    ]
    for name, info in patches.items():
        ptype = getattr(info, "type", "patch")
        if ptype == "empty":
            lines += [f"    {name}", "    {", "        type            empty;", "    }"]
        else:
            lines += [f"    {name}", "    {", "        type            calculated;",
                      f"        value           uniform {zero};", "    }"]
    lines += ["}", "", "// ************************************************************************* //"]
    return lines


def _same(path: Path, values: np.ndarray, n_cells: int) -> bool:
    text = path.read_text(encoding="utf-8", errors="replace")
    if GENERATOR not in text:
        raise ConversionError(
            f"{path} exists and was not written by {GENERATOR}; refusing to "
            f"replace it")
    try:
        existing = np.asarray(read_internal(path).over(n_cells, str(path)), dtype=np.float64)
    except FoamDataError as exc:
        raise ConversionError(f"{path} exists and cannot be read: {exc}") from exc
    if existing.shape != values.shape or not np.array_equal(existing, values):
        raise ConversionError(
            f"{path} exists with different values; the converter never "
            f"overwrites. Remove it by hand if the source fields really changed.")
    return True


def _n_cells(case_dir: Path) -> int:
    owner = case_dir / "constant" / "polyMesh" / "owner"
    match = re.search(r"nCells:\s*(\d+)", owner.read_text(encoding="utf-8", errors="replace"))
    if not match:
        raise ConversionError(f"{owner}: no 'nCells:' note in the header")
    return int(match.group(1))


def convert_time(time_dir: Path, *, species_mass_kg: float, average_start_s: float,
                 fields, patches: dict, n_cells: int | None = None) -> ConversionResult:
    """Write the dsmcFoam fields of one reconstructed time directory.

    ``species_mass_kg`` is recorded for symmetry with ``read_case``; the
    conversion itself needs no mass because uniGasFoam writes ``rhoM``.
    """
    del species_mass_kg
    time_dir = Path(time_dir)
    if n_cells is None:
        n_cells = _n_cells(time_dir.parent)
    result = ConversionResult(time=time_dir.name)
    sampled = float(time_dir.name) > average_start_s * (1.0 + TIME_TOLERANCE)
    cache: dict = {}

    def read(name: str) -> np.ndarray:
        if name not in cache:
            cache[name] = _values(time_dir, name, n_cells)
        return cache[name]

    for target, sources in _sources(fields).items():
        if MAPPING[target][2] and not sampled:
            continue
        if not all((time_dir / s).is_file() for s in sources):
            result.skipped.append(target)
            continue
        values = _compute(target, fields, read)
        path = time_dir / target
        if path.exists():
            _same(path, values, n_cells)
            result.kept.append(target)
            continue
        write_lines(path, _render(target, values, sources, patches))
        result.written.append(path)
    return result


def time_directories(case_dir: Path) -> list[Path]:
    """Reconstructed numeric time directories, except ``0``, oldest first."""
    found = []
    for path in Path(case_dir).iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0.0:
            found.append((value, path))
    return [p for _, p in sorted(found)]


def convert_case(case_dir: Path) -> list[ConversionResult]:
    """Convert every reconstructed time directory of a cai2012-unigas case."""
    from plumetools.cai2012ext import direct_case
    from plumetools.cai2012ugf.settings import load_settings
    from plumetools.mesh.boundary import read_boundary

    case_dir = Path(case_dir)
    settings = load_settings(case_dir)
    case = direct_case.derive(case_dir)
    patches = read_boundary(case_dir)
    n_cells = _n_cells(case_dir)
    return [convert_time(t, species_mass_kg=case.cfg.gas.mass_kg,
                         average_start_s=case.run.average_start_s,
                         fields=settings.fields, patches=patches, n_cells=n_cells)
            for t in time_directories(case_dir)]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m plumetools.cai2012ugf.convert",
        description="Write dsmcFoam-named fields from a uniGasFoam case's output.")
    parser.add_argument("case_dir", type=Path)
    args = parser.parse_args(argv)
    try:
        results = convert_case(args.case_dir)
    except ConversionError as exc:
        print(f"convert: {exc}", file=sys.stderr)
        return 1
    for r in results:
        print(f"  {r.time:>14}  wrote {len(r.written)}  kept {len(r.kept)}"
              + (f"  skipped {', '.join(r.skipped)}" if r.skipped else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
