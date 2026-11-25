"""Decomposition weights: balance the expected particles, not the cells.

uniGasFoam's cost follows particles and cells. The cells crowd near the nozzle
(they are sized by the mean free path) while the particles spread through the
plume and pile up in front of the plate, so scotch balancing CELLS put 90% of
the smoke run's particles on one rank. ``decomposePar`` balances a cell field
instead when ``decomposeParDict`` names one (``weightField``); this writes it:

    weight = CELL_COST + n V / W

with ``n`` the collisionless plume-plus-plate-gas density
(:mod:`plumetools.plateugf.collisionless`), ``V`` the cell volume and ``W``
``nEquivalentParticles``. On the Kn_e = 100 prototypes the prediction matched
the measured particles per cell (median ratio 1.00-1.03), and the split it gave
was as even as one from the measured counts. Where cells outnumber particles
(Kn_e = 1e-4) the weights fall back towards a plain cell split by themselves.

The weights do not depend on time: ``Allmesh`` writes them at ``0`` and
``Allrun`` copies them to any later time it re-decomposes.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from plumetools.foamio.volfield import read_internal_field
from plumetools.nozzleplume.config import load_case_config
from plumetools.plateugf.collisionless import Model

FIELD = "decompositionWeight"
#: One cell's work per step, in particles' worth (a call, not a measurement).
CELL_COST = 1.0
_DIVIDER = "// " + "*" * 73 + " //"


class WeightError(RuntimeError):
    """The weights cannot be written as asked."""


def particle_weight(case_dir: Path) -> float:
    """``nEquivalentParticles`` from ``constant/uniGasProperties``."""
    path = Path(case_dir) / "constant" / "uniGasProperties"
    match = re.search(r"^\s*nEquivalentParticles\s+([^;\s]+)\s*;",
                      path.read_text(encoding="utf-8"), re.M)
    if not match:
        raise WeightError(f"no nEquivalentParticles in {path}")
    return float(match.group(1))


def cell_weights(case_dir: Path, time: str = "0") -> np.ndarray:
    """``CELL_COST`` + the collisionless particles of every cell, from ``<time>/C`` and ``V``."""
    case_dir = Path(case_dir)
    t = case_dir / time
    for name in ("C", "V"):
        if not (t / name).is_file():
            raise WeightError(f"no {t / name}: run 'postProcess -func writeCellCentres -time {time}' "
                              f"and 'postProcess -func writeCellVolumes -time {time}' first")
    centres = np.asarray(read_internal_field(t / "C"), float).reshape(-1, 3)
    volumes = np.asarray(read_internal_field(t / "V"), float).ravel()
    density = Model(load_case_config(case_dir)).density(centres)
    return CELL_COST + density * volumes / particle_weight(case_dir)


def write_field(path: Path, values: np.ndarray) -> Path:
    """A dimensionless ``volScalarField`` named after ``path``, ending on OpenFOAM's divider."""
    path = Path(path)
    head = (
        "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
        f"    class       volScalarField;\n    object      {path.name};\n}}\n\n"
        "// Written by python -m plumetools.plateugf weights: decomposePar's weightField.\n\n"
        f"dimensions      [0 0 0 0 0 0 0];\n\ninternalField   nonuniform List<scalar>\n{len(values)}\n(\n")
    tail = (")\n;\n\nboundaryField\n{\n    #includeEtc \"caseDicts/setConstraintTypes\"\n"
            "    \".*\" { type calculated; value uniform 1; }\n}\n\n" + _DIVIDER + "\n")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(head)
        np.savetxt(handle, np.asarray(values, float), fmt="%.6g")
        handle.write(tail)
    return path


def write_weights(case_dir: Path, time: str = "0") -> Path:
    """Write ``<time>/decompositionWeight``; print what it balances."""
    case_dir = Path(case_dir)
    weights = cell_weights(case_dir, time)
    path = write_field(case_dir / time / FIELD, weights)
    particles = float((weights - CELL_COST).sum())
    print(f"    {path}: {weights.size:,} cells, {particles:,.0f} expected particles "
          f"(collisionless plume + plate gas), weight {weights.min():.3g} .. {weights.max():.3g}")
    return path
