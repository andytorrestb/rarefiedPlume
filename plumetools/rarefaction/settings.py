"""How the rarefaction diagnostics are run: thresholds, lengths, field names.

One schema, read from two places. A nozzle-plume case keeps it as the optional
``diagnostics:`` section of its ``case.yaml``
(:mod:`plumetools.nozzleplume.config`); any other dsmcFoam case can hand the
same keys, plus a ``gas:`` block, to ``python -m plumetools.rarefaction``.
Every key has a default, so an absent section means "run them all with the
defaults below".

Unknown keys are an error, as everywhere else in this repository: a mistyped
threshold should stop the run, not fall back to a default nobody chose.
"""

from __future__ import annotations

import math
import re
from dataclasses import MISSING, dataclass, field, fields
from typing import Any

#: The averaged ``dsmcFoam`` fields each input is read from, and the cell
#: geometry ``postProcess`` writes. ``field_names`` in the settings overrides
#: any of these -- a run whose ``fieldAverage`` was named differently, or a
#: solver that writes other names, needs no code change.
DEFAULT_FIELD_NAMES = {
    "number_density": "rhoNMean",           # n [1/m^3]
    "mass_density": "rhoMMean",             # rho [kg/m^3]
    "momentum_density": "momentumMean",     # rho u [kg/(m^2 s)]
    "kinetic_energy_density": "linearKEMean",   # (1/2) rho <c^2> [J/m^3]
    "parcels_per_step": "dsmcRhoNMean",     # simulated parcels in the cell [-]
    "cell_centres": "C",                    # postProcess -func writeCellCentres
    "cell_volumes": "V",                    # postProcess -func writeCellVolumes
}

#: What ``min_sampled_particles`` is compared against.
#:
#: ``per_step``     the sampled -- time-averaged -- number of parcels in the
#:                  cell, ``dsmcRhoNMean``: how many DISTINCT parcels its
#:                  moments come from at any step. The default, and what
#:                  "particles per cell" means everywhere else in this
#:                  repository (``dsmc.target_particles_per_exit_cell``).
#: ``accumulated``  parcels per step x steps averaged: every parcel sample
#:                  behind the means. Far larger, and far weaker: successive
#:                  steps see the same parcels, so it counts a slow parcel
#:                  once per step it lingers. On a real run, cells sampled by
#:                  ONE parcel reached 26-149 "samples" this way.
PARTICLE_COUNT_BASES = ("per_step", "accumulated")

_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


class SettingsError(ValueError):
    """The diagnostics settings are malformed or inconsistent."""


@dataclass(frozen=True)
class DiagnosticSettings:
    """Settings for :mod:`plumetools.rarefaction`.

    Attributes:
        enabled: run the diagnostics at all.
        min_sampled_particles: a cell's diagnostics are VALID only with **more**
            than this many sampled parcels. A confidence flag on the numbers
            computed for that cell -- not a refinement trigger and not a
            statement that the run has converged.
        particle_count_basis: what the threshold counts; see
            :data:`PARTICLE_COUNT_BASES`.
        cells_per_mean_free_path: ``k`` in the target cell size ``lambda / k``.
            A cell is reported as under-resolved where ``Delta / lambda > 1/k``.
        practical_min_cell_size_m: the smallest nominal cell the next mesh may
            use -- a planning constraint. Targets below it are REPORTED, never
            clipped.
        mesh_budget_cells: the cell count a planned mesh is measured against;
            a planning constraint, not a limit anything enforces.
        reference_lengths_m: ``{name: L}``; each gives a field ``Kn<name>`` =
            ``lambda / L``. A nozzle-plume case adds ``D``, the exit diameter.
        select: output field names to write, or ``None`` for all of them.
        field_names: overrides of :data:`DEFAULT_FIELD_NAMES`.
        invalid_value: what a cell with nothing to compute from holds in a
            written field (an under-sampled cell keeps its value). Negative,
            because every diagnostic is >= 0, so a reader can never mistake it
            for a value -- and finite, because a ``nan`` is not something every
            OpenFOAM reader parses.
        n_worst_cells: how many worst resolution violations the summary lists.
        worst_cells_separation: listed violations are at least this many of
            their own cell sizes apart, so the list is not ten neighbours.
    """

    enabled: bool = True
    min_sampled_particles: float = 20.0
    particle_count_basis: str = "per_step"
    cells_per_mean_free_path: float = 3.0
    practical_min_cell_size_m: float = 0.010
    mesh_budget_cells: int = 10_000_000
    reference_lengths_m: dict = field(default_factory=dict)
    select: tuple | None = None
    field_names: dict = field(default_factory=dict)
    invalid_value: float = -1.0
    n_worst_cells: int = 10
    worst_cells_separation: float = 4.0

    @property
    def target_ratio(self) -> float:
        """``1 / cells_per_mean_free_path``: the largest acceptable ``Delta / lambda``."""
        return 1.0 / float(self.cells_per_mean_free_path)

    def input_field(self, role: str) -> str:
        """The file name an input is read from."""
        return {**DEFAULT_FIELD_NAMES, **self.field_names}[role]

    def as_dict(self) -> dict:
        """Plain values, for the summary's record of how it was computed."""
        return {
            "min_sampled_particles": float(self.min_sampled_particles),
            "particle_count_basis": self.particle_count_basis,
            "cells_per_mean_free_path": float(self.cells_per_mean_free_path),
            "target_cell_over_mean_free_path": self.target_ratio,
            "practical_min_cell_size_m": float(self.practical_min_cell_size_m),
            "mesh_budget_cells": int(self.mesh_budget_cells),
            "reference_lengths_m": {k: float(v)
                                    for k, v in self.reference_lengths_m.items()},
            "select": None if self.select is None else list(self.select),
            "field_names": {**DEFAULT_FIELD_NAMES, **self.field_names},
            "invalid_value": float(self.invalid_value),
        }

    @classmethod
    def from_mapping(cls, data: Any, where: str = "diagnostics"
                     ) -> "DiagnosticSettings":
        """Build and validate from a YAML mapping (``None`` -> all defaults).

        Raises:
            SettingsError: naming the key at fault.
        """
        if data is None:
            data = {}
        if not isinstance(data, dict):
            raise SettingsError(f"{where}: expected a mapping, got "
                                f"{type(data).__name__}")
        known = {f.name for f in fields(cls)}
        unknown = sorted(set(data) - known)
        if unknown:
            raise SettingsError(f"{where}: unknown key(s) {unknown}; valid keys "
                                f"are {sorted(known)}")
        values = {}
        for f in fields(cls):
            if data.get(f.name) is None:        # null -> the default
                continue
            value = data[f.name]
            default = f.default if f.default is not MISSING else f.default_factory()
            if isinstance(default, bool) != isinstance(value, bool):
                raise SettingsError(
                    f"{where}: {f.name} is {value!r}; "
                    + ("expected true or false" if isinstance(default, bool)
                       else "a boolean is not a number"))
            if isinstance(default, (int, float)) and not isinstance(default, bool):
                try:
                    value = type(default)(float(value)) if isinstance(default, int) \
                        else float(value)
                except (TypeError, ValueError):
                    raise SettingsError(
                        f"{where}: {f.name} is {value!r}, which is not a number. "
                        f"YAML reads 1.0e7 as a string; write 1.0e+7.") from None
            values[f.name] = value
        settings = cls(**values)
        settings._validate(where)
        return settings

    def _validate(self, where: str) -> None:
        def fail(message: str):
            raise SettingsError(f"{where}: {message}")

        for name in ("min_sampled_particles",):
            if not (math.isfinite(getattr(self, name)) and getattr(self, name) >= 0):
                fail(f"{name} must be a number >= 0")
        for name in ("cells_per_mean_free_path", "practical_min_cell_size_m",
                     "worst_cells_separation"):
            value = getattr(self, name)
            if not (math.isfinite(value) and value > 0.0):
                fail(f"{name} must be finite and positive, got {value}")
        for name in ("mesh_budget_cells", "n_worst_cells"):
            if int(getattr(self, name)) < 0:
                fail(f"{name} must be >= 0")
        if self.particle_count_basis not in PARTICLE_COUNT_BASES:
            fail(f"particle_count_basis is {self.particle_count_basis!r}; valid: "
                 f"{list(PARTICLE_COUNT_BASES)}")
        if not (math.isfinite(self.invalid_value) and self.invalid_value < 0.0):
            fail(f"invalid_value is {self.invalid_value}; it must be finite and "
                 f"negative, so it cannot be mistaken for a diagnostic (all of "
                 f"which are >= 0)")

        if not isinstance(self.reference_lengths_m, dict):
            fail("reference_lengths_m must be a mapping of name: metres")
        for name, length in self.reference_lengths_m.items():
            if not _NAME.match(str(name)):
                fail(f"reference_lengths_m: {name!r} becomes the field name "
                     f"Kn{name}; use letters, digits and underscores")
            try:
                length = float(length)
            except (TypeError, ValueError):
                fail(f"reference_lengths_m.{name} is {length!r}, not a length")
            if not (math.isfinite(length) and length > 0.0):
                fail(f"reference_lengths_m.{name} must be positive")

        if self.select is not None:
            if not isinstance(self.select, (list, tuple)) or not all(
                    isinstance(name, str) for name in self.select):
                fail("select must be a list of output field names, or null for all")
            object.__setattr__(self, "select", tuple(self.select))

        if not isinstance(self.field_names, dict):
            fail("field_names must be a mapping of input: field name")
        unknown = sorted(set(self.field_names) - set(DEFAULT_FIELD_NAMES))
        if unknown:
            fail(f"field_names: unknown input(s) {unknown}; valid: "
                 f"{sorted(DEFAULT_FIELD_NAMES)}")
