"""One case, one sampled time: read, compute, write, summarise.

The whole pipeline in one call, shared by ``python -m plumetools.rarefaction``
and the nozzle-plume family's ``diagnostics`` command. A case family adds what
only it knows -- its reference lengths, its background cell, a centreline --
around this, not inside it.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from plumetools.rarefaction import core, foam, summary


@dataclass
class Diagnosis:
    """What :func:`diagnose` produced.

    Attributes:
        data: the inputs read (:class:`~plumetools.rarefaction.foam.CaseData`).
        context: the :class:`~plumetools.rarefaction.core.Context`.
        outputs: ``{field name: Field}`` computed.
        skipped: ``{field name: reason}`` not computed.
        document: the summary.
        written: field files written into the time directory.
    """

    data: object
    context: object
    outputs: dict
    skipped: dict
    document: dict
    written: list


def file_sha256(path: Path) -> str | None:
    """The SHA-256 of a file, for provenance; ``None`` if it is absent."""
    path = Path(path)
    if not path.is_file():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def diagnose(case_dir: Path, species, settings, *, convention: str = "vhs",
             time: str | None = None, reference_lengths: dict | None = None,
             background_cell_size_m: float | None = None,
             header: dict | None = None, write_fields: bool = True) -> Diagnosis:
    """Compute, write and summarise the diagnostics of one sampled time.

    Args:
        case_dir: a reconstructed ``dsmcFoam`` case.
        species: the gas -- ``diameter_m``, ``omega``, ``t_ref_K`` and
            ``gas_constant_j_per_kg_k``, as
            :class:`plumetools.nozzleplume.gas.VhsSpecies` has.
        settings: :class:`~plumetools.rarefaction.settings.DiagnosticSettings`.
        convention: the mean-free-path convention of the case.
        time: the time directory; the latest by default.
        reference_lengths: ``{name: metres}`` beyond
            ``settings.reference_lengths_m`` -- what the case family supplies.
        background_cell_size_m: the level-0 cell of an octree mesh, for the
            refinement-level planning.
        header: case identification for the top of the summary.
        write_fields: write the volume fields into the time directory.

    Raises:
        plumetools.rarefaction.foam.PostError: no time, no mesh, mixed
            averaging windows, or a field from another mesh.
        ValueError: a ``select`` entry that is not an output.
    """
    case_dir = Path(case_dir)
    data = foam.load(case_dir, settings, time=time)
    lengths = {**dict(reference_lengths or {}), **settings.reference_lengths_m}
    extra = {f"Kn{name}": core.knudsen_on(name, length)
             for name, length in lengths.items()}
    ctx = core.Context(data.inputs, species=species, settings=settings,
                       convention=convention, unavailable=data.unavailable,
                       extra=extra)
    names = core.output_names(ctx, settings.select)
    outputs, skipped = core.evaluate(ctx, names)

    written = (foam.write_fields(data, outputs, ctx.registry,
                                 settings.invalid_value) if write_fields else [])
    document = summary.build(
        outputs, ctx, data, settings, skipped=skipped,
        background_cell_size_m=background_cell_size_m,
        header={
            "generator": "plumetools.rarefaction",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "case_dir": str(case_dir.resolve()),
            "mesh_cells": data.n_cells,
            "fields_written": [path.name for path in written],
            **(header or {}),
        })
    return Diagnosis(data=data, context=ctx, outputs=outputs, skipped=skipped,
                     document=document, written=written)
