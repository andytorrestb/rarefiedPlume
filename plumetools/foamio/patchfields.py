"""Patch entries of a field file: which boundary condition each patch type takes.

Shared by every writer of ``0/`` and of diagnostic fields. Four identical copies
of :data:`CONSTRAINT_PATCH_TYPES` and :func:`patch_field_type` lived in
``foamio.fields``, ``cai2012.foamfields``, ``markelov1999.foamfields`` and
``nozzleplume.foamfields``; each of those names now refers to these.
"""

from __future__ import annotations

from typing import Iterable

#: Geometric patch types whose patchField type must match exactly. OpenFOAM
#: rejects anything else with "inconsistent patch and patchField types".
CONSTRAINT_PATCH_TYPES = frozenset({
    "symmetry", "symmetryPlane", "empty", "wedge", "cyclic", "cyclicAMI",
    "processor", "processorCyclic", "nonuniformTransformCyclic",
})


def patch_field_type(geometric_type: str) -> str:
    """The patchField type a given geometric patch type requires.

    Constraint patches -- ``symmetry``, ``empty``, ``wedge``, ``cyclic`` and
    friends -- must carry a patchField of the same name; OpenFOAM fails with
    "inconsistent patch and patchField types" otherwise. Every other patch gets
    the neutral ``zeroGradient``.
    """
    return geometric_type if geometric_type in CONSTRAINT_PATCH_TYPES else "zeroGradient"


def nonuniform_value_entry(patch: str, list_type: str, rows: Iterable[str]) -> list[str]:
    """A ``fixedValue`` patch entry carrying one value per face.

    Args:
        patch: the patch name.
        list_type: ``scalar`` or ``vector``.
        rows: one formatted value per face, in patch-face order, each already
            indented as it is to appear.

    Returns:
        The entry's lines, ready for a ``boundaryField`` block.
    """
    rows = list(rows)
    return [
        f"    {patch}",
        "    {",
        "        type            fixedValue;",
        f"        value           nonuniform List<{list_type}>",
        f"        {len(rows)}",
        "        (",
        *rows,
        "        );",
        "    }",
    ]
