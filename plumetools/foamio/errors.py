"""The error the shared OpenFOAM readers raise."""

from __future__ import annotations


class FoamDataError(ValueError):
    """Required OpenFOAM data is missing, unreadable, or does not fit the mesh.

    Raised by :mod:`plumetools.foamio.volfield` and :mod:`plumetools.foamio.times`
    for a missing time directory or field file, an entry that cannot be parsed,
    and a field whose length is not the mesh's. It is a ``ValueError``, so a
    caller written against the older readers keeps catching it.

    ``nozzleplume.post.PostError`` and ``rarefaction.foam.PostError`` are this
    class under the names those modules have always raised;
    ``cai2012.post.PostError`` is a subclass of it.
    """
