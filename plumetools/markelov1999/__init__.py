"""The AIAA 99-3455 flat-plate-in-a-cylinder-wake case family.

Lumpkin, F. E., Stewart, B. D., and Markelov, G. N., "Study of 3D Rarefied Flow
on a Flat Plate in the Wake of a Cylinder," AIAA 99-3455, 1999.

Everything in this package is **new work**, not a refactor of the legacy plume
model. It shares :mod:`plumetools.mesh`, :mod:`plumetools.geometry` and the
generic dictionary writers in :mod:`plumetools.foamio`, but it deliberately does
**not** reuse :mod:`plumetools.sourceflow`: that module is frozen bit-for-bit to
reproduce known defects (SM-01 .. SM-10), and the paper-faithful equations
disagree with it in five separate places. See
:mod:`plumetools.markelov1999.sourceflow` for the itemised list.

Module map
----------

===========================  ==================================================
``constants``                CODATA 2018 constants -- *not* the legacy 4-5 s.f.
                             values in :mod:`plumetools.constants`
``sourceflow``               the equations as printed in AIAA 99-3455
``geometry``                 cylinder / plate / gap / inflow coordinates
``mesh``                     blockMesh background + snappyHexMesh primitives
``foamfields``               the ``0/`` field files this case family needs
``foamdicts``                ``dsmcProperties``, ``controlDict``, and friends
``inflow``                   evaluate the model over the meshed inflow patch
``flux``                     mass / momentum / energy conservation checks
``resolution``               particles-per-cell estimate and post-run audit
``checks``                   DSMC mesh and time-step quality checks
``postprocess``              mesh stats, particle stats, wall pressures
``study``                    ``study.yaml`` parsing, case matrix, manifest
``verify``                   verify a generated mesh against ``case.yaml``
===========================  ==================================================
"""

from __future__ import annotations

#: Bibliographic record, carried into every generated manifest and case summary.
PAPER = {
    "reference": "AIAA 99-3455",
    "title": "Study of 3D Rarefied Flow on a Flat Plate in the Wake of a Cylinder",
    "authors": "Lumpkin, Stewart, Markelov",
    "year": 1999,
}
