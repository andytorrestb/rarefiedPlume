"""The Cai & Wang 2012 collisionless circular-plume validation family.

Cai, C., and Wang, L., "Numerical Validations for a Set of Collisionless Rocket
Plume Solutions," *Journal of Spacecraft and Rockets*, Vol. 49, No. 1, 2012,
pp. 59-68.

The problem is a **round free jet expanding from a circular nozzle exit directly
into vacuum**. The exit is a uniform drifting Maxwellian at ``(n0, U0, T0)`` on
the disk ``r <= D/2`` in the plane ``x = 0``; the plume axis is ``+x``.

This is not the repository's hemispherical analytical source flow. There is no
source-flow boundary anywhere in this package: the *physical nozzle exit* is the
DSMC inlet, which is the whole point of the case.

What is new here, and what is borrowed
--------------------------------------
Everything in this package is new work. It reuses, unchanged:

* :mod:`plumetools.mesh` -- ``constant/polyMesh`` parsing;
* :mod:`plumetools.foamio.primitives` -- the ``FoamFile`` banner, the LF-only
  writer and the length formatter;
* ``applications/dsmcBoundaryModels/plumeFieldInflow`` -- the inflow model.

It deliberately does **not** reuse :mod:`plumetools.sourceflow` (frozen
bit-for-bit against a regression golden, and a different physical model
entirely), :mod:`plumetools.config` (whose schema cross-validates a
source-flow sphere radius and a throat radius this case does not have), or
:mod:`plumetools.foamio.fields` (pinned to the legacy formatting).

Module map
----------

====================  =========================================================
``constants``         CODATA 2018 / SI-2019 exact constants
``gas``               argon VHS properties; Kn <-> lambda <-> n0; S0 <-> U0
``analytical``        Cai's collisionless circular-exit solution
``config``            ``case.yaml`` schema, loading and validation
``geometry``          nozzle disk, domain extents, X/D normalisation
``mesh``              graded multi-block Cartesian ``blockMeshDict``,
                      ``topoSetDict`` / ``createPatchDict`` for the nozzle disk,
                      and the mesh-cost diagnostics
``inflow``            the derived exit state, and the flux it implies
``foamfields``        the ``0/`` fields this family needs
``dictionaries``      ``dsmcProperties``, ``controlDict`` and friends
``study``             ``study.yaml`` parsing, the Kn matrix crossed with the
                      numerical-particle matrix, and the manifest
``checks``            DSMC mesh / time-step / occupancy quality checks
``verify``            whether the GENERATED tree is the study that describes it
``post``              centreline extraction, density plane, error metrics
====================  =========================================================

Two axes
--------
The family sweeps ``Kn`` (the physical matrix, **[PAPER]**) crossed with a
numerical-particle multiplier (1x, 2x, 5x, 10x parcels). A level changes
``nEquivalentParticles = baseline / multiplier`` and nothing else -- not the
mesh, not ``deltaT``, not the density -- so what it measures is statistical
resolution rather than an accident. :mod:`plumetools.cai2012.verify` enforces
that on the generated dictionaries before the cases are submitted.

Provenance labels
-----------------
Every physical value carries one of three labels, in ``case.yaml``, in the
docstrings and in ``docs/cai2012-case.md``:

``[PAPER]``
    printed in Cai & Wang 2012.
``[DERIVED]``
    computed from ``[PAPER]`` values; never written as a literal.
``[ASSUMPTION]``
    not in the paper. An implementation choice, with its reasoning recorded.

Cai used an **axisymmetric** DSMC code (GRASP) on a uniform grid. This is
**3-D Cartesian OpenFOAM** ``dsmcFoam``. That difference is structural, not
cosmetic, and is documented rather than hidden.
"""

from __future__ import annotations

#: Bibliographic record, carried into every manifest and case summary.
PAPER = {
    "reference": "Cai and Wang 2012",
    "title": "Numerical Validations for a Set of Collisionless Rocket Plume Solutions",
    "authors": "Cai, Wang",
    "journal": "Journal of Spacecraft and Rockets",
    "volume": 49,
    "number": 1,
    "year": 2012,
}

#: The model identifier a ``case.yaml`` in this family must declare.
MODEL = "cai2012_circular_plume"

#: Maximum centreline density errors Cai reports between his analytical solution
#: and his DSMC results, in percent. **Reference values, not pass/fail
#: tolerances**: this implementation uses a different DSMC code and a Cartesian
#: mesh, so reproducing the numbers exactly is not the claim. The *trend* is what
#: is being validated -- see ``docs/cai2012-case.md``.
CAI_REPORTED_MAX_DENSITY_ERROR_PERCENT = {
    100.0: 0.07,
    0.1: 3.63,
    0.01: 12.03,
}
