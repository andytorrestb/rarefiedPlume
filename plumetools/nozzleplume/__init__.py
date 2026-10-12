"""Nozzle-exit plumes into vacuum, with or without a plate in the way.

The general form of the problem ``plumetools.cai2012`` solves for one paper: a
**uniform drifting Maxwellian on a circular nozzle exit**, expanding into a
vacuum box. The exit is the disk ``r <= D/2`` in the plane ``x = 0``, the plume
axis is ``+x``, and the *physical nozzle exit* is the DSMC inlet -- there is no
analytical source surface standing in for it.

What is general about it
------------------------
``plumetools.cai2012`` is welded to its paper in four places, and each one is a
thing a second nozzle would need changed:

==========================  ===================================================
``cai2012``                 here
==========================  ===================================================
exit state from ``Kn``,     exit state given **directly** -- ``n`` (or ``rho``),
``S0`` and an assumed       ``U`` and ``T`` as a nozzle code or a data sheet
``T0``                      prints them. ``Kn`` and ``S`` are still accepted,
                            so a Cai case is expressible, and a unit test holds
                            the two packages to the same numbers for one.
argon, monatomic, VHS       any single VHS species; internal degrees of freedom
                            and Larsen-Borgnakke exchange are configuration.
vacuum on every boundary    an optional **flat plate** normal to the axis: a
                            ``wall`` the plume impinges on, with its own
                            temperature and gas-surface model.
a graded ``blockMesh``      a uniform cubic background refined by
                            ``snappyHexMesh`` -- which is what makes a body in
                            the plume possible at all.
==========================  ===================================================

``plumetools.cai2012`` itself is untouched and remains the record of that
validation; nothing here imports it. The duplication between the two is
deliberate, in the same way that package restates two physical constants rather
than import another family's: a paper-specific package that a general one
depended on could never be frozen.

Module map
----------

====================  =========================================================
``constants``         CODATA 2018 / SI-2019 exact constants
``gas``               one VHS species; mean free path, speeds, wall fluxes
``config``            ``case.yaml`` schema, loading and validation
``exitstate``         the exit state, and everything it implies
``geometry``          nozzle disk, plate, domain box, symmetry planes
``mesh``              background ``blockMeshDict``, ``snappyHexMeshDict``,
                      the nozzle-disk ``topoSetDict`` / ``createPatchDict``,
                      and the mesh plan they follow from
``run``               time step, particle weight, durations, parcel estimate
``dictionaries``      ``dsmcProperties``, ``controlDict`` and friends
``foamfields``        the ``0/`` fields
``checks``            DSMC mesh / time-step / occupancy quality checks
``verify``            whether the built mesh is the mesh that was planned
``post``              plate loads, centreline profiles
``workflow``          the steps a case's ``Allmesh`` / ``Allrun`` call
====================  =========================================================

A case is a ``case.yaml`` and four thin scripts. Everything OpenFOAM reads is
generated from that one file::

    python -m plumetools.nozzleplume describe    what case.yaml implies; writes nothing
    python -m plumetools.nozzleplume setup       system/ and constant/ dictionaries
    python -m plumetools.nozzleplume verify      the mesh against case.yaml
    python -m plumetools.nozzleplume fields      0/
    python -m plumetools.nozzleplume post        results/

Provenance labels
-----------------
``case.yaml`` tags every physical value:

``[GIVEN]``
    supplied for the case being modelled -- a data sheet, a nozzle code.
``[DERIVED]``
    computed from ``[GIVEN]`` values; printed, never typed in.
``[ASSUMPTION]``
    not supplied. A modelling choice, with its reasoning recorded.
"""

from __future__ import annotations

#: The model identifier a ``case.yaml`` in this family must declare.
MODEL = "nozzle_exit_plume"
