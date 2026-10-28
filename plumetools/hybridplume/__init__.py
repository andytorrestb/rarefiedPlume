"""Steady one-way hybrid CFD -> DSMC free plume, full 3-D (Cai & Wang 2012 source).

A new case family: it imports the shared OpenFOAM readers (``plumetools.mesh``,
``plumetools.foamio``), the generic diagnostics (``plumetools.rarefaction.core``),
the block-segment solver of ``plumetools.cai2012.mesh`` and, for cross-checks
only, ``plumetools.cai2012.analytical``; it edits none of them.

Modules, in the order a case uses them::

    physics        the source state, the gas, VHS-consistent transport
    config         hybrid.yaml: profiles, validation, provenance labels
    mesh, cfd      the 3-D box mesh and the rhoCentralFoam precursor
    resources      cells, parcels, memory, steps, runtime -- before launching
    meshgeom       whole-mesh geometry from polyMesh
    breakdown      Kn_rho, Kn_T, Kn_U, Kn_GLL, B_exp; the retained region
    interface      the fixed coupling surface and its suitability verdict
    halfrange      exact Maxwellian half-range fluxes
    fluxes         Navier-Stokes-Fourier flux densities
    mapping        CFD interface state -> DSMC inflow faces, verified
    dsmc           standard dsmcFoam + plumeFieldInflow dictionaries and fields
    particles      parcel snapshots: tensor moments, exact crossings
    post, budgets  profiles, azimuthal statistics, conservation
    collisionless  the independent collisionless reference
    report, hpc    summary, figures, batch scripts
    workflow       the explicit steps; ``python -m plumetools.hybridplume``

See ``docs/hybridplume.md``.
"""
