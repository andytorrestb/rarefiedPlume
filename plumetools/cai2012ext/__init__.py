r"""The Cai & Wang 2012 Knudsen extension and explicit seven-decade suite.

Cai, C., and Wang, L., "Numerical Validations for a Set of Collisionless Rocket
Plume Solutions", *J. Spacecraft and Rockets*, Vol. 49, No. 1, 2012.

``cases/cai2012`` runs Cai's own three Knudsen numbers. This extension runs the
same circular argon nozzle -- ``D = 0.2 m``, ``S0 = 2``, ``T0 = 300 K``, a uniform
drifting Maxwellian into vacuum, collisions on -- at four Knudsen numbers Cai did
**not** report, each at ONE numerical-particle level, the baseline:

.. code-block:: text

    Kn          1e-4   1e-3   1   10          (EXTENSION_KNUDSEN)
    parcels     1x -- 40 in the exit cell     (PARTICLE_MULTIPLIERS)

These four configurations remain the legacy profile. There is no parcel-population sweep: the
baseline is cases/cai2012-health's highest-resolution weight, and the statistical
health of each run is judged from its own sampled output (occupancy audit,
stationarity, the rarefaction diagnostics) rather than from sibling runs at other
weights. Cai's three -- 0.01, 0.1 and 100 -- are
:data:`EXCLUDED_KNUDSEN`: they are the conditions his paper already compares,
and the legacy profile does not generate, mesh, schedule or post-process a case
at one of them. The explicit ``meta.knudsen_suite: full_decades`` profile adds
all seven decade cases, including those paper Kn values, using study-local
generated cases. It does not import other studies' results or relax numerical
feasibility requirements.

Nothing physical lives here
---------------------------
Everything that decides what a case *is* comes from :mod:`plumetools.cai2012`
unchanged: the config schema and its validation, ``Kn -> n0``, the mesh plan, the
particle weight, the time step, the dictionaries, the checks, the analytical
solution and the centreline post-processing. ``cases/cai2012-kn-extension``
writes its ``study.yaml`` in that family's own schema and its cases are ordinary
``cai2012_circular_plume`` cases.

What this package adds is what an extension needs and the parent study did not:

==================  =========================================================
:mod:`.policy`      the fixed matrix, the excluded Kn, the inherited settings
:mod:`.provenance`  every derived exit condition, recomputed independently
:mod:`.feasibility` whether a case is fit for production, and what it costs --
                    the nozzle-plume family's pre-run resolution checks
:mod:`.generation`  idempotent generation that never overwrites results
:mod:`.hpc`         what a case's PBS job does around ``./Allrun``: the Kn
                    guard, the NON-PRODUCTION label, the locked wall clock
:mod:`.diagnostics` the nozzle-plume family's post-run rarefaction and
                    mesh-resolution diagnostics (:mod:`plumetools.rarefaction`)
:mod:`.comparison`  angular, off-axis velocity, temperature and pressure
:mod:`.convergence` stationarity of sampled quantities across frames
==================  =========================================================

and the **interface track** (``docs/cai2012-kn-interface.md``), which removes the
dense core and injects a continuum continuation of the exit on a downstream
sphere:

=========================  ==================================================
:mod:`.interface`          the Markelov-source mapping, the closure, Bird's
                           B, the sphere, the per-face state (pure)
:mod:`.interface_mesh`     the parent's lattice with the sphere carved out
:mod:`.interface_case`     one derivation path for a case: assess, mesh,
                           inflow, identity, the PBS guard and label
:mod:`.interface_study`    ``interfaces.yaml`` and the ``Interface/`` tree
:mod:`.overlap`            truncated run against the full-domain benchmark
=========================  ==================================================

``plumetools.cai2012`` is not edited by any of it, nor is
``plumetools.markelov1999``, whose angular law the interface calls. Those
packages and their case families are the record of what was run, so this one
imports from them and leaves them alone.
"""

from __future__ import annotations

from plumetools.cai2012 import CAI_REPORTED_MAX_DENSITY_ERROR_PERCENT

#: Cai's own Knudsen numbers, excluded by the legacy extension profile.
#:
#: Taken from the table of his reported centreline errors rather than typed out
#: again, so the list the extension refuses and the list the parent study
#: compares against cannot drift apart.
EXCLUDED_KNUDSEN = tuple(sorted(float(k)
                                for k in CAI_REPORTED_MAX_DENSITY_ERROR_PERCENT))

#: The Knudsen numbers this extension runs -- exactly these, on the nozzle
#: diameter. The requested scientific range is 1e-4 .. 10 by decades; the two
#: members of it Cai already reports (0.01 and 0.1) are excluded above.
EXTENSION_KNUDSEN = (1.0e-4, 1.0e-3, 1.0, 10.0)

#: Explicit full-decade profile, ordered most rarefied first.
FULL_SUITE_KNUDSEN = (100.0, 10.0, 1.0, 0.1, 0.01, 0.001, 0.0001)

#: The numerical-particle levels, as multipliers on the baseline parcel
#: population. ``nEquivalentParticles = baseline / multiplier``.
#:
#: The baseline only. The parent study's machinery still names a case by its
#: level, which is why the directories are ``<Kn>_np1x``.
PARTICLE_MULTIPLIERS = (1.0,)

#: Cai refers his ``dx/lambda0``, ``dt/t0`` and his ``10^4 t0`` startup to the
#: Kn = 0.01 exit state. The extension keeps that reference for EVERY case.
#: Computing reference quantities at this Kn does not generate a case there.
REFERENCE_KNUDSEN = 0.01

#: Where the extension's cases live, for messages and provenance.
FAMILY_DIRECTORY = "cases/cai2012-kn-extension"

#: Recorded as ``meta.generated_by`` in every generated ``case.yaml``.
GENERATOR = f"{FAMILY_DIRECTORY}/generate_cases.py"

__all__ = [
    "EXCLUDED_KNUDSEN",
    "EXTENSION_KNUDSEN",
    "FULL_SUITE_KNUDSEN",
    "FAMILY_DIRECTORY",
    "GENERATOR",
    "PARTICLE_MULTIPLIERS",
    "REFERENCE_KNUDSEN",
]
