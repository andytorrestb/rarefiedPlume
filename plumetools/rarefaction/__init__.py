r"""Local mean free path, Knudsen numbers and mesh resolution, from a DSMC sample.

What this is for
----------------
Deciding, after a run, where the **next** mesh needs to be finer. It reads the
time-averaged fields one ``dsmcFoam`` run sampled and writes, cell by cell:

=====================  ====================================================
field                  meaning
=====================  ====================================================
``mfp``                mean free path ``lambda`` [m], equation (G1)
``Kn<L>``              ``lambda / L`` for each reference length -- ``KnD``,
                       on the nozzle exit diameter, for a nozzle-plume case
``cellSizeV``          ``V^(1/3)`` [m]
``cellSizeMax``        largest bounding-box dimension ``Delta_max`` [m]
``KnCellV``,           ``lambda / Delta``: resolution, as a Knudsen number
``KnCellMax``
``cellOverMfpV``,      ``Delta / lambda``: the same, inverted; above ``1/k``
``cellOverMfpMax``     (``k = 3``) the cell is coarser than the target
``KnGradN``            ``lambda |grad n| / n``
``KnGradT``            ``lambda |grad T| / T``
``cellSizeTarget``     ``lambda / k`` [m], never clipped
``mfpResolved``        the criterion as a map: 1 where ``lambda / Delta_max
                       >= k``, 0 where not, -1 where there is no value
``parcelsPerStep``     simulated parcels in the cell, per step (dsmcRhoNMean)
``sampledParcels``     parcel samples over the averaging window
``diagStatus``         0 valid, else why not (flags; see :mod:`.core`)
``diagValid``,         1 where the cell's values, and its gradients, can be
``diagValidGradN``,    trusted
``diagValidGradT``
=====================  ====================================================

and a summary -- valid fractions, resolution statistics, where the worst
violations are, which octree levels a practical minimum cell size allows, and
what a mesh that resolved ``lambda / k`` would cost.

The **cell** numbers measure the mesh: how many mean free paths one cell spans.
``KnD`` and the gradient numbers measure the **flow**: how rarefied it is on the
nozzle's length scale and on the scales its own density and temperature change
over. None of them, alone or against a threshold, says a DSMC solution is
valid.

Nothing here changes a mesh or a run. It is the "inspect" step of
run -> inspect -> remesh -> start fresh; :mod:`plumetools.iteration` is the
"start fresh".

Layout
------

====================  ==========================================================
:mod:`.core`          the numerics and the diagnostic registry. numpy only
:mod:`.settings`      thresholds, lengths, field names: the ``diagnostics:``
                      schema
:mod:`.foam`          read a sampled time directory; write the fields back
:mod:`.summary`       statistics, violations, refinement levels, the budget
:mod:`.run`           the three together, for one case and one time
``__main__``          ``python -m plumetools.rarefaction`` for any dsmcFoam case
====================  ==========================================================

Enabling it for a case
----------------------
**A nozzle-plume case** (``model: nozzle_exit_plume``) has it already: its
``./Allpost`` runs ``python -m plumetools.nozzleplume diagnostics``, which
takes the gas, the mean-free-path convention and the exit diameter from
``case.yaml``. An optional ``diagnostics:`` section there changes the
thresholds; leaving it out means the defaults in :mod:`.settings`.

**Any other dsmcFoam case** needs only a settings file::

    # diagnostics.yaml
    gas: {molecular_mass_kg: 6.63e-26, molecular_diameter_m: 4.17e-10,
          omega: 0.81, t_ref_K: 273.0, mean_free_path_convention: vhs}
    background_cell_size_m: null        # set it for an octree mesh
    diagnostics:
      reference_lengths_m: {D: 0.2}     # -> KnD

    postProcess -func writeCellCentres -latestTime
    postProcess -func writeCellVolumes -latestTime
    python -m plumetools.rarefaction <case> --settings diagnostics.yaml

or, from Python, :func:`plumetools.rarefaction.run.diagnose`.

Adding a diagnostic
-------------------
Register it in :mod:`.core` -- one decorated function naming its inputs; see
"Adding a diagnostic" there. It is then computed, written, summarised (a name
starting ``Kn`` lands in the summary's ``rarefaction`` block) and selectable
through ``diagnostics.select`` with no other change. A new input is read in
:func:`plumetools.rarefaction.foam.load`; until it is, the diagnostic is
skipped with the reason, and the others run.

See ``docs/rarefaction-diagnostics.md``.
"""
