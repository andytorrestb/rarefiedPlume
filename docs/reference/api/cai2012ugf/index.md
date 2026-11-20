# `plumetools.cai2012ugf`

The Cai Knudsen sweep on uniGasFoam's hybrid USP-DSMC solver
(`cases/cai2012-unigas`): the kn-extension's seven cases, case for case, on
another solver. It writes the uniGasFoam dictionaries, converts uniGasFoam's
fields to the dsmcFoam names the analysis reads (losslessly), judges each run
after it ends, runs the two-leg, cloud-safe restarts, and adapts the scaling
inventory and reports. It imports `cai2012` and `cai2012ext` and edits neither.

Depends on: `numpy`, `PyYAML`; `plumetools.cai2012`, `plumetools.cai2012ext`,
`plumetools.foamio`, `plumetools.mesh`, `plumetools.rarefaction`.
Workflow: [Kn sweep on uniGasFoam](../../../tutorials/cai2012-unigas.md).

::: plumetools.cai2012ugf
    options:
      show_root_heading: false
      heading_level: 3
