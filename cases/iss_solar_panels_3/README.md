# cases/iss_solar_panels_3 — legacy solar-panel impingement case

A copy of the 3d-inflow source-flow case set up for a notional panel
(commit `115d4f5`, 2023-02-28). It still runs the pre-refactor
`processInflowData.py` with inputs typed into the script, names the MNF fork's
`dsmcFoam+` in `system/controlDict`, and has no run script, so **it cannot be
run with the supported toolchain** as it stands. No result here is validated.

> **Do not run `clean.sh`.** It is `rm -r 0* boundaries fieldMeasurements
> postProcessing processor*` with no `cd` guard: run from anywhere else, it
> deletes matching paths there (finding RP-01).

Variants, blockers and what can safely be inspected:
[`docs/tutorials/iss-solar-panels.md`](../../docs/tutorials/iss-solar-panels.md).
