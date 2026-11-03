# Tutorials

One tutorial per case family. Each starts with an *at a glance* box —
purpose, model, evidence, runtime, starting directory, blockers — and then
follows the same order: purpose and status, geometry and coordinates,
prerequisites, inputs, cheap checks that need no OpenFOAM, the production run
(generate, mesh, run, resume, post-process), outputs and how to judge them,
resources, safe customisation, a figure placeholder, and links.

| Family | Runs today? | Cheapest meaningful check | Tutorial |
|---|---|---|---|
| `cai2012` | yes | `generate_cases.py --dry-run` | [cai2012](cai2012.md) |
| `cai2012-health` | yes, with CaseFoam | `generate_cases.py --dry-run` (no CaseFoam) | [cai2012-health](cai2012-health.md) |
| `cai2012-kn-extension` | yes; three Kn unsuitable for production | `generate_cases.py --dry-run`, `validate_cases.py --inputs-only` | [direct cases](cai2012-kn-extension.md) |
| — interface track | smoke case only | `generate_interfaces.py --dry-run` | [interface track](cai2012-interfaces.md) |
| — field scaling | yes (Python) | `Allscale.py inventory`, `smoke` | [field scaling](cai2012-scaling.md) |
| `cai2012-virtual-collision` | yes (Python) | `assess.py inventory`, `smoke` into a scratch directory | [virtual collisions](cai2012-virtual-collision.md) |
| `cai2012-hybrid` | yes; all recorded interfaces unsuitable | `python -m plumetools.hybridplume setup …` | [hybrid](cai2012-hybrid.md) |
| `markelov1999` | yes, with CaseFoam | `Allmesh --dict-only` in a copy of `baseCase` | [markelov1999](markelov1999.md) |
| `r-4d` | yes | `./Allmesh --describe` | [r-4d](r-4d.md) |
| `3d-inflow` | after obtaining a mesh | `./Allmesh --dict-only` | [3d-inflow](3d-inflow.md) |
| `iss_solar_panels{,_2,_3}` | **no** | read-only inspection | [ISS solar panels](iss-solar-panels.md) |
| `solar_panel_particle_resolution_study` | **no** | read-only inspection | [particle resolution](particle-resolution.md) |

"Runs today" means the supported stack (OpenFOAM v2512, the custom libraries,
Python) can execute the workflow as committed. It never means the result is
production-quality: each tutorial's *Judging the output* section says what to
check, and the [case catalog](../reference/case-catalog.md) lists the evidence
and blockers in one table.

Conventions in the commands: a block marked *WSL / Linux shell* needs bash;
`python script.py …` works with native Windows Python; anything that meshes or
solves needs a sourced OpenFOAM; `pvpython` is ParaView's interpreter. Commands
that delete results are shown only inside warning boxes, never in a first-run
sequence.
