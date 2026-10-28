# plumetools: who owns what, and how data flows

How the package is divided, which module owns each responsibility, and where a
new diagnostic, mesh or study variation goes. Written 2026-10-27 with the
maintainability refactor; [`refactor-2026-10-27.md`](refactor-2026-10-27.md) is
the change record and migration note.

## Layers

```
case scripts        All*, runCase.sh, generate_*.py         thin adapters: arguments,
  (cases/)                                                   logging, OpenFOAM commands
        |
families            cai2012  cai2012ext  markelov1999        the science: configuration,
  (plumetools/)     nozzleplume  virtualcollision            derivation, assessment,
                    sourceflow/inflow/config (3d-inflow)     dictionaries, post-processing
        |
generic tools       rarefaction  viz  iteration  casemanifest
        |
shared OpenFOAM     foamio/  errors times volfield patchfields checkpoints
  I/O                        primitives blockmesh snappy fields
                    mesh/    boundary polymesh patches sets
```

Each layer imports only from the layers below it. A family may import another
family's **public** model (the Kn extension builds on `cai2012`, the interface
closure calls `markelov1999.sourceflow`), never its private helpers; generic
tools import no family's post-processing. A unit test pins that the extracted
interface modules do not import their coordinator.

## The shared OpenFOAM layer

| module | owns |
|---|---|
| `foamio/errors.py` | `FoamDataError`, the one error the readers raise (a `ValueError`) |
| `foamio/times.py` | which subdirectories are time directories; `latest_time` |
| `foamio/volfield.py` | reading a field's `internalField` or one patch's `value`, with its metadata (`FieldValues`: uniform or list, value type, declared count); `FieldValues.over(n, what)` is the only place a value is broadcast |
| `foamio/patchfields.py` | the patchField type each patch type takes; the `fixedValue` list entry |
| `foamio/checkpoints.py` | restart checkpoints across ranks: complete or not, what a resume must do |
| `foamio/primitives.py`, `blockmesh.py`, `snappy.py` | dictionary rendering shared by the meshing families |
| `foamio/fields.py` | the legacy 3d-inflow `0/` writer (frozen by the regression golden) |
| `mesh/boundary.py`, `sets.py` | `polyMesh/boundary`, face sets |
| `mesh/polymesh.py` | points, faces (whole, by range, by label -- streamed), owner/neighbour, cell count |
| `mesh/patches.py` | `PatchGeometry`: centres, area vectors, normals of boundary patches |

Rules the readers keep, so callers do not have to:

- A list must hold the count it declares and, through `over`, one value per
  cell or face; only a `uniform` entry is expanded. A binary file, a
  non-number, or a list from another mesh is an error that says which.
- Readers stream and close: a `faces` file is passed over once per call and
  only the requested faces are kept; nothing is cached between calls, so a
  remeshed case is never read through a stale copy.
- Time directory names are kept as OpenFOAM wrote them; comparison is numeric,
  identity is by name.

## Data flow, per family

**Direct-exit Cai extension case** (`cai2012ext`):

```
case.yaml --direct_case.derive--> DirectCase(cfg, geom, exit_state, plan, run,
                                             nozzle, particles, report)
          --feasibility.assess_case--> verdict dict
          --> generation.build_entry (manifest row)  /  hpc.label (NON-PRODUCTION)
          --> interface_case.post_context (Allcompare, Allconverge)
```

**Interface-mode case** (`cai2012ext`, an enabled `interface.yaml`):

```
interface.yaml --interface_config.load_interface--> mode (the one place it is decided)
case.yaml + interface.yaml --interface_case.derive--> Derived
        Derived.planned: InterfacePlan | None      all planned quantities, or none
        --interface_assessment.assess--> verdict   (writes nothing)
              uses interface_evidence.coupling_validation   overlap.yaml contract
        --interface_report.report / manifest_entry / describe--> files and lines
        --interface_identity--> interface-identity.yaml       restart compatibility
        --interface_case.write_inflow--> 0/ (foamio.patchfields, read back by volfield)
```

`interface.py` is the pure model (equivalent source, Bird's parameter, sphere
selection, inlet faces); `interface_mesh.py` plans and renders the mesh.

**Nozzle-plume case** (`nozzleplume`): `workflow.load(case_dir)` gives a typed
`Case`; `setup`, `fields`, `post` and `diagnostics` consume it.

**Running** (every family's `Allrun`, and the PBS steps):

```
Allrun --python -m plumetools.foamio.checkpoints plan . --np N-->
        RESUME / DECOMPOSE (none|latest|all) / RECONSTRUCT_FIRST, or exit 3 with
        what is incomplete;  Allrun then runs the OpenFOAM tools accordingly
hpc finished, health AllrunCases --checkpoints.reached_end-->
        the newest checkpoint complete on every rank is at endTime
family runners (AllmeshCases, ...) --python -m plumetools.casemanifest-->
        every case path, after the whole manifest is checked
```

## Extending it

**A diagnostic.** One registered function in `rarefaction/core.py`; a new
input is read in `rarefaction/foam.py` with the shared readers
([`rarefaction-diagnostics.md`](rarefaction-diagnostics.md#adding-a-diagnostic)).
A family-specific post-processing quantity goes in that family's `post.py`,
reading through `foamio.volfield` and `mesh.patches`, and keeps that family's
error class by translating `FoamDataError` (as `cai2012.post.read_cell_values`
does).

**Mesh support.** A family's mesh module has three parts, and a new mesh type
should keep them: `plan(cfg, geom, exit_state) -> MeshPlan` (a frozen
dataclass, nothing written), `render_*(cfg, geom, plan) -> str` per dictionary
(no planning recomputed), and `write_mesh_setup` (the only writer). Reading a
built mesh uses `mesh/`; nothing new parses `polyMesh` itself.

**A study variation.** A new Knudsen number or parcel level is data
(`study.yaml`, `interfaces.yaml`). A new derived quantity belongs in the
derivation (`direct_case.derive_config`, `interface_case.derive`), a new
criterion in the assessment, a new column in the report -- each in one place,
so generation, the PBS steps and post-processing see the same number.

## What is deliberately not shared

- The paper families' models stay separate: `cai2012`, `markelov1999` and
  `nozzleplume` each keep their own gas, exit-state and dictionary writers.
  Their outputs are pinned byte-for-byte; merging them would make one
  family's change silently another's.
- Independent checks stay independent: `cai2012.audit`, the verification
  modules and the tests' oracles compute what they check without calling the
  code under check.
- `interface_case.openfoam_normals` and the legacy `plumetools.geometry`
  centroids are not replaced by `mesh.patches.face_geometry`: different
  summation order could change written fields in the last bit.
