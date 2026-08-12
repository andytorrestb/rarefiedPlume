# Slice imagery

Survey pictures of any field the solver wrote, on any cutting plane, driven by a
YAML file. Implemented in [`plumetools/viz/`](../plumetools/viz/) on top of
[VifPara](https://github.com/virtual-vehicle/VifPara), which automates ParaView.

This is **not** the validation path. `plumetools.cai2012.post` extracts the four
outputs that get scored against Cai and nothing else, which §9 of the case
specification asks for deliberately. This produces images to look at. The two
share no code and answer different questions.

## Running it

```bash
vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100
```

**Not `python`.** `paraview.simple` exists only inside ParaView's own
interpreter; the `vifpara` launcher re-runs the script under `pvpython` with this
virtual environment injected via `PV_VENV`. `python render_slices.py` fails at
the import, and says so.

```bash
# one field, one plane
vifpara plumetools/viz/render_slices.py <case> --field rhoN --plane xz

# what does this case actually have?
vifpara plumetools/viz/render_slices.py <case> --list

# which images would be drawn? (no ParaView needed)
vifpara plumetools/viz/render_slices.py <case> --dry-run

# a spec of your own, every written time
vifpara plumetools/viz/render_slices.py <case> --spec mine.yaml --time all
```

Output lands in `<case>/results/viz/`, alongside a `manifest.yaml` recording
which field name, which time and which colour range produced each PNG — a
directory of images otherwise records none of that.

## In the study loops

Both `AllpostCases` scripts render the centre plane of every case at the end of
post-processing:

```bash
cases/cai2012/AllpostCases                  # validation outputs, then imagery
cases/markelov1999/AllpostCases

./AllpostCases --no-viz        # numbers only
./AllpostCases --viz-only      # imagery only, skip everything else
./AllpostCases --viz-spec F    # a different sampling configuration
```

The imagery is a **survey, not a result**. It runs *after* the validation
outputs and the study table are on disk, and a render that fails is reported
without costing anyone the numbers. If `vifpara` is not on `PATH` the step prints
one line saying so and the rest is unaffected — the studies still work on a
machine with no ParaView.

The whole study goes through **one** `vifpara` invocation rather than one per
case. ParaView's startup and the reader's mesh scan cost more than the rendering
does for a small case, and `render_slices.py` therefore takes several case
directories and resets the session between them.

### The mid-plane spec

[`plumetools/viz/midplane.yaml`](../plumetools/viz/midplane.yaml) is what those
loops use: one cut, `y = 0`, for both families. One file covers both because
both put the plume on `+x` and their centre plane at `y = 0`:

| | domain | `y = 0` is |
|---|---|---|
| `cases/cai2012` | full 3-D box, `y ∈ [-10 D, +10 D]` | the middle |
| `cases/markelov1999` | **half** domain, `y ∈ [0, 0.3 m]` | the symmetry plane, on the mesh edge |

They differ in which fields mean anything, and that is settled by measurement
rather than by a second config file: `cai2012` has no solid surface, so `q` and
`fD` are identically zero and get skipped; `markelov1999` has a cylinder and a
plate, the plume strikes both, and the same two entries are drawn.

Why one plane and not the default three: rendering is roughly half a minute an
image, so a four-case study on three planes is an hour bolted onto a step that
otherwise takes seconds. Drop `--spec` to get all three on demand.

## The configuration file

[`plumetools/viz/slices.yaml`](../plumetools/viz/slices.yaml) is the default and
is heavily commented. Its shape:

```yaml
version: 1
extends: default          # optional

sampling: {...}           # which time, averaged or instantaneous
image:    {...}           # pixel size, colour preset, colour bar
output:   {...}           # where the PNGs go and what they are called
planes:   [...]           # the cutting planes
derived:  [...]           # fields computed from other fields
fields:   [...]           # what to draw
```

`fields` × `planes` is a cross product. Unknown keys are rejected everywhere — a
misspelled `prefere_mean` that was silently ignored would produce instantaneous
images that look exactly like averaged ones.

### Adding a field

```yaml
fields:
  - name: rhoM
    log: true
```

Everything unstated comes from [`catalog.py`](../plumetools/viz/catalog.py),
which knows what each `dsmcFoam` field is, its units, and whether it is a volume
field at all. A bare string (`- rhoM`) is a valid entry.

### Adding a case-specific spec

```yaml
# cases/cai2012/Cases/Kn0p01/viz.yaml
extends: default
fields:
  - name: rhoN
    range: [1e14, 1e20]     # pinned, so the three Knudsen cases compare
```

Mappings merge key by key; lists (`planes`, `fields`, `derived`) are **replaced**
wholesale. A spec listing three fields draws three fields, not three plus
whatever the parent carried.

### Derived fields

`dsmcFoam` writes moments, not the quantities anyone wants to look at, so
velocity and temperature are Calculator expressions rather than hard-coded:

```yaml
derived:
  - name: U
    expression: "momentum{mean}/rhoM{mean}"
    component: Magnitude
    units: m/s
  - name: Ttra
    expression: "(2*linearKE{mean}/rhoM{mean} - mag(momentum{mean}/rhoM{mean})^2)/(3*{R})"
    units: K
```

`{mean}` becomes `Mean` or nothing, following `sampling.prefer_mean`. `{R}` is
the specific gas constant `k_B/m` and `{m}` the molecular mass, both read from
the case's `case.yaml`. A spec using `{R}` against a case whose gas is unknown is
an error, not a substituted default: a temperature computed with the wrong gas
constant is wrong by a factor nobody would spot in a picture.

## Six things that are easy to get wrong

### A symmetry plane is the edge of the mesh, and cuts nothing

`cases/markelov1999` models `y >= 0`. Its centre plane, `y = 0`, is therefore a
**boundary face**, and VTK's cutter finds no cell interiors there: the slice
comes back empty, every field reads as "not present on this geometry", and no
image is drawn — with nothing to suggest the plane was the problem.

The renderer detects an origin sitting on a bounding face and moves it a
ten-thousandth of the domain inward, printing a note. At that distance the
picture is the symmetry plane to far better than one cell. Only axis-aligned
normals are adjusted; a tilted plane through a corner could be nudged several
ways and none is obviously right, so it is reported instead of guessed at.

### Surface fields are drawn on walls only

The other patches are the vacuum box enclosing the whole domain. Drawing that as
well puts an opaque slab in front of every body — for `markelov1999` it is the
difference between seeing the heat flux on the cylinder and seeing the outside
of a grey box. Patch types come from `constant/polyMesh/boundary`, since
ParaView's reader does not expose them; a mesh with no wall falls back to every
patch.

## Four more

### The sign of the plane normal mirrors the image

The camera looks **along** `normal` with `camera_up` pointing up, so screen-right
is `normal × camera_up`. For Cai's X–Z figures, flow to the right and Z up needs

```yaml
normal: [0, 1, 0]        # [0,1,0] x [0,0,1] = [1,0,0]
camera_up: [0, 0, 1]
```

`[0, -1, 0]` puts the nozzle on the right instead. The picture is otherwise
perfectly correct, which is why nobody notices.

(VifPara's own `camera_up × normal` sizes the viewport. It is not the screen axis
and it points the other way.)

### `prefer_mean` is on, and should stay on

A single-timestep `rhoN` is one step's worth of parcels — `DSMCCloud::resetFields`
zeroes it every step — so an instantaneous image is a picture of shot noise that
looks entirely plausible. `--instantaneous` turns it off; it is useful for
judging particle statistics and for nothing else.

Before `fieldAverage`'s `timeStart` there are no `*Mean` fields at all. The
renderer falls back to the instantaneous field and prints a `note` line saying
so, rather than drawing nothing.

### `q` and `fD` cannot be sliced

Four fields live entirely in their `boundaryField`; `DSMCCloud` accumulates them
when a molecule *strikes a wall*:

| field | | |
|---|---|---|
| `q` | heat flux into a boundary | W/m² |
| `fD` | force density on a boundary | Pa |
| `boundaryT` | wall temperature | K |
| `boundaryU` | wall velocity | m/s |

Their internal field is `uniform 0`. Cutting a plane through the domain and
colouring by `q` gives a uniformly zero image every time. They are tagged
`SURFACE` in the catalogue and drawn on the boundary patches instead.

`cases/cai2012` has **no wall** — its three patches are two vacuums and an inlet
— so `q` and `fD` are identically zero everywhere, patches included, and there is
nothing to draw. `skip_constant_fields` (on by default) drops them with a message
and records the reason in `manifest.yaml`, rather than writing blank PNGs. The
same spec against an impingement case would draw them.

### Colour ranges are set explicitly, never autoscaled

A DSMC plume has cells that are *exactly* zero — no parcel ever reached them —
and a logarithmic scale has no colour for zero. ParaView's response is to
silently substitute a range starting at 1.0, which for a number density in m⁻³
throws away most of the plume. So the range is read off the cut itself and the
logarithmic floor comes from `log_decades` below the maximum:

```yaml
- name: rhoN
  log: true
  log_decades: 8      # Kn = 100 empties the far field completely
```

Pin `range: [min, max]` when comparing cases — autoscaled images of two Knudsen
numbers use two different scales and cannot be read side by side.

## Implementation notes

### The reader's array list is stale

ParaView's OpenFOAM reader builds `CellArrays` from the **first** time directory.
`fieldAverage` does not start writing until `timeStart`, so no `*Mean` field ever
appears in that property — not after `UpdatePipelineInformation()`, not after
updating at a later time, not for a freshly constructed reader.

The data is fine. Left at its default (everything selected) the reader loads
whatever is actually present at the requested time. Calling
`Case.set_cell_arrays()` with a list built from that stale property would
silently *deselect* every averaged field, so `render.py` never restricts the
selection and probes `GetDataInformation().GetCellDataInformation()` at the
resolved time instead — what arrived, not what was advertised.

### The exit code

On a headless/WSLg display, `pvpython --force-offscreen-rendering` dies at exit
releasing its GL contexts during interpreter finalisation:

```
X Error of failed request:  GLXBadContext
Minor opcode of failed request:  5 (X_GLXMakeCurrent)
```

Every image and the manifest are already written and closed by then, but the
process ends with status 1 — backwards for anything looping over cases and
checking `$?`, which is exactly what `AllpostCases` does.

`render_slices.py` calls `pv.Disconnect()` before exiting, which releases the
contexts in a defined order; measured here, that gives exit 0 with the full log
intact. `PLUMETOOLS_VIZ_RAW_EXIT=1` skips it and lets the crash happen, for when
you are debugging ParaView rather than using it.

> The obvious shortcut — `os._exit` to skip teardown entirely — is **wrong**.
> VifPara imports colorama, which replaces `sys.stdout` with a `StreamWrapper`;
> flushing that does not push the underlying buffer out, so `os._exit` silently
> discards every line the run printed. Both behaviours were measured, not
> assumed.

### The layering

| module | ParaView? | |
|---|---|---|
| `spec.py` | no | schema, YAML, validation |
| `catalog.py` | no | what the `dsmcFoam` fields are |
| `resolve.py` | no | which time, which field, what colour scale |
| `render.py` | **yes** | turns a spec into PNGs |
| `render_slices.py` | **yes** | the command-line entry point |

The split is drawn there on purpose: the decisions most likely to be quietly
wrong are all on the side that `tests/unit/test_viz.py` can reach with no
ParaView installed.

## Installation

See [`environment.md`](environment.md). In short: ParaView 5.11–5.13 installed
separately, its internal Python matching the virtual environment's minor version,
then `pip install -e ".[viz]"`.
