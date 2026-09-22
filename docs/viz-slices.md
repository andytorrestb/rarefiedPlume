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
vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100_np1x
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

### Each study has its own sampling settings

| study | spec |
|---|---|
| `cases/cai2012` | [`cases/cai2012/viz.yaml`](../cases/cai2012/viz.yaml) |
| `cases/markelov1999` | [`cases/markelov1999/viz.yaml`](../cases/markelov1999/viz.yaml) |
| `cases/cai2012-health` | [`cases/cai2012-health/viz.yaml`](../cases/cai2012-health/viz.yaml) |

Each sits beside the `AllpostCases` that uses it, and each is what that study
draws by default. They are **not** shared: the two families have different
geometry, different boundaries and different fields worth looking at.

Both `extends: default`, which inherits
[`plumetools/viz/slices.yaml`](../plumetools/viz/slices.yaml) — image size,
colour preset, output paths, and the `U`/`Ttra` derivations. Mappings merge key
by key; `planes` and `fields` are replaced wholesale, so what a study file lists
is exactly what it draws.

What actually differs:

| | `cai2012` | `markelov1999` |
|---|---|---|
| domain | full 3-D box, `y ∈ [-10 D, +10 D]` | **half** domain, `y ∈ [0, 0.3 m]` |
| `y = 0` is | the middle of the mesh | the symmetry plane, **on the mesh edge** |
| solid bodies | none — two vacuums and an inlet | a cylinder and a plate |
| `q`, `fD` | not asked for; there is no wall to strike | drawn, on the wall patches |
| `*Mean` fields | written from `timeStart` onward | none yet — runs stop before `timeStart` |

Both cut one plane, not the default three: rendering is roughly half a minute an
image, so a four-case study on three planes is an hour bolted onto a step that
otherwise takes seconds. Pass `--viz-spec` (or `--spec` directly) for anything
else on demand.

To change what a study draws, edit its `viz.yaml`. To change it for one run
only, `./AllpostCases --viz-spec mine.yaml`.

## Time series

Both studies render **a frame per written time**, each field in its own folder:

```
Cases/Kn100_np1x/results/viz/
    rhoN/rhoN_0000.png  rhoN_0001.png  ...  rhoN_0004.png
    U/U_0000.png        ...
    Ttra/ ...
    manifest.yaml
```

```yaml
sampling:
  time: all                     # 'latest' for one frame; ~5x quicker
output:
  subdirectory: "{requested}"
  filename: "{requested}_{index}"
```

### The videos

Each folder is encoded into a video beside its own frames, automatically, at the
end of the render:

```
results/viz/rhoN/rhoN_0000.png … rhoN_0004.png
                 rhoN.mp4
```

```yaml
video:
  enabled: true
  framerate: 4
  quality: 18       # x264 CRF, lower is better
```

Encoding needs **ffmpeg but not ParaView**, which is the point:
[`plumetools/viz/video.py`](../plumetools/viz/video.py) is importable and runs
under plain `python`, so changing the frame rate is seconds rather than a
re-render. Fifteen cai2012 videos re-encode in about six seconds; rendering
their frames again takes twenty minutes.

```bash
python plumetools/viz/video.py <viz-dir> --framerate 10
python plumetools/viz/video.py cases/cai2012/Cases/*/results/viz \
    --spec cases/cai2012/viz.yaml        # the study's own settings
```

A missing ffmpeg prints one line and is never fatal — the frames are the output,
the video is a convenience over them.

## Dissecting a sweep

[`plumetools/viz/dissect.py`](../plumetools/viz/dissect.py) sits beside
`video.py` under the same contract — ffmpeg, no ParaView, importable, runnable
under plain `python` — and does two things to frames that already exist.

```bash
python plumetools/viz/dissect.py sheet --matrix matrix.yaml --output sheet.png
python plumetools/viz/dissect.py converge <viz-dir> ... --field dsmcRhoN
```

**The contact sheet** puts one case per grid cell — `scale`/`pad` to a common
size, `drawtext` for the label, `xstack` at explicit pixel offsets. A slice's
pixel width follows its geometry, so cells are not naturally the same size;
offsets are computed in pixels rather than `xstack`'s `w0`/`h0` arithmetic,
which silently overlaps cells when an input is not the size assumed.
`--trim-bottom 110` crops the colour bar off each cell, which is worth doing
only when the range is pinned across the matrix and all the bars are therefore
the same legend. A `.yaml` legend is written beside the image: a PNG records
nothing about where its cells came from, and the sheet is the artefact most
likely to be pasted somewhere on its own.

**The convergence curve** compares successive frames of a running average with
`-lavfi psnr` — once the mean has converged, successive frames stop differing.
Each frame is also compared with the *final* one, because successive
differences go to zero for a stalled series as readily as a converged one.

> **PSNR between two colour-mapped PNGs is not a physical error.** The scale is
> logarithmic and clipped to the pinned range, the 8-bit palette hides any
> change below one colour step (and reports `inf` — recorded as `saturated`,
> not as a number), and the vacuum counts as many pixels as the plume. It is a
> perceptual proxy for *has this stopped changing*. The caveat is attached to
> the written data, not just the report, and `convergence_report` prints a
> warning in place of the physical column when `results/metrics.yaml` is
> missing. Where the image metric and the centreline error disagree, that
> disagreement is the finding.

Both are used by [`cases/cai2012-health`](../cases/cai2012-health/).

Three flags in that command are not optional, and each fails in a way that is
easy to miss:

| | |
|---|---|
| `pad=ceil(iw/2)*2:ceil(ih/2)*2` | `libx264` **refuses** odd dimensions — `width not divisible by 2 (345x800)`. A slice's width follows the geometry, so it is odd about half the time. |
| `-f concat` with an explicit list | `-pattern_type glob` sorts lexicographically, which is right only while the names happen to be zero-padded. The list is built from the manifest's frame numbers. |
| `-pix_fmt yuv420p` | not the default for a PNG input; without it the file plays in ffplay and VLC but not in a browser or QuickTime. |

Three more things about a series are not the same as for a single image.

### Order frames by `{index}`, never `{time}`

Times are formatted with `%g`, so `0.0021211` and `0.00990533` come out
different widths and sort into the wrong order in every shell glob and file
browser — which is the order the encoder would then use. `{index}` is
zero-padded and always sorts. The manifest records the real time of every frame.

### The colour range is pinned across the whole series, automatically

Autoscaling each frame to its own data is right for one survey image and ruinous
for a series: the scale moves underneath the flow, so a plume filling a vacuum
looks like a plume that is not changing at all, and no two frames can be
compared. The renderer therefore scans every time first and colours the whole
series on the union. It costs one slice update per plane per time — one probe
reports every array at once — which is small beside the rendering it corrects.

Pinning still has to be explicit *across cases*: each case gets its own union,
so the Kn = 100 and Kn = 0.01 series use different scales. Set `range:` on the
field to compare them side by side.

### A series wants the instantaneous fields, not the averaged ones

Both study specs set `prefer_mean: false`, which is the opposite of the advice
for a single image. Two reasons:

- `rhoNMean` is a **running average** accumulated from `fieldAverage`'s
  `timeStart`. A series of it shows the average converging, not the flow
  evolving — by construction it barely moves.
- It keeps the series homogeneous. With `prefer_mean` on, early times have no
  `*Mean` field and fall back while later ones do not, so the series would
  silently change quantity part way through.

The cost is real: these frames are one timestep's parcels and the speckle is
shot noise. For a single averaged survey image, set `time: latest` and
`prefer_mean: true`.

### A window, for the series that wants the averaged ones after all

The second reason above — the series changing quantity part way through — is a
property of *when the frames were written*, not of averaging. Restrict the
series to the times where the average exists and it goes away:

```yaml
sampling:
  time: all
  time_min: 5.658e-3      # dsmc.average_start_s: before this there is no *Mean
  prefer_mean: true
```

`time_min`/`time_max` filter the written times **before** `time` selects from
them, so `latest` means the latest frame inside the window. Both default to
unbounded, and `manifest.yaml` records them either way — which frames were left
out is not recoverable from a directory of images.

This is what [`cases/cai2012-health`](../cases/cai2012-health/) uses, and it is
the only reason a study can run a series on `prefer_mean: true` at all. That
study is measuring how the running average converges, so the convergence the
section above calls a defect is exactly its subject. Both settings are right;
they answer different questions.

A window that matches no written time is **not** an error — a case still inside
its transient is a case part way through, and `AllpostCases` must not exit
non-zero for one.

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

### Adding a study, or a one-off

A new study family gets a `viz.yaml` beside its `AllpostCases`, the same way
the two existing ones do:

```yaml
# cases/<study>/viz.yaml
extends: default
planes:
  - name: midplane
    normal: [0, 1, 0]
    camera_up: [0, 0, 1]
    origin: [0, 0, 0]
fields:
  - name: rhoN
    log: true
    range: [1e14, 1e20]     # pinned, so the cases compare side by side
```

Mappings merge key by key; lists (`planes`, `fields`, `derived`) are **replaced**
wholesale. A spec listing three fields draws three fields, not three plus
whatever the parent carried.

For a one-off, put the file anywhere and pass it — nothing has to be installed
or registered:

```bash
./AllpostCases --viz-spec /tmp/mine.yaml
```

Do **not** put a spec inside a generated case directory (`Cases/Kn100_np1x/`,
`Cases/gap06in/p005psi/`). Those are gitignored build products that
`generate_cases.py` deletes and rewrites, so edits there are lost on the next
regeneration.

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

## Seven things that are easy to get wrong

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

### A 3-D view is not sized for you

A `Slice` derives its viewport from the cut's own bounding box, so slice images
are already undistorted. A `Visualization3D` takes an **explicit width**, and
anything fixed there stretches the geometry by whatever the difference is —
`height * 4 / 3` drew the markelov cylinder, which is half again as tall as it
is wide, at 1.33 against a true 0.69. Nearly twice too wide.

The width now comes from the bounding box projected onto the camera's own axes
(`normal × camera_up` and `camera_up`), so one world unit is the same number of
pixels across as down. That is the whole of what "represents the geometry" means
here, and it is worth checking on any new view type: measure the PNG and compare
its ratio against the extents it is supposed to show.

Note that a colour bar adds a row *below* the image. It changes the file's
overall ratio without distorting anything inside it.

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

### `prefer_mean` depends on whether you are drawing one image or a series

For **one image**, leave it on. A single-timestep `rhoN` is one step's worth of
parcels — `DSMCCloud::resetFields` zeroes it every step — so an instantaneous
image is a picture of shot noise that looks entirely plausible. `--instantaneous`
turns it off; on a single image that is useful for judging particle statistics
and for nothing else.

For a **series**, turn it off, as both study specs do — see
[Time series](#a-series-wants-the-instantaneous-fields-not-the-averaged-ones)
above. The averaged field is a running average and barely moves; the
instantaneous one is the only one that evolves.

Before `fieldAverage`'s `timeStart` there are no `*Mean` fields at all. With
`prefer_mean` on, the renderer falls back to the instantaneous field and prints
a `note` saying so, rather than drawing nothing — which is also why a series
left on `prefer_mean` would change quantity part way through.

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
| `video.py` | no | encodes frames to video; needs ffmpeg, runs standalone |
| `render.py` | **yes** | turns a spec into PNGs |
| `render_slices.py` | **yes** | the command-line entry point |

The split is drawn there on purpose: the decisions most likely to be quietly
wrong are all on the side that `tests/unit/test_viz.py` can reach with no
ParaView installed.

## Installation

See [`environment.md`](environment.md). In short: ParaView 5.11–5.13 installed
separately, its internal Python matching the virtual environment's minor version,
then `pip install -e ".[viz]"`.
