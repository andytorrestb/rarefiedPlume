r"""Turn a :class:`~plumetools.viz.spec.VizSpec` into PNGs, through VifPara.

**This module imports ParaView.** It can only be imported by ``pvpython``, which
in practice means it can only be imported by a script launched with the
``vifpara`` command. Everything importable without ParaView lives in
:mod:`plumetools.viz.spec` and :mod:`plumetools.viz.catalog`; see
:mod:`plumetools.viz` for why the split is drawn there.

Three things here are not obvious from VifPara's API, and all three are the
difference between a correct image and a plausible one.

1. The reader's array list is stale, so never restrict it
---------------------------------------------------------
ParaView's OpenFOAM reader builds ``CellArrays`` from the **first** time
directory. ``fieldAverage`` does not start writing until ``timeStart``, so at
the first time directory there are no ``*Mean`` fields at all, and
``reader.CellArrays`` therefore never lists them -- not after
``UpdatePipelineInformation``, not after updating at a later time, not for a
freshly constructed reader.

The data is fine. Left at its default (everything selected) the reader loads
whatever is actually present at the requested time, ``rhoNMean`` included.
Calling ``Case.set_cell_arrays`` with a list built from that stale property
would silently *deselect* every averaged field. So this module never restricts
the selection, and probes availability from
``GetDataInformation().GetCellDataInformation()`` at the resolved time instead
-- what arrived, not what was advertised.

2. Colour ranges are set explicitly, never autoscaled
------------------------------------------------------
``RescaleTransferFunctionToDataRange`` would autoscale to the whole reader
output rather than to the cut, and a logarithmic scale would then be handed a
minimum of exactly zero -- a DSMC plume has cells no parcel ever reached. So the
range is read off the slice itself, the logarithmic floor is set from
``log_decades``, and :class:`~vifpara.ColorMap` is given both ends.

3. VifPara slices a ``Case``; derived fields are not a ``Case``
---------------------------------------------------------------
``Slice(case=...)`` feeds ``case.get_case()`` into ``pv.Slice`` and touches
nothing else on it, so :class:`_SourceAsCase` presents a Calculator's output
under the same one-method interface. That is what lets ``U`` and ``Ttra`` --
which ``dsmcFoam`` never writes -- be sliced like any other field.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from pathlib import Path
from typing import Any, Iterable, Sequence

import paraview.simple as pv
import yaml
from vifpara import (
    Case,
    CaseType,
    ColorMap,
    Config,
    Exporter,
    Layout,
    Slice,
    Vector3,
    Visualization3D,
    logger,
)

from plumetools.viz.catalog import FieldKind
from plumetools.viz.resolve import (
    NothingToRender,
    RenderError,
    colour_range,
    component_range,
    decade_labels,
    expression_fields,
    is_constant,
    nudge_plane_origin,
    resolve_field_name,
    resolve_times,
    scalar_bar_label_budget,
    view_shape,
)
from plumetools.viz.spec import PlaneSpec, RenderTask, VizSpec, VizSpecError

__all__ = [
    "CaseFacts",
    "NothingToRender",
    "RenderError",
    "RenderedImage",
    "Renderer",
    "SkippedImage",
    "array_ranges",
    "colour_range",
    "component_range",
    "decade_labels",
    "discover_case",
    "render_case",
    "resolve_field_name",
    "resolve_times",
]


# --------------------------------------------------------------------------- #
# the case
# --------------------------------------------------------------------------- #

@dataclass
class CaseFacts:
    """What the renderer needs to know about a case beyond its file path.

    Attributes:
        directory: the case directory.
        name: its basename, for ``{case}`` and the manifest.
        foam_file: the ``.foam`` stub ParaView opens.
        length_scale_m: the ``L`` of ``origin_over_L``. The nozzle diameter for
            a Cai case; ``None`` when the case is of a family this cannot read,
            in which case a spec using ``origin_over_L`` fails rather than
            guessing.
        substitutions: ``{R}`` and ``{m}`` for derived-field expressions.
    """

    directory: Path
    name: str
    foam_file: Path
    length_scale_m: float | None = None
    substitutions: dict[str, str] = dataclass_field(default_factory=dict)


def _cai2012_facts(directory: Path) -> tuple[float, float] | None:
    """``(length_scale_m, mass_kg)`` if this is a Cai case, else ``None``.

    The Cai family's characteristic length is the nozzle **diameter** -- the
    same convention its Knudsen number uses, so ``origin_over_L: [5, 0, 0]``
    means the same "5 D" the paper's axes do.
    """
    from plumetools.cai2012.config import load_case_config

    cfg = load_case_config(directory)
    return float(cfg.nozzle.diameter_m), float(cfg.gas.mass_kg)


def _source_flow_facts(directory: Path) -> tuple[float, float] | None:
    """``(length_scale_m, mass_kg)`` if this is a source-flow case, else ``None``.

    The Markelov family has no nozzle in the domain -- the physical orifice is
    0.83 mm and sits outside it. Its length scale is the **inflow hemisphere
    radius**, which is the only geometric scale a plane origin could sensibly be
    quoted in.
    """
    from plumetools.config import load_case_config
    from plumetools.markelov1999.constants import molecular_mass_kg

    cfg = load_case_config(directory)
    mass_kg = cfg.dsmc.species.mass_kg
    if mass_kg is None:
        mass_kg = molecular_mass_kg(cfg.dsmc.species.molar_mass_g_per_mol)
    return float(cfg.geometry.sphere_radius_m), float(mass_kg)


#: Case-family readers, tried in order. Each raises if the case is not its own.
CASE_READERS = (_cai2012_facts, _source_flow_facts)


def discover_case(case_dir: str | Path) -> CaseFacts:
    """Locate the ``.foam`` stub and read what the case knows about itself.

    The gas and the length scale come from the case's own ``case.yaml``, via
    whichever family loader recognises it -- the two schemas reject each other
    on their ``model:`` key, so trying both in turn is unambiguous.

    Both are optional. A case of a family neither loader knows still renders; it
    just cannot use ``origin_over_L`` or the ``{R}`` and ``{m}`` placeholders,
    and the specs that need those fail with their own message rather than
    silently substituting a wrong gas constant.

    Args:
        case_dir: an OpenFOAM case directory.

    Returns:
        The assembled :class:`CaseFacts`.

    Raises:
        RenderError: if the directory or the ``.foam`` stub is missing.
    """
    directory = Path(case_dir).resolve()
    if not directory.is_dir():
        raise RenderError(f"no such case directory: {directory}")

    stubs = sorted(directory.glob("*.foam"))
    if not stubs:
        raise RenderError(
            f"{directory} has no *.foam stub for ParaView to open.\n"
            f"       Create one with:  touch {directory / 'open.foam'}")
    facts = CaseFacts(directory=directory, name=directory.name,
                      foam_file=stubs[0])

    from plumetools.cai2012.constants import BOLTZMANN_J_PER_K

    for reader in CASE_READERS:
        try:
            resolved = reader(directory)
        except Exception:
            # Not this family's case, or an unreadable one. Try the next.
            continue
        if resolved is None:                              # pragma: no cover
            continue
        length_scale_m, mass_kg = resolved
        facts.length_scale_m = length_scale_m
        facts.substitutions = {
            "R": repr(BOLTZMANN_J_PER_K / mass_kg),
            "m": repr(mass_kg),
        }
        break
    return facts


class _SourceAsCase:
    """Presents a ParaView source where VifPara expects a :class:`vifpara.Case`.

    ``Slice`` and ``Visualization3D`` call exactly one method on the object
    handed to them -- ``get_case()`` -- and use its return value as a pipeline
    input. Wrapping a Calculator's output in this adapter is what lets a derived
    field be sliced; see the module docstring.
    """

    def __init__(self, source: Any) -> None:
        self._source = source

    def get_case(self) -> Any:
        return self._source


# --------------------------------------------------------------------------- #
# time
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# arrays
# --------------------------------------------------------------------------- #

def array_ranges(source: Any, time: float,
                 association: str) -> dict[str, tuple[float, ...]]:
    """Component ranges of every array actually present on *source* at *time*.

    Args:
        source: any ParaView pipeline object.
        time: the time to update at.
        association: ``"CELLS"`` or ``"POINTS"``.

    Returns:
        ``{name: (min_c0, max_c0, ..., min_mag, max_mag)}`` -- flattened per
        component with the magnitude last, indexed by
        :func:`component_range`.
    """
    source.UpdatePipeline(time)
    information = source.GetDataInformation()
    dataset = (information.GetCellDataInformation() if association == "CELLS"
               else information.GetPointDataInformation())

    ranges: dict[str, tuple[float, ...]] = {}
    for index in range(dataset.GetNumberOfArrays()):
        array = dataset.GetArrayInformation(index)
        components = array.GetNumberOfComponents()
        flattened: list[float] = []
        for component in range(components):
            flattened.extend(array.GetComponentRange(component))
        flattened.extend(array.GetComponentRange(-1))     # magnitude, or scalar
        ranges[array.GetName()] = tuple(flattened)
    return ranges


# --------------------------------------------------------------------------- #
# rendering
# --------------------------------------------------------------------------- #

@dataclass
class RenderedImage:
    """One PNG that was written, for the manifest."""

    path: Path
    field: str
    requested_field: str
    plane: str
    time: float
    value_min: float
    value_max: float
    log: bool
    component: str
    kind: str


@dataclass
class SkippedImage:
    """One PNG that was not written, and why."""

    field: str
    plane: str
    time: float
    reason: str


class Renderer:
    """Draws one case, over however many times and fields the spec asks for.

    The reader is opened once and reused; only the derived-field Calculators are
    rebuilt per time, because which fields exist changes with the time (before
    ``fieldAverage``'s ``timeStart`` there are no averages to compute from).
    """

    def __init__(self, spec: VizSpec, facts: CaseFacts, *,
                 allow_constant: bool = False, verbose: bool = False) -> None:
        self.spec = spec
        self.facts = facts
        self.allow_constant = allow_constant
        self.verbose = verbose

        self.length_scale_m = (spec.sampling.length_scale_m
                               or facts.length_scale_m)
        self.output_dir = self._resolve_output_dir()
        log_dir = Path(spec.output.log_directory)
        self.log_dir = (log_dir if log_dir.is_absolute()
                        else self.output_dir / log_dir)

        self.written: list[RenderedImage] = []
        self.skipped: list[SkippedImage] = []
        self._origins: dict[str, tuple[float, float, float]] = {}
        self._surface_regions: list[str] | None = None

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.config = Config(custom_config={
            "case_path": str(facts.foam_file),
            "plot_path": str(self.output_dir),
            "log_path": str(self.log_dir),
        })
        logger.set_log_path(self.config.get_log_path())

        self.case = Case(config=self.config, case_type=CaseType.RECONSTRUCTED)
        self.reader = self.case.get_case()

        # internalMesh only, for two reasons: the lagrangian block carries its
        # own 'U' and 'origId' arrays that would collide by name with a derived
        # field, and slicing a point cloud produces nothing anyway.
        self.case.set_mesh_regions(list(spec.sampling.mesh_regions))

    # -- paths -------------------------------------------------------------- #

    def _resolve_output_dir(self) -> Path:
        directory = Path(self.spec.output.directory)
        if not directory.is_absolute():
            directory = self.facts.directory / directory
        return directory

    def _filename(self, task: RenderTask, resolved_field: str, time: float,
                  many_times: bool) -> str:
        name = self.spec.output.filename.format(
            field=resolved_field,
            plane=task.plane_name,
            case=self.facts.name,
            time=format(time, "g"),
        )
        # Rendering several times into one template that does not mention the
        # time would write each image over the last.
        if many_times and "{time}" not in self.spec.output.filename:
            name = f"{name}_t{format(time, 'g')}"
        return name

    # -- the pipeline ------------------------------------------------------- #

    def _build_derived(self, time: float, available: Iterable[str]) -> Any:
        """Chain a Calculator per derived field onto the reader.

        Returns the tail of the chain, which is the reader itself when the spec
        defines no derived fields or none of them can be computed.
        """
        source = self.reader
        present = set(available)
        prefer_mean = self.spec.sampling.prefer_mean

        for entry in self.spec.derived:
            # {mean} has to follow what the case WROTE, not what the spec asked
            # for. A run that stopped before fieldAverage's timeStart has no
            # *Mean fields at all, and an expression built from the setting
            # alone would name arrays that do not exist -- so try the preferred
            # expansion, then the other one.
            expression = None
            missing: list[str] = []
            for use_mean in ((True, False) if prefer_mean else (False, True)):
                try:
                    candidate = entry.resolved_expression(
                        use_mean, self.facts.substitutions)
                except VizSpecError as exc:
                    self._skip(entry.name, "-", time, str(exc))
                    expression = None
                    missing = []
                    break
                # An expression naming a field this time has no value for
                # produces a Calculator that fails at render with a ParaView
                # dialog's worth of C++ trace. Check first, say something useful.
                absent = sorted(expression_fields(candidate) - present)
                if not absent:
                    expression = candidate
                    if use_mean is not prefer_mean:
                        print(f"    note  {entry.name}: no averaged fields at "
                              f"this time; derived from the instantaneous ones, "
                              f"which are one timestep's parcels")
                    break
                missing = missing or absent
            if expression is None:
                if missing:
                    self._skip(entry.name, "-", time,
                               f"expression needs {missing}, which this time "
                               f"does not have")
                continue

            calculator = pv.Calculator(Input=source)
            calculator.AttributeType = (
                "Cell Data" if self.spec.sampling.field_type == "CELLS"
                else "Point Data")
            calculator.ResultArrayName = entry.name
            calculator.Function = expression
            # Empty cells have rhoM = 0, so every expression that divides by a
            # density hits 0/0 out in the vacuum. Replace rather than propagate
            # a NaN, which would poison the colour range for the whole image.
            calculator.ReplaceInvalidResults = 1
            calculator.ReplacementValue = 0.0
            source = calculator

        return source

    def _skip(self, field: str, plane: str, time: float, reason: str) -> None:
        self.skipped.append(SkippedImage(field=field, plane=plane, time=time,
                                         reason=reason))
        print(f"    skip  {field:<20} {plane:<10} {reason}")

    # -- the run ------------------------------------------------------------ #

    def run(self, only_fields: Sequence[str] | None = None,
            only_planes: Sequence[str] | None = None) -> list[RenderedImage]:
        """Render everything the spec asks for.

        Args:
            only_fields: restrict to these field names.
            only_planes: restrict to these plane names.

        Returns:
            The images written, in spec order.
        """
        tasks = self.spec.tasks(only_fields=only_fields, only_planes=only_planes)
        if not tasks:
            raise RenderError("the spec produced no images to draw "
                              "(every field disabled, or filtered out)")

        times = resolve_times(self.reader.TimestepValues, self.spec.sampling.time)
        scene = pv.GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()

        print(f"### {self.facts.name}")
        print(f"    case        {self.facts.foam_file}")
        print(f"    output      {self.output_dir}")
        print(f"    times       {', '.join(format(t, 'g') for t in times)}")
        print(f"    tasks       {len(tasks)}")
        print()

        for time in times:
            scene.AnimationTime = time
            self._render_time(tasks, time, many_times=len(times) > 1)

        if self.spec.output.manifest:
            self._write_manifest(times)
        return self.written

    def _render_time(self, tasks: Sequence[RenderTask], time: float,
                     many_times: bool) -> None:
        association = self.spec.sampling.field_type

        available = array_ranges(self.reader, time, association)
        if not available:
            # resolve_times cannot catch this. ParaView's OpenFOAM reader
            # reports a timestep even for a case with no time directories at
            # all, so a meshed-but-never-run case looks like a case at t = 0
            # with every field missing -- which would otherwise print one
            # baffling "not written at this time" per field and an empty
            # manifest.
            raise NothingToRender(
                f"no {association.lower()} field data at t = {time:g}. This "
                f"case looks meshed but never run -- there is nothing to draw.\n"
                f"       Run ./Allrun in {self.facts.directory} first.")

        source = self._build_derived(time, available)
        if source is not self.reader:
            available = array_ranges(source, time, association)

        if self.verbose:
            print(f"    arrays at t = {format(time, 'g')}: "
                  f"{', '.join(sorted(available))}")
            print()

        derived = self.spec.derived_by_name
        volume_source = _SourceAsCase(source)

        # One probe slice per plane, reused for every field on it: the ranges of
        # all arrays come out of the same GetDataInformation call.
        probes: dict[str, dict[str, tuple[float, ...]]] = {}

        for task in tasks:
            field = task.field
            prefer_mean = self.spec.prefer_mean_for(field)
            candidates = field.candidate_names(prefer_mean, derived)
            resolved = resolve_field_name(candidates, available)

            if resolved is None:
                self._skip(field.name, task.plane_name, time,
                           f"not written at this time (looked for "
                           f"{', '.join(candidates)})")
                continue
            if resolved != candidates[0] and prefer_mean:
                print(f"    note  {field.name}: no {candidates[0]} at this time; "
                      f"using the instantaneous {resolved}, which is one "
                      f"timestep's parcels")

            kind = field.resolve_kind(derived)
            plane = task.plane
            if plane is None:
                self._skip(field.name, task.plane_name, time,
                           "no plane to view it from")
                continue

            # -- the value range --------------------------------------------- #
            if kind is FieldKind.VOLUME:
                key = plane.name
                if key not in probes:
                    probes[key] = self._probe_slice(source, plane, time,
                                                    association)
                geometry_ranges = probes[key]
            else:
                geometry_ranges = self._probe_surface(time, association)

            if resolved not in geometry_ranges:
                self._skip(field.name, task.plane_name, time,
                           f"{resolved} is not present on this geometry")
                continue

            component = field.resolve_component(derived)
            low, high = component_range(geometry_ranges[resolved], component)
            log = field.resolve_log(derived)

            if field.range is not None:
                scale = field.range
            else:
                scale = colour_range(low, high, log, field.log_decades)
                if scale is None:
                    reason = (f"constant at {low:.4g} over this geometry"
                              if is_constant(low, high)
                              else f"no positive values for a log scale "
                                   f"(max {high:.4g})")
                    if not self.allow_constant:
                        self._skip(field.name, task.plane_name, time, reason)
                        continue
                    scale = (low, high if high > low else low + 1.0)

            # -- draw it ------------------------------------------------------ #
            filename = self._filename(task, resolved, time, many_times)
            colour_map = self._colour_map(field, resolved, component, log, scale,
                                          derived)
            show_bar = (field.color_bar if field.color_bar is not None
                        else self.spec.image.color_bar)
            if kind is FieldKind.VOLUME:
                self._render_slice(volume_source, plane, colour_map, filename,
                                   show_bar, log, scale)
            else:
                self._render_surface(plane, colour_map, filename, show_bar,
                                     log, scale)

            path = self.output_dir / f"{filename}.png"
            self.written.append(RenderedImage(
                path=path, field=resolved, requested_field=field.name,
                plane=task.plane_name, time=time,
                value_min=scale[0], value_max=scale[1], log=log,
                component=component, kind=kind.value))
            print(f"    wrote {path.name:<44} "
                  f"[{scale[0]:.3e}, {scale[1]:.3e}]{'  log' if log else ''}")

        for probe in getattr(self, "_probe_objects", []):
            pv.Delete(probe)
        self._probe_objects = []

    # -- geometry probes ----------------------------------------------------- #

    def plane_origin(self, plane: PlaneSpec) -> tuple[float, float, float]:
        """Where this plane actually cuts, once moved clear of a bounding face.

        A half domain's symmetry plane sits on the edge of the mesh, where the
        cut comes back empty; see :func:`~plumetools.viz.resolve.nudge_plane_origin`.
        Resolved once and cached, so the probe and the drawn slice cannot end up
        on different planes.
        """
        if plane.name in self._origins:
            return self._origins[plane.name]

        requested = plane.origin_m(self.length_scale_m)
        bounds = self.reader.GetDataInformation().GetBounds()
        origin, moved = nudge_plane_origin(requested, plane.normal, bounds)
        if moved:
            print(f"    note  plane '{plane.name}' lies on the domain boundary "
                  f"(a symmetry plane); moved it "
                  f"{max(abs(a - b) for a, b in zip(origin, requested)):.3g} m "
                  f"inside, or the cut would be empty")
        self._origins[plane.name] = origin
        return origin

    def _probe_slice(self, source: Any, plane: PlaneSpec, time: float,
                     association: str) -> dict[str, tuple[float, ...]]:
        """Ranges of every array on this cut, so the colour scale fits the image.

        The whole-domain range would be wrong for any plane that misses the
        plume's peak, and autoscaling inside ParaView uses exactly that.
        """
        probe = pv.Slice(Input=source)
        probe.SliceType = "Plane"
        probe.SliceType.Origin = list(self.plane_origin(plane))
        probe.SliceType.Normal = list(plane.normal)
        ranges = array_ranges(probe, time, association)

        if not ranges:
            print(f"    note  plane '{plane.name}' cuts nothing: the origin "
                  f"{list(self.plane_origin(plane))} is outside the mesh, whose "
                  f"bounds are "
                  f"{[round(b, 4) for b in self.reader.GetDataInformation().GetBounds()]}")

        if not hasattr(self, "_probe_objects"):
            self._probe_objects: list[Any] = []
        self._probe_objects.append(probe)
        return ranges

    def _probe_surface(self, time: float, association: str
                       ) -> dict[str, tuple[float, ...]]:
        """Ranges over the boundary patches, for a surface field."""
        self.case.set_mesh_regions(self._patch_regions())
        try:
            return array_ranges(self.reader, time, association)
        finally:
            self.case.set_mesh_regions(list(self.spec.sampling.mesh_regions))

    def _patch_regions(self) -> list[str]:
        """The reader blocks to draw for a surface field.

        **Walls only**, when the mesh has any. A surface field is accumulated
        where molecules *strike something solid*, and the domain's other patches
        are the vacuum box enclosing the whole case -- drawing that too puts an
        opaque slab in front of every body in the picture. For
        ``cases/markelov1999`` this is the difference between seeing the heat
        flux on the cylinder and seeing the outside of a grey box.

        Patch types come from ``constant/polyMesh/boundary``; ParaView's reader
        does not expose them. A mesh with no wall at all falls back to every
        patch, which is the only thing left to draw and is what a case with an
        interesting *inflow* boundary field would want anyway.
        """
        if self._surface_regions is not None:
            return self._surface_regions

        blocks = [entry["name"] for entry in self.case.get_patch_array_info()
                  if entry["name"].startswith("patch/")]
        walls: list[str] = []
        try:
            from plumetools.mesh.boundary import read_boundary
            patches = read_boundary(self.facts.directory)
            wall_names = {name for name, info in patches.items()
                          if info.type.lower() == "wall"}
            walls = [block for block in blocks
                     if block.split("/", 1)[-1] in wall_names]
        except Exception as exc:
            print(f"    note  could not read constant/polyMesh/boundary "
                  f"({exc}); drawing every patch")

        if not walls:
            print("    note  this mesh has no wall patch; drawing every "
                  "boundary for the surface fields")
        self._surface_regions = walls or blocks
        return self._surface_regions

    # -- drawing ------------------------------------------------------------- #

    def _colour_map(self, field: Any, resolved: str, component: str, log: bool,
                    scale: tuple[float, float], derived: dict) -> ColorMap:
        image = self.spec.image
        return ColorMap(
            field=resolved,
            legend_title=field.resolve_legend_title(resolved, derived),
            component_title=component,
            preset=field.preset or image.preset,
            min_value=scale[0],
            max_value=scale[1],
            use_log_scale=log,
            field_type=field.field_type or self.spec.sampling.field_type,
            orientation=image.orientation,
            location=image.location,
            size=image.bar_size,
            legend_format_type=image.legend_format_type,
            legend_format_digits_after=image.legend_digits_after,
        )

    def _render_slice(self, source: _SourceAsCase, plane: PlaneSpec,
                      colour_map: ColorMap, filename: str, show_bar: bool,
                      log: bool, scale: tuple[float, float]) -> None:
        """Cut the plane, colour it, and save it.

        A ``Slice`` never draws its own colour bar -- it applies the colour map
        with the legend explicitly hidden and offers a *separate*
        ``ColorBarView`` instead. So a slice with a legend is a two-cell layout:
        the image on top, the bar underneath.
        """
        image = self.spec.image
        layout = Layout([1, 1] if show_bar else [1])
        view = Slice(
            case=source,
            origin=Vector3.from_list(list(self.plane_origin(plane))),
            normal=Vector3.from_list(list(plane.normal)),
            camera_up=Vector3.from_list(list(plane.camera_up)),
            height=plane.height if plane.height is not None else image.height,
            color_map=colour_map,
            zoom=plane.zoom if plane.zoom is not None else image.zoom,
            margin_x=plane.margin_x if plane.margin_x is not None else image.margin_x,
            margin_y=plane.margin_y if plane.margin_y is not None else image.margin_y,
            offset_x=plane.offset_x if plane.offset_x is not None else image.offset_x,
            offset_y=plane.offset_y if plane.offset_y is not None else image.offset_y,
            show_orientation_axis=image.show_orientation_axis,
        )
        view.render(layout, row=0, col=0)
        extra: list[Any] = []
        if show_bar:
            view.set_color_bar_size(height=self.spec.image.bar_height)
            view.render_color_bar(layout, row=1, col=0)
            colour_bar = view.get_color_bar()
            extra.append(colour_bar)
            self._tune_scalar_bar(
                colour_bar.get_render_view(), log, scale,
                bar_width_px=colour_bar.get_width() * image.bar_size)
        self._save(view, layout, filename, extra_views=extra)

    def _tune_scalar_bar(self, render_view: Any, log: bool,
                         scale: tuple[float, float],
                         bar_width_px: float = 0.0) -> None:
        """Stop a wide range's tick labels from printing on top of each other.

        ParaView labels a logarithmic scale at whatever interval it likes, which
        over the eight decades of a plume's density puts a dozen labels into a
        bar a few hundred pixels wide -- they overlap into an unreadable smear.
        Decade ticks, thinned to at most :data:`MAX_SCALAR_BAR_LABELS`, are
        legible and are what a reader of a log plot expects anyway.

        Every property is set only if this ParaView build has it: the scalar bar
        widget's property set has changed across versions, and a missing one is
        a cosmetic loss, not a reason to fail the render.
        """
        image = self.spec.image
        low, high = scale
        label_format = f"%-#6.{image.legend_digits_after}{image.legend_format_type}"

        for representation in render_view.Representations:
            if representation.GetXMLName() != "ScalarBarWidgetRepresentation":
                continue
            properties = set(representation.ListProperties())

            def put(name: str, value: Any) -> None:
                if name in properties:
                    try:
                        setattr(representation, name, value)
                    except (AttributeError, ValueError):    # pragma: no cover
                        pass

            put("AutomaticLabelFormat", 0)
            put("LabelFormat", label_format)
            put("RangeLabelFormat", label_format)
            put("AddRangeLabels", 1)

            if log:
                budget = scalar_bar_label_budget(bar_width_px)
                labels = decade_labels(float(low), float(high), budget) if budget else []
                put("UseCustomLabels", 1)
                put("CustomLabels", labels)

    def _render_surface(self, plane: PlaneSpec, colour_map: ColorMap,
                        filename: str, show_bar: bool, log: bool,
                        scale: tuple[float, float]) -> None:
        """Draw a boundary-only field on the patches, viewed along *plane*."""
        image = self.spec.image
        self.case.set_mesh_regions(self._patch_regions())
        try:
            bounds = self.reader.GetDataInformation().GetBounds()
            centre = Vector3(
                0.5 * (bounds[0] + bounds[1]),
                0.5 * (bounds[2] + bounds[3]),
                0.5 * (bounds[4] + bounds[5]),
            )
            span = max(bounds[1] - bounds[0], bounds[3] - bounds[2],
                       bounds[5] - bounds[4]) or 1.0
            normal = Vector3.from_list(list(plane.normal)).normalized()
            camera_up = Vector3.from_list(list(plane.camera_up))
            height = plane.height if plane.height is not None else image.height
            zoom = plane.zoom if plane.zoom is not None else image.zoom
            # The viewport is sized from the geometry, not fixed at 4:3, so one
            # world unit is the same number of pixels across as down.
            width, half_height = view_shape(bounds, plane.normal,
                                            plane.camera_up, height)

            layout = Layout([1])
            view = Visualization3D(
                case=self.case,
                cam_position=centre - normal * span,
                cam_up=camera_up,
                focal_point=centre,
                width=width,
                height=height,
                zoom=zoom,
                show_orientation_axis=image.show_orientation_axis,
                # Unlike Slice, a 3D view draws its own legend inline, so this
                # one needs no second layout cell.
                show_color_bar=show_bar,
            )
            view.add_case_or_clip_to_view(
                self.case, color_map=colour_map,
                field_type=("Cells" if self.spec.sampling.field_type == "CELLS"
                            else "Points"))
            # Visualization3D sets CameraParallelScale to 1/zoom, which frames a
            # 2 m box regardless of how big the geometry is. The scale that
            # actually fits it has to be measured against the camera's own axes;
            # half the largest Cartesian extent clips a body that is not square
            # to the camera.
            view.get_render_view().CameraParallelScale = half_height / zoom
            view.render(layout, row=0, col=0)
            if show_bar:
                self._tune_scalar_bar(view.get_render_view(), log, scale,
                                      bar_width_px=width * image.bar_size)
            self._save(view, layout, filename)
        finally:
            self.case.set_mesh_regions(list(self.spec.sampling.mesh_regions))

    def _save(self, view: Any, layout: Layout, filename: str,
              extra_views: Sequence[Any] = ()) -> None:
        """Save an already-rendered layout, then tear all of it down.

        Args:
            view: the primary view.
            layout: the layout it was rendered into.
            filename: without the ``.png``.
            extra_views: further views occupying cells -- the colour bar --
                that also have to be deleted.
        """
        image = self.spec.image
        if image.background is not None:
            for target in (view, *extra_views):
                render_view = target.get_render_view()
                render_view.UseColorPaletteForBackground = 0
                render_view.Background = list(image.background)
        try:
            Exporter(config=self.config, layout=layout).save_snapshot(
                filename=filename)
        finally:
            # Views and layouts are global ParaView state. Left behind, the
            # next image renders into a grid of every previous one.
            for target in (view, *extra_views):
                target.delete_view()
            layout.delete_layout()

    def close(self) -> None:
        """Drop this case's ParaView state, so the next one starts clean.

        A study renders several cases in one ``pvpython`` process -- starting a
        new one per case would pay ParaView's startup and the reader's mesh scan
        again each time. But readers, views and layouts are global session
        state, and a second case built on top of the first inherits its
        pipeline, its animation timesteps and its render views.

        ``ResetSession`` clears the lot. It is safe here precisely because
        nothing is meant to survive: every image is already written to disk.
        """
        for probe in getattr(self, "_probe_objects", []):
            try:
                pv.Delete(probe)
            except Exception:                             # pragma: no cover
                pass
        self._probe_objects = []
        try:
            pv.ResetSession()
        except Exception:                                 # pragma: no cover
            pass

    # -- the manifest -------------------------------------------------------- #

    def _write_manifest(self, times: Sequence[float]) -> Path:
        path = self.output_dir / "manifest.yaml"
        document = {
            "case": {
                "name": self.facts.name,
                "directory": str(self.facts.directory),
                "length_scale_m": self.length_scale_m,
            },
            "spec": str(self.spec.source) if self.spec.source else None,
            "sampling": {
                "times": [float(time) for time in times],
                "prefer_mean": self.spec.sampling.prefer_mean,
                "field_type": self.spec.sampling.field_type,
            },
            "images": [
                {
                    "file": entry.path.name,
                    "field": entry.field,
                    "requested": entry.requested_field,
                    "plane": entry.plane,
                    "time": float(entry.time),
                    "min": float(entry.value_min),
                    "max": float(entry.value_max),
                    "log": bool(entry.log),
                    "component": entry.component or None,
                    "kind": entry.kind,
                }
                for entry in self.written
            ],
            "skipped": [
                {
                    "field": entry.field,
                    "plane": entry.plane,
                    "time": float(entry.time),
                    "reason": entry.reason,
                }
                for entry in self.skipped
            ],
        }
        path.write_text(
            yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
            encoding="utf-8", newline="\n")
        print(f"    wrote {path.name}")
        return path


def render_case(spec: VizSpec, case_dir: str | Path, *,
                only_fields: Sequence[str] | None = None,
                only_planes: Sequence[str] | None = None,
                allow_constant: bool = False,
                verbose: bool = False) -> Renderer:
    """Render one case end to end.

    Args:
        spec: the validated sampling configuration.
        case_dir: the OpenFOAM case directory.
        only_fields: restrict to these field names.
        only_planes: restrict to these plane names.
        allow_constant: draw constant fields instead of skipping them.
        verbose: list every array the case actually has at each time.

    Returns:
        The :class:`Renderer`, carrying ``written`` and ``skipped``.
    """
    facts = discover_case(case_dir)
    renderer = Renderer(spec, facts, allow_constant=allow_constant,
                        verbose=verbose)
    renderer.run(only_fields=only_fields, only_planes=only_planes)
    return renderer
