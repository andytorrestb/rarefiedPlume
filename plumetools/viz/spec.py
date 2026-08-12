r"""The sampling configuration: schema, loading, validation.

A *spec* is a YAML file that says which fields to draw, on which planes, with
what colour scale, and where to put the PNGs. It is the whole user interface of
this package -- adding a field is a list entry, adding a case is a path
argument, and neither requires touching Python.

Nothing here imports ParaView. This module is pure schema, so the test suite can
check the validation without a ParaView installation; see
:mod:`plumetools.viz`.

The shape of a spec
-------------------
.. code-block:: yaml

    version: 1
    extends: default          # optional: inherit from another spec

    sampling: {...}           # which time, averaged or instantaneous
    image:    {...}           # pixel size, colour preset, colour bar
    output:   {...}           # where the PNGs go and what they are called
    planes:   [...]           # the cutting planes
    derived:  [...]           # fields computed from other fields
    fields:   [...]           # what to draw

``fields`` x ``planes`` is a cross product: every field is drawn on every plane,
unless the field names a ``planes:`` subset of its own. A surface field is not
drawn on a plane at all -- see :mod:`plumetools.viz.catalog`.

Unknown keys are rejected, everywhere and always. A misspelled ``prefere_mean``
that was silently ignored would produce instantaneous single-timestep images
that look exactly like averaged ones, which is the sort of mistake that is only
caught months later.

Inheritance
-----------
``extends:`` names either the shipped default (``default``) or a path relative to
the file doing the extending. Mappings (``sampling``, ``image``, ``output``) are
merged key by key; lists (``planes``, ``fields``, ``derived``) are **replaced**
wholesale when present. Replacement rather than concatenation is deliberate: a
case spec that lists three fields should draw three fields, not three plus
whatever the parent happened to carry.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field as dataclass_field
from pathlib import Path
from typing import Any, Iterable, Sequence

import yaml

from plumetools.viz.catalog import FieldKind, catalog_entry, mean_name

#: The spec shipped with the package. ``extends: default`` and a bare
#: ``load_spec()`` both resolve to this file.
DEFAULT_SPEC_PATH = Path(__file__).resolve().parent / "slices.yaml"

#: Schema version this module understands.
SPEC_VERSION = 1

#: Selectors for :attr:`SamplingSpec.time` that are resolved against the case's
#: actual time directories rather than used as a number.
TIME_SELECTORS = ("latest", "first", "all")

#: ParaView associations a colour map may read from.
FIELD_TYPES = ("CELLS", "POINTS")

#: Components a vector field may be coloured by. ``Magnitude`` is the default
#: for every vector in the catalogue.
COMPONENTS = ("", "Magnitude", "X", "Y", "Z")


class VizSpecError(ValueError):
    """A spec file is missing, malformed, or internally inconsistent."""


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def _reject_unknown(raw: dict, allowed: Iterable[str], where: str) -> None:
    """Raise unless every key of *raw* is in *allowed*."""
    unknown = sorted(set(raw) - set(allowed))
    if unknown:
        raise VizSpecError(
            f"{where}: unknown key(s) {unknown}; valid: {sorted(allowed)}")


def _mapping(raw: Any, where: str) -> dict:
    """Coerce *raw* to a mapping, treating a missing section as empty."""
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise VizSpecError(f"{where}: expected a mapping, got {type(raw).__name__}")
    return raw


def _vector3(raw: Any, where: str) -> tuple[float, float, float]:
    """Coerce *raw* to a 3-vector of floats."""
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        raise VizSpecError(f"{where}: expected a list of three numbers, got {raw!r}")
    try:
        return (float(raw[0]), float(raw[1]), float(raw[2]))
    except (TypeError, ValueError):
        raise VizSpecError(f"{where}: expected three numbers, got {raw!r}") from None


def _positive(value: Any, where: str) -> float:
    """Coerce *raw* to a strictly positive float."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise VizSpecError(f"{where}: expected a number, got {value!r}") from None
    if not number > 0.0:
        raise VizSpecError(f"{where}: must be positive, got {number}")
    return number


# --------------------------------------------------------------------------- #
# sections
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class SamplingSpec:
    """Which data gets read out of the case.

    Attributes:
        time: ``"latest"``, ``"first"``, ``"all"``, or an explicit time in
            seconds. Resolved against the case's real time directories, so a
            time that was never written is an error rather than a silent
            snap to the nearest one.
        prefer_mean: reach for ``<field>Mean`` whenever ``fieldAverage`` wrote
            one. On by default, and it should stay on -- an instantaneous
            ``rhoN`` image is a picture of shot noise.
        mesh_regions: the reader's blocks. ``internalMesh`` alone for slices;
            the renderer adds the patch blocks itself for surface fields.
        field_type: ``CELLS`` reproduces the piecewise-constant DSMC samples
            honestly. ``POINTS`` uses ParaView's cell-to-point interpolation,
            which is prettier and smooths over exactly the sampling noise you
            would want an image to show you.
        skip_constant_fields: do not write a PNG for a field that turned out to
            be constant over the geometry being drawn. This is what stops
            ``q`` from producing a uniformly zero image in a case with no wall.
        length_scale_m: the ``L`` of ``origin_over_L``. ``None`` means take it
            from the case -- the nozzle diameter, for the Cai family.
    """

    time: str | float = "latest"
    prefer_mean: bool = True
    mesh_regions: tuple[str, ...] = ("internalMesh",)
    field_type: str = "CELLS"
    skip_constant_fields: bool = True
    length_scale_m: float | None = None

    _KEYS = ("time", "prefer_mean", "mesh_regions", "field_type",
             "skip_constant_fields", "length_scale_m")

    @classmethod
    def from_raw(cls, raw: Any, where: str = "sampling") -> "SamplingSpec":
        raw = _mapping(raw, where)
        _reject_unknown(raw, cls._KEYS, where)

        time = raw.get("time", "latest")
        if isinstance(time, str):
            if time not in TIME_SELECTORS:
                raise VizSpecError(
                    f"{where}.time: {time!r} is not one of {list(TIME_SELECTORS)}, "
                    f"and is not a number")
        else:
            try:
                time = float(time)
            except (TypeError, ValueError):
                raise VizSpecError(f"{where}.time: expected a number or one of "
                                   f"{list(TIME_SELECTORS)}, got {time!r}") from None

        field_type = str(raw.get("field_type", "CELLS")).upper()
        if field_type not in FIELD_TYPES:
            raise VizSpecError(
                f"{where}.field_type: {field_type!r} is not one of {list(FIELD_TYPES)}")

        regions = raw.get("mesh_regions", ["internalMesh"])
        if isinstance(regions, str) or not isinstance(regions, (list, tuple)):
            raise VizSpecError(f"{where}.mesh_regions: expected a list of names")

        length_scale = raw.get("length_scale_m")
        if length_scale is not None:
            length_scale = _positive(length_scale, f"{where}.length_scale_m")

        return cls(
            time=time,
            prefer_mean=bool(raw.get("prefer_mean", True)),
            mesh_regions=tuple(str(r) for r in regions),
            field_type=field_type,
            skip_constant_fields=bool(raw.get("skip_constant_fields", True)),
            length_scale_m=length_scale,
        )


@dataclass(frozen=True)
class ImageSpec:
    """Pixel geometry and colour-bar styling shared by every image.

    Every one of these can be overridden per plane (the geometric ones) or per
    field (the colour ones); this is only the fallback.

    Attributes:
        height: image height in pixels. The width follows from the aspect
            ratio of the slice, so a wide domain gives a wide image.
        zoom: 2.0 fills twice as much of the frame; 0.5 half.
        margin_x, margin_y: pixels added to the computed size. Negative crops.
        offset_x, offset_y: shift the camera within the slice plane, in model
            units, for framing a region rather than the whole cut.
        preset: a ParaView colour preset name. ``"Viridis (matplotlib)"`` is
            the default because it is perceptually uniform and survives being
            printed in greyscale, which ``"Rainbow Desaturated"`` does not.
        color_bar: draw the legend. A VifPara ``Slice`` applies its colour map
            with the legend explicitly hidden and exposes a separate colour-bar
            view instead, so this adds a second row to the layout rather than
            drawing over the image.
        bar_height: height in pixels of that second row.
        orientation, location, bar_size: legend placement within it.
        legend_format_type: ``e`` for exponential -- the right choice when a
            number density runs from 1e10 to 1e20 -- or ``f``, ``g``, ``d``.
        legend_digits_after: digits after the decimal point in the legend.
        show_orientation_axis: draw ParaView's little XYZ widget.
        background: ``[r, g, b]`` in 0..1, or ``None`` for the ParaView default.
    """

    height: int = 700
    zoom: float = 1.0
    margin_x: int = 0
    margin_y: int = 0
    offset_x: float = 0.0
    offset_y: float = 0.0
    preset: str = "Viridis (matplotlib)"
    color_bar: bool = True
    bar_height: int = 110
    orientation: str = "Horizontal"
    location: str = "Lower Center"
    bar_size: float = 0.6
    legend_format_type: str = "e"
    legend_digits_after: int = 1
    show_orientation_axis: bool = False
    background: tuple[float, float, float] | None = None

    _KEYS = ("height", "zoom", "margin_x", "margin_y", "offset_x", "offset_y",
             "preset", "color_bar", "bar_height", "orientation", "location",
             "bar_size", "legend_format_type", "legend_digits_after",
             "show_orientation_axis", "background")

    @classmethod
    def from_raw(cls, raw: Any, where: str = "image") -> "ImageSpec":
        raw = _mapping(raw, where)
        _reject_unknown(raw, cls._KEYS, where)

        background = raw.get("background")
        if background is not None:
            background = _vector3(background, f"{where}.background")

        legend_format = str(raw.get("legend_format_type", "e"))
        if legend_format not in ("e", "f", "g", "d"):
            raise VizSpecError(
                f"{where}.legend_format_type: expected one of 'e', 'f', 'g', 'd', "
                f"got {legend_format!r}")

        return cls(
            height=int(raw.get("height", 700)),
            zoom=_positive(raw.get("zoom", 1.0), f"{where}.zoom"),
            margin_x=int(raw.get("margin_x", 0)),
            margin_y=int(raw.get("margin_y", 0)),
            offset_x=float(raw.get("offset_x", 0.0)),
            offset_y=float(raw.get("offset_y", 0.0)),
            preset=str(raw.get("preset", "Viridis (matplotlib)")),
            color_bar=bool(raw.get("color_bar", True)),
            bar_height=int(raw.get("bar_height", 110)),
            orientation=str(raw.get("orientation", "Horizontal")),
            location=str(raw.get("location", "Lower Center")),
            bar_size=_positive(raw.get("bar_size", 0.6), f"{where}.bar_size"),
            legend_format_type=legend_format,
            legend_digits_after=int(raw.get("legend_digits_after", 1)),
            show_orientation_axis=bool(raw.get("show_orientation_axis", False)),
            background=background,
        )


@dataclass(frozen=True)
class OutputSpec:
    """Where the PNGs land.

    Attributes:
        directory: relative to the case directory unless absolute.
        filename: a template. ``{field}`` is the resolved field name (so
            ``rhoNMean``, not ``rhoN``, when the average was used),
            ``{plane}`` the plane name or ``surface``, ``{case}`` the case
            directory name, ``{time}`` the time as written by the solver.
        log_directory: VifPara's log file. Relative to *directory* unless
            absolute.
        manifest: write ``manifest.yaml`` next to the images, recording what
            was drawn from what. Absent it, a directory of PNGs says nothing
            about which time or which field name produced them.
    """

    directory: str = "results/viz"
    filename: str = "{field}_{plane}"
    log_directory: str = "logs"
    manifest: bool = True

    _KEYS = ("directory", "filename", "log_directory", "manifest")

    #: Placeholders :attr:`filename` may use.
    PLACEHOLDERS = ("field", "plane", "case", "time")

    @classmethod
    def from_raw(cls, raw: Any, where: str = "output") -> "OutputSpec":
        raw = _mapping(raw, where)
        _reject_unknown(raw, cls._KEYS, where)

        filename = str(raw.get("filename", "{field}_{plane}"))
        try:
            probe = filename.format(**{key: "x" for key in cls.PLACEHOLDERS})
        except (KeyError, IndexError) as exc:
            raise VizSpecError(
                f"{where}.filename: {filename!r} uses an unknown placeholder "
                f"{exc}; valid: {list(cls.PLACEHOLDERS)}") from None
        if not probe:
            raise VizSpecError(f"{where}.filename: must not be empty")

        return cls(
            directory=str(raw.get("directory", "results/viz")),
            filename=filename,
            log_directory=str(raw.get("log_directory", "logs")),
            manifest=bool(raw.get("manifest", True)),
        )


@dataclass(frozen=True)
class PlaneSpec:
    """One cutting plane.

    The camera looks **along** ``normal`` with ``camera_up`` pointing up, so
    screen-right is ``normal x camera_up`` -- which means **the sign of the
    normal decides which way the picture reads**, and getting it wrong produces
    a mirror image that is otherwise entirely correct and very easy to miss.
    For Cai's X-Z figures, flow to the right and Z up needs ``normal:
    [0, 1, 0]`` with ``camera_up: [0, 0, 1]``: ``[0,1,0] x [0,0,1] = [1,0,0]``.
    Flipping the normal to ``[0, -1, 0]`` puts the nozzle on the right.

    (VifPara's own ``camera_up x normal`` sizes the viewport. It is not the
    screen axis, and it points the other way.)

    Attributes:
        name: goes into the filename.
        normal: plane normal, and the direction the camera looks along.
        camera_up: which way is up in the image.
        origin: a point on the plane, in metres.
        origin_over_L: the same point in units of
            :attr:`SamplingSpec.length_scale_m`. Mutually exclusive with
            ``origin``.
        height, zoom, margin_x, margin_y, offset_x, offset_y: per-plane
            overrides of :class:`ImageSpec`. ``None`` inherits.
    """

    name: str
    normal: tuple[float, float, float] = (0.0, -1.0, 0.0)
    camera_up: tuple[float, float, float] = (0.0, 0.0, 1.0)
    origin: tuple[float, float, float] | None = None
    origin_over_L: tuple[float, float, float] | None = None
    height: int | None = None
    zoom: float | None = None
    margin_x: int | None = None
    margin_y: int | None = None
    offset_x: float | None = None
    offset_y: float | None = None

    _KEYS = ("name", "normal", "camera_up", "origin", "origin_over_L",
             "height", "zoom", "margin_x", "margin_y", "offset_x", "offset_y")

    @classmethod
    def from_raw(cls, raw: Any, index: int) -> "PlaneSpec":
        where = f"planes[{index}]"
        raw = _mapping(raw, where)
        _reject_unknown(raw, cls._KEYS, where)

        name = raw.get("name")
        if not name:
            raise VizSpecError(f"{where}: needs a 'name'")
        where = f"planes[{name}]"

        if "origin" in raw and "origin_over_L" in raw:
            raise VizSpecError(
                f"{where}: set 'origin' (metres) or 'origin_over_L' (diameters), "
                f"not both")

        origin = _vector3(raw["origin"], f"{where}.origin") if "origin" in raw else None
        origin_over_L = (_vector3(raw["origin_over_L"], f"{where}.origin_over_L")
                         if "origin_over_L" in raw else None)
        if origin is None and origin_over_L is None:
            origin = (0.0, 0.0, 0.0)

        normal = _vector3(raw.get("normal", [0.0, -1.0, 0.0]), f"{where}.normal")
        camera_up = _vector3(raw.get("camera_up", [0.0, 0.0, 1.0]),
                             f"{where}.camera_up")

        if all(component == 0.0 for component in normal):
            raise VizSpecError(f"{where}.normal: must not be the zero vector")
        if all(component == 0.0 for component in camera_up):
            raise VizSpecError(f"{where}.camera_up: must not be the zero vector")

        cross = (
            camera_up[1] * normal[2] - camera_up[2] * normal[1],
            camera_up[2] * normal[0] - camera_up[0] * normal[2],
            camera_up[0] * normal[1] - camera_up[1] * normal[0],
        )
        if all(abs(component) < 1e-12 for component in cross):
            raise VizSpecError(
                f"{where}: camera_up {list(camera_up)} is parallel to normal "
                f"{list(normal)}, so the image has no horizontal axis. Pick an "
                f"up-vector lying in the plane.")

        zoom = raw.get("zoom")
        return cls(
            name=str(name),
            normal=normal,
            camera_up=camera_up,
            origin=origin,
            origin_over_L=origin_over_L,
            height=None if raw.get("height") is None else int(raw["height"]),
            zoom=None if zoom is None else _positive(zoom, f"{where}.zoom"),
            margin_x=None if raw.get("margin_x") is None else int(raw["margin_x"]),
            margin_y=None if raw.get("margin_y") is None else int(raw["margin_y"]),
            offset_x=None if raw.get("offset_x") is None else float(raw["offset_x"]),
            offset_y=None if raw.get("offset_y") is None else float(raw["offset_y"]),
        )

    def origin_m(self, length_scale_m: float | None) -> tuple[float, float, float]:
        """The plane origin in metres.

        Args:
            length_scale_m: the ``L`` of ``origin_over_L``; only consulted when
                this plane used that key.

        Returns:
            ``(x, y, z)`` in metres.

        Raises:
            VizSpecError: if ``origin_over_L`` was used and no length scale is
                available -- guessing 1 m would put the plane somewhere
                arbitrary and still draw a plausible-looking picture.
        """
        if self.origin_over_L is None:
            return self.origin or (0.0, 0.0, 0.0)
        if not length_scale_m:
            raise VizSpecError(
                f"planes[{self.name}] uses origin_over_L, but no length scale is "
                f"known for this case. Set sampling.length_scale_m explicitly.")
        return tuple(component * length_scale_m  # type: ignore[return-value]
                     for component in self.origin_over_L)


@dataclass(frozen=True)
class DerivedSpec:
    """A field computed from other fields by a ParaView Calculator.

    ``dsmcFoam`` writes moments, not the quantities anyone wants to look at:
    velocity is ``momentum/rhoM`` and temperature comes out of ``linearKE``.
    Rather than hard-code those two, the spec carries expressions, so a new
    derived quantity is a YAML entry.

    Expressions are ParaView Calculator syntax over the *resolved* field names,
    which means they must be written against whatever ``prefer_mean`` will
    actually pick. ``{mean}`` expands to ``Mean`` or to nothing accordingly, so
    ``momentum{mean}/rhoM{mean}`` works either way.

    Two further placeholders are filled in from the case rather than the spec,
    so that an expression stays gas-independent:

    ==========  ===========================================================
    ``{R}``     specific gas constant ``k_B / m``, in J/(kg K)
    ``{m}``     molecular mass, in kg
    ==========  ===========================================================

    A spec that uses one of them against a case whose gas is unknown is an
    error, not a substitution of some default -- a temperature computed with
    the wrong gas constant is wrong by a factor nobody would notice in a
    picture.

    Attributes:
        name: the new array's name, and what ``fields:`` refers to.
        expression: the Calculator expression.
        kind: volume or surface, as for a written field.
        component: ``"Magnitude"`` if the result is a vector.
        units: for the colour-bar title.
        log: logarithmic colour scale.
    """

    name: str
    expression: str
    kind: FieldKind = FieldKind.VOLUME
    component: str = ""
    units: str = ""
    log: bool = False

    _KEYS = ("name", "expression", "kind", "component", "units", "log")

    @classmethod
    def from_raw(cls, raw: Any, index: int) -> "DerivedSpec":
        where = f"derived[{index}]"
        raw = _mapping(raw, where)
        _reject_unknown(raw, cls._KEYS, where)

        name = raw.get("name")
        if not name:
            raise VizSpecError(f"{where}: needs a 'name'")
        expression = raw.get("expression")
        if not expression:
            raise VizSpecError(f"derived[{name}]: needs an 'expression'")

        component = str(raw.get("component", ""))
        if component not in COMPONENTS:
            raise VizSpecError(
                f"derived[{name}].component: {component!r} is not one of "
                f"{[c for c in COMPONENTS if c]}")

        kind = str(raw.get("kind", FieldKind.VOLUME.value)).lower()
        try:
            kind_enum = FieldKind(kind)
        except ValueError:
            raise VizSpecError(
                f"derived[{name}].kind: {kind!r} is not 'volume' or 'surface'"
            ) from None

        return cls(
            name=str(name),
            expression=str(expression),
            kind=kind_enum,
            component=component,
            units=str(raw.get("units", "")),
            log=bool(raw.get("log", False)),
        )

    def resolved_expression(self, prefer_mean: bool,
                            substitutions: dict[str, str] | None = None) -> str:
        """The expression with its placeholders expanded.

        Args:
            prefer_mean: fills ``{mean}``.
            substitutions: the case-derived placeholders, e.g. ``{"R": "208.2"}``.

        Returns:
            A ParaView Calculator expression.

        Raises:
            VizSpecError: if a placeholder survives, which means the case could
                not supply it.
        """
        expression = self.expression.replace("{mean}", "Mean" if prefer_mean else "")
        for key, value in (substitutions or {}).items():
            expression = expression.replace(f"{{{key}}}", value)

        leftover = re.findall(r"\{(\w+)\}", expression)
        if leftover:
            raise VizSpecError(
                f"derived[{self.name}]: placeholder(s) {sorted(set(leftover))} "
                f"could not be filled in for this case. '{{R}}' and '{{m}}' need "
                f"the gas to be known -- they are read from the case's case.yaml.")
        return expression


@dataclass(frozen=True)
class FieldSpec:
    """One field to draw, and how.

    Everything but ``name`` is optional; unset keys fall back to
    :mod:`plumetools.viz.catalog` and then to :class:`ImageSpec`.

    Attributes:
        name: the field as the solver writes it, without any ``Mean`` suffix --
            ``prefer_mean`` adds that. Naming ``rhoNMean`` directly also works
            and pins the averaged field regardless of the setting.
        enabled: keep the entry in the file but skip it.
        kind: override the catalogue's volume/surface classification.
        planes: draw only on these planes. ``None`` means all of them.
        prefer_mean: per-field override of :attr:`SamplingSpec.prefer_mean`.
        component: for vectors. Defaults to the catalogue's choice.
        log: logarithmic colour scale. Defaults to the catalogue's choice.
        log_decades: when the scale is logarithmic and the range is
            autoscaled, clamp the bottom to ``max / 10**log_decades``. A DSMC
            plume has cells with *exactly* zero density -- no parcel ever
            reached them -- and a log scale has no colour for zero, so
            something has to set the floor. Six decades is a lot of plume.
        range: ``[min, max]`` to pin the colour scale. Pin it when comparing
            cases: autoscaled images of two Knudsen numbers use two different
            scales and cannot be read side by side.
        preset, legend_title, color_bar, field_type: per-field overrides.
    """

    name: str
    enabled: bool = True
    kind: FieldKind | None = None
    planes: tuple[str, ...] | None = None
    prefer_mean: bool | None = None
    component: str | None = None
    log: bool | None = None
    log_decades: float = 6.0
    range: tuple[float, float] | None = None
    preset: str | None = None
    legend_title: str | None = None
    color_bar: bool | None = None
    field_type: str | None = None

    _KEYS = ("name", "enabled", "kind", "planes", "prefer_mean", "component",
             "log", "log_decades", "range", "preset", "legend_title",
             "color_bar", "field_type")

    @classmethod
    def from_raw(cls, raw: Any, index: int) -> "FieldSpec":
        where = f"fields[{index}]"
        if isinstance(raw, str):          # a bare name is a valid entry
            raw = {"name": raw}
        raw = _mapping(raw, where)
        _reject_unknown(raw, cls._KEYS, where)

        name = raw.get("name")
        if not name:
            raise VizSpecError(f"{where}: needs a 'name'")
        where = f"fields[{name}]"

        kind = raw.get("kind")
        if kind is not None:
            try:
                kind = FieldKind(str(kind).lower())
            except ValueError:
                raise VizSpecError(
                    f"{where}.kind: {kind!r} is not 'volume' or 'surface'") from None

        component = raw.get("component")
        if component is not None:
            component = str(component)
            if component not in COMPONENTS:
                raise VizSpecError(
                    f"{where}.component: {component!r} is not one of "
                    f"{[c for c in COMPONENTS if c]}")

        field_type = raw.get("field_type")
        if field_type is not None:
            field_type = str(field_type).upper()
            if field_type not in FIELD_TYPES:
                raise VizSpecError(
                    f"{where}.field_type: {field_type!r} is not one of "
                    f"{list(FIELD_TYPES)}")

        value_range = raw.get("range")
        if value_range is not None:
            if (not isinstance(value_range, (list, tuple))
                    or len(value_range) != 2):
                raise VizSpecError(
                    f"{where}.range: expected [min, max], got {value_range!r}")
            low, high = float(value_range[0]), float(value_range[1])
            if not high > low:
                raise VizSpecError(
                    f"{where}.range: max ({high}) must exceed min ({low})")
            value_range = (low, high)

        planes = raw.get("planes")
        if planes is not None:
            if isinstance(planes, str) or not isinstance(planes, (list, tuple)):
                raise VizSpecError(f"{where}.planes: expected a list of plane names")
            planes = tuple(str(plane) for plane in planes)

        log_decades = float(raw.get("log_decades", 6.0))
        if not log_decades > 0.0:
            raise VizSpecError(f"{where}.log_decades: must be positive")

        return cls(
            name=str(name),
            enabled=bool(raw.get("enabled", True)),
            kind=kind,
            planes=planes,
            prefer_mean=(None if raw.get("prefer_mean") is None
                         else bool(raw["prefer_mean"])),
            component=component,
            log=None if raw.get("log") is None else bool(raw["log"]),
            log_decades=log_decades,
            range=value_range,
            preset=None if raw.get("preset") is None else str(raw["preset"]),
            legend_title=(None if raw.get("legend_title") is None
                          else str(raw["legend_title"])),
            color_bar=(None if raw.get("color_bar") is None
                       else bool(raw["color_bar"])),
            field_type=field_type,
        )

    # -- resolution against the catalogue ---------------------------------- #

    def resolve_kind(self, derived: dict[str, DerivedSpec]) -> FieldKind:
        """Volume or surface, taking the spec's override then the catalogue."""
        if self.kind is not None:
            return self.kind
        if self.name in derived:
            return derived[self.name].kind
        return catalog_entry(self.name).kind

    def resolve_component(self, derived: dict[str, DerivedSpec]) -> str:
        """The ParaView component to colour by; empty for a scalar."""
        if self.component is not None:
            return self.component
        if self.name in derived:
            return derived[self.name].component
        return catalog_entry(self.name).component

    def resolve_log(self, derived: dict[str, DerivedSpec]) -> bool:
        """Whether to use a logarithmic colour scale."""
        if self.log is not None:
            return self.log
        if self.name in derived:
            return derived[self.name].log
        return catalog_entry(self.name).log

    def resolve_legend_title(self, resolved_name: str,
                             derived: dict[str, DerivedSpec]) -> str:
        """The colour-bar title: the field name, plus units when known."""
        if self.legend_title is not None:
            return self.legend_title
        units = (derived[self.name].units if self.name in derived
                 else catalog_entry(self.name).units)
        return f"{resolved_name}  [{units}]" if units else resolved_name

    def candidate_names(self, prefer_mean: bool,
                        derived: dict[str, DerivedSpec]) -> tuple[str, ...]:
        """Field names to look for in the case, best first.

        A derived field is produced by name and has no averaged counterpart, so
        it is returned as-is. For a written field the averaged name is tried
        first when ``prefer_mean`` holds, and the instantaneous one is the
        fallback -- which is what makes a case that has not reached
        ``fieldAverage``'s ``timeStart`` yet still produce images, with a
        warning, rather than nothing at all.

        Args:
            prefer_mean: the effective setting for this field.
            derived: derived fields by name.

        Returns:
            Names in preference order.
        """
        if self.name in derived:
            return (self.name,)
        averaged = mean_name(self.name)
        if averaged == self.name:      # the spec pinned the averaged field
            return (self.name,)
        if prefer_mean:
            return (averaged, self.name)
        return (self.name, averaged)


@dataclass(frozen=True)
class RenderTask:
    """One PNG: a field, and the plane it is drawn against.

    What the plane means depends on the field. For a volume field it is the
    **cut**. For a surface field nothing is cut -- the boundary patches are
    drawn as they are -- and the plane supplies only the **camera orientation**,
    which is why a surface field takes just one plane by default instead of the
    whole list.
    """

    field: FieldSpec
    plane: PlaneSpec | None

    @property
    def plane_name(self) -> str:
        """The ``{plane}`` placeholder's value."""
        return self.plane.name if self.plane is not None else "surface"


# --------------------------------------------------------------------------- #
# the whole spec
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class VizSpec:
    """A complete, validated sampling configuration."""

    sampling: SamplingSpec = dataclass_field(default_factory=SamplingSpec)
    image: ImageSpec = dataclass_field(default_factory=ImageSpec)
    output: OutputSpec = dataclass_field(default_factory=OutputSpec)
    planes: tuple[PlaneSpec, ...] = ()
    fields: tuple[FieldSpec, ...] = ()
    derived: tuple[DerivedSpec, ...] = ()
    source: Path | None = None

    _KEYS = ("version", "extends", "sampling", "image", "output",
             "planes", "fields", "derived")

    @property
    def derived_by_name(self) -> dict[str, DerivedSpec]:
        """Derived fields keyed by the array name they produce."""
        return {entry.name: entry for entry in self.derived}

    @property
    def planes_by_name(self) -> dict[str, PlaneSpec]:
        """Planes keyed by name."""
        return {plane.name: plane for plane in self.planes}

    def prefer_mean_for(self, field: FieldSpec) -> bool:
        """The effective ``prefer_mean`` for one field."""
        if field.prefer_mean is not None:
            return field.prefer_mean
        return self.sampling.prefer_mean

    def tasks(self, only_fields: Sequence[str] | None = None,
              only_planes: Sequence[str] | None = None) -> list[RenderTask]:
        """Expand the spec into one task per image.

        The cross product of enabled fields and planes, minus the pairs a field
        excluded with its own ``planes:`` list. A surface field is not cut by
        anything, so it gets **one** task -- viewed along the first plane --
        unless it names its own ``planes:``, in which case it gets one view per
        named plane.

        Args:
            only_fields: restrict to these field names (the ``--field`` flag).
            only_planes: restrict to these plane names (the ``--plane`` flag).

        Returns:
            Tasks in spec order, so the images come out in a reproducible order.

        Raises:
            VizSpecError: if a filter or a field's ``planes:`` names something
                the spec does not define.
        """
        derived = self.derived_by_name
        known_planes = self.planes_by_name

        if only_fields is not None:
            known_fields = {field.name for field in self.fields}
            unknown = sorted(set(only_fields) - known_fields)
            if unknown:
                raise VizSpecError(
                    f"no such field(s) in the spec: {unknown}; "
                    f"it has {sorted(known_fields)}")
        if only_planes is not None:
            unknown = sorted(set(only_planes) - set(known_planes))
            if unknown:
                raise VizSpecError(
                    f"no such plane(s) in the spec: {unknown}; "
                    f"it has {sorted(known_planes)}")

        tasks: list[RenderTask] = []
        for field in self.fields:
            if not field.enabled:
                continue
            if only_fields is not None and field.name not in only_fields:
                continue

            if field.planes is not None:
                unknown = sorted(set(field.planes) - set(known_planes))
                if unknown:
                    raise VizSpecError(
                        f"fields[{field.name}].planes names undefined plane(s) "
                        f"{unknown}; the spec defines {sorted(known_planes)}")
                wanted = [known_planes[name] for name in field.planes]
            elif field.resolve_kind(derived) is FieldKind.SURFACE:
                # Not a cut -- the plane only says where the camera stands, and
                # one view of a surface is the sane default.
                wanted = list(self.planes[:1])
            else:
                wanted = list(self.planes)

            for plane in wanted:
                if only_planes is not None and plane.name not in only_planes:
                    continue
                tasks.append(RenderTask(field=field, plane=plane))
        return tasks

    @classmethod
    def from_raw(cls, raw: dict, source: Path | None = None) -> "VizSpec":
        """Validate an already-merged mapping."""
        where = str(source) if source is not None else "spec"
        _reject_unknown(raw, cls._KEYS, where)

        version = raw.get("version", SPEC_VERSION)
        if int(version) != SPEC_VERSION:
            raise VizSpecError(
                f"{where}: version {version} is not supported; this is "
                f"schema version {SPEC_VERSION}")

        planes_raw = raw.get("planes") or []
        if not isinstance(planes_raw, list):
            raise VizSpecError(f"{where}.planes: expected a list")
        planes = tuple(PlaneSpec.from_raw(entry, index)
                       for index, entry in enumerate(planes_raw))
        names = [plane.name for plane in planes]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise VizSpecError(f"{where}.planes: duplicate name(s) {duplicates}")

        fields_raw = raw.get("fields") or []
        if not isinstance(fields_raw, list):
            raise VizSpecError(f"{where}.fields: expected a list")
        fields = tuple(FieldSpec.from_raw(entry, index)
                       for index, entry in enumerate(fields_raw))
        names = [entry.name for entry in fields]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise VizSpecError(f"{where}.fields: duplicate name(s) {duplicates}")

        derived_raw = raw.get("derived") or []
        if not isinstance(derived_raw, list):
            raise VizSpecError(f"{where}.derived: expected a list")
        derived = tuple(DerivedSpec.from_raw(entry, index)
                        for index, entry in enumerate(derived_raw))
        names = [entry.name for entry in derived]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise VizSpecError(f"{where}.derived: duplicate name(s) {duplicates}")

        spec = cls(
            sampling=SamplingSpec.from_raw(raw.get("sampling")),
            image=ImageSpec.from_raw(raw.get("image")),
            output=OutputSpec.from_raw(raw.get("output")),
            planes=planes,
            fields=fields,
            derived=derived,
            source=source,
        )

        # Every field needs a plane -- a volume field to be cut on, a surface
        # field to be looked at from. A spec that quietly produces zero images
        # is worse than one that fails.
        if not spec.planes and any(entry.enabled for entry in spec.fields):
            raise VizSpecError(
                f"{where}: fields are enabled but no planes are defined, so "
                f"nothing would be drawn")
        return spec


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #

def _merge(parent: dict, child: dict) -> dict:
    """Merge *child* over *parent*: mappings key by key, lists replaced."""
    merged = dict(parent)
    for key, value in child.items():
        if (key in merged and isinstance(merged[key], dict)
                and isinstance(value, dict)):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _read(path: Path, seen: tuple[Path, ...]) -> dict:
    """Read one spec file and fold in whatever it extends."""
    if path in seen:
        chain = " -> ".join(str(item) for item in (*seen, path))
        raise VizSpecError(f"circular 'extends': {chain}")
    if not path.is_file():
        raise VizSpecError(f"no spec file at {path}")

    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise VizSpecError(f"{path}: expected a mapping at the top level")

    extends = raw.pop("extends", None)
    if extends is None:
        return raw

    if extends == "default":
        parent_path = DEFAULT_SPEC_PATH
    else:
        parent_path = Path(extends)
        if not parent_path.is_absolute():
            parent_path = (path.parent / parent_path).resolve()
    return _merge(_read(parent_path, (*seen, path)), raw)


def load_spec(path: str | Path | None = None) -> VizSpec:
    """Read, merge and validate a spec file.

    Args:
        path: the YAML file, or ``None`` for the one shipped with the package.

    Returns:
        A validated :class:`VizSpec`.

    Raises:
        VizSpecError: for a missing file, an unknown key, a circular
            ``extends``, or any failed cross-check.
    """
    resolved = DEFAULT_SPEC_PATH if path is None else Path(path)
    return VizSpec.from_raw(_read(resolved, ()), source=resolved)
