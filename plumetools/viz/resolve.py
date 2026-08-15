r"""Deciding what to draw and how to scale it -- the parts with no ParaView in them.

:mod:`plumetools.viz.render` does the drawing and therefore imports ParaView,
which puts it out of reach of the test suite. The decisions it makes on the way
-- which time, which of ``rhoN``/``rhoNMean``, what the colour range should be
given that a DSMC far field is exactly zero -- are ordinary arithmetic, and they
are the parts most likely to be quietly wrong. So they live here, where they can
be tested with no ParaView, no OpenFOAM and no case.
"""

from __future__ import annotations

import math
import re
from typing import Iterable, Sequence

#: Relative width of a range below which the field is called constant. Not
#: exactly zero: a field uniform to twelve significant figures is uniform, and
#: colouring it stretches floating-point noise across the whole palette.
CONSTANT_TOLERANCE = 1e-12

#: Most tick labels to put on a colour bar, however wide it is.
MAX_SCALAR_BAR_LABELS = 5

#: Pixels an exponential tick label needs before it runs into its neighbour.
#: Measured against the 16 pt default: "1.0e+14" is about 60 px wide, and the
#: gap has to be bigger than the label.
SCALAR_BAR_LABEL_PITCH = 110


def scalar_bar_label_budget(bar_width_px: float,
                            limit: int = MAX_SCALAR_BAR_LABELS) -> int:
    """How many decade labels fit on a colour bar *bar_width_px* wide.

    A tall narrow slice -- and a plume down its own axis is exactly that -- gets
    a colour bar only a couple of hundred pixels across, where the decade ticks
    that read perfectly on a wide one overlap into ``1.0e+12.0e+14``. The
    minimum and maximum are drawn separately as range labels and are the two
    that actually have to be legible, so on a narrow bar the answer is none.

    Args:
        bar_width_px: drawn width of the bar itself, not of the view.
        limit: ceiling however wide it gets.

    Returns:
        A label count, possibly 0.
    """
    if bar_width_px <= 0:
        return limit
    return max(0, min(limit, int(bar_width_px // SCALAR_BAR_LABEL_PITCH) - 1))

#: Names a ParaView Calculator expression may use that are not field arrays.
CALCULATOR_FUNCTIONS = frozenset({
    "mag", "abs", "sqrt", "exp", "ln", "log", "log10", "sin", "cos", "tan",
    "asin", "acos", "atan", "sinh", "cosh", "tanh", "ceil", "floor", "min",
    "max", "sign", "cross", "dot", "norm", "iHat", "jHat", "kHat", "coords",
    "coordsX", "coordsY", "coordsZ", "e", "Pi",
})


class RenderError(RuntimeError):
    """The case cannot be rendered: no mesh, no such time, no such field."""


class NothingToRender(RenderError):
    """The case holds no data to draw, because it has not been run yet.

    Deliberately *not* a failure. A study is routinely part-way through -- some
    cases meshed, some solved -- and both families' per-case ``Allpost`` report
    what they can and carry on. A survey image that made ``AllpostCases`` exit
    non-zero for a case nobody has run yet would be stricter than the validation
    step it is bolted onto, which is the wrong way round for a picture.
    """


def resolve_times(available: Sequence[float], wanted: str | float, *,
                  time_min: float | None = None,
                  time_max: float | None = None) -> list[float]:
    """Turn ``sampling.time`` into actual time values.

    An explicit time is matched against what the solver wrote rather than
    snapped to the nearest one: silently rendering 0.0085 s when 0.0099 s was
    asked for would put the wrong time in the manifest and the right-looking
    picture on the screen.

    Args:
        available: the reader's ``TimestepValues``, ascending.
        wanted: ``"latest"``, ``"first"``, ``"all"``, or a time in seconds.
        time_min: drop times below this, or ``None`` for no lower bound.
        time_max: drop times above this, or ``None`` for no upper bound.

    Returns:
        The times to render, ascending.

    Raises:
        RenderError: if the requested time is absent, or is outside the window.
        NothingToRender: if the case has no times at all, or none inside the
            window.

    The window is applied **before** the selector, so ``latest`` means the
    latest time inside it and ``all`` means all of them inside it.

    Why a window exists at all
    --------------------------
    ``fieldAverage`` writes no ``*Mean`` field before its ``timeStart``, so a
    run that writes throughout its transient produces a series whose early
    frames have only the instantaneous field to offer. With
    ``sampling.prefer_mean`` on, those frames silently fall back and the series
    changes quantity part way through -- the inhomogeneity ``docs/viz-slices.md``
    warns about, which is invisible in the resulting pictures because both
    quantities look like a plume.

    ``time_min`` set to ``dsmc.average_start_s`` restricts the series to the
    frames where the average actually exists. That is what
    ``cases/cai2012-health`` does, and it is the only reason a study can turn
    ``prefer_mean`` on over a series at all.
    """
    times = [float(value) for value in available]
    if not times:
        raise NothingToRender(
            "the case has no time directories past 0 -- the solver wrote "
            "nothing. Run ./Allrun first.")

    windowed = [value for value in times
                if (time_min is None or value >= time_min - 1.0e-15)
                and (time_max is None or value <= time_max + 1.0e-15)]
    if not windowed:
        # Not an error: a case part way through its transient has written
        # times but none in the sampling window yet, and a study routinely
        # holds cases at different stages. Same reasoning as NothingToRender.
        raise NothingToRender(
            f"no time in [{_bound(time_min)}, {_bound(time_max, '+inf')}]; the case "
            f"wrote {format(times[0], 'g')} .. {format(times[-1], 'g')}. If "
            f"this case is still running, its sampling window has not been "
            f"reached yet.")

    if wanted == "latest":
        return [windowed[-1]]
    if wanted == "first":
        return [windowed[0]]
    if wanted == "all":
        return windowed

    target = float(wanted)
    for value in windowed:
        if math.isclose(value, target, rel_tol=1e-9, abs_tol=1e-15):
            return [value]
    excluded = any(math.isclose(value, target, rel_tol=1e-9, abs_tol=1e-15)
                   for value in times)
    if excluded:
        raise RenderError(
            f"time {target:g} was written, but the sampling window "
            f"[{_bound(time_min)}, {_bound(time_max, '+inf')}] excludes it.")
    raise RenderError(
        f"no time {target:g} in this case. It wrote: "
        f"{', '.join(format(value, 'g') for value in windowed)}")


def _bound(value: float | None, unbounded: str = "-inf") -> str:
    """A window edge, for a message. ``None`` is unbounded, not zero."""
    return unbounded if value is None else format(float(value), "g")


def resolve_field_name(candidates: Iterable[str],
                       available: Iterable[str]) -> str | None:
    """The first candidate the case actually wrote, or ``None``.

    Args:
        candidates: names in preference order, from
            :meth:`~plumetools.viz.spec.FieldSpec.candidate_names`.
        available: the arrays the case has at this time.

    Returns:
        The name to colour by, or ``None`` if none of them are there.
    """
    present = set(available)
    for candidate in candidates:
        if candidate in present:
            return candidate
    return None


def component_range(flattened: Sequence[float], component: str
                    ) -> tuple[float, float]:
    """Pick one component's ``(min, max)`` out of a flattened range tuple.

    Args:
        flattened: per-component ranges with the magnitude range last, as
            :func:`plumetools.viz.render.array_ranges` stores them.
        component: ``""`` or ``"Magnitude"`` for the magnitude -- the same
            thing as the value range for a scalar -- or ``"X"``, ``"Y"``,
            ``"Z"``.

    Returns:
        ``(min, max)``.
    """
    index = {"X": 0, "Y": 1, "Z": 2}.get(component)
    if index is None or 2 * index + 1 >= len(flattened) - 2:
        return float(flattened[-2]), float(flattened[-1])
    return float(flattened[2 * index]), float(flattened[2 * index + 1])


def is_constant(low: float, high: float) -> bool:
    """Whether a range is too narrow to be worth colouring."""
    return high - low <= CONSTANT_TOLERANCE * max(abs(high), abs(low), 1.0)


def colour_range(low: float, high: float, log: bool, decades: float
                 ) -> tuple[float, float] | None:
    """The range to hand the colour map, or ``None`` if there is nothing to draw.

    On a logarithmic scale the lower end has to be strictly positive, and a DSMC
    plume's far field is *exactly* zero -- no parcel ever reached those cells --
    so the floor comes from *decades* below the maximum instead. Without that,
    ParaView silently substitutes a range starting at 1.0, which for a number
    density in ``m^-3`` throws away most of the plume.

    Args:
        low: the data minimum over the geometry being drawn.
        high: the data maximum.
        log: whether the scale is logarithmic.
        decades: how many decades to show below *high* when the data minimum
            cannot be used.

    Returns:
        ``(min, max)``, or ``None`` when the field is constant (nothing to
        colour) or wholly non-positive on a log scale (nothing to show).
    """
    if not math.isfinite(low) or not math.isfinite(high):
        return None
    if is_constant(low, high):
        return None
    if not log:
        return low, high

    if high <= 0.0:
        return None
    floor = high * 10.0 ** (-decades)
    return (max(low, floor) if low > 0.0 else floor), high


def decade_labels(low: float, high: float,
                  limit: int = MAX_SCALAR_BAR_LABELS) -> list[float]:
    """Powers of ten inside ``[low, high]``, thinned to at most *limit* of them.

    ParaView labels a logarithmic scale at whatever interval it likes, which
    over eight decades puts a dozen labels into a bar a few hundred pixels wide.
    Decade ticks are legible and are what a reader of a log plot expects.

    Args:
        low: bottom of the colour range; must be positive.
        high: top of the colour range.
        limit: how many labels the bar can fit.

    Returns:
        Ascending decade values, or ``[]`` if the range spans no whole decade.
    """
    if low <= 0.0 or high <= low:
        return []
    first = math.ceil(math.log10(low))
    last = math.floor(math.log10(high))
    if last < first:
        return []
    stride = max(1, math.ceil((last - first + 1) / limit))
    return [10.0 ** exponent for exponent in range(first, last + 1, stride)]


#: How far to move a cutting plane that lands exactly on a bounding face, as a
#: fraction of the domain's extent along the plane normal. Small enough to stay
#: well inside the first cell of any mesh this repository generates, large
#: enough to clear floating-point coincidence.
PLANE_NUDGE_FRACTION = 1e-4


def nudge_plane_origin(origin: Sequence[float], normal: Sequence[float],
                       bounds: Sequence[float],
                       fraction: float = PLANE_NUDGE_FRACTION
                       ) -> tuple[tuple[float, float, float], bool]:
    """Move a cutting plane off a bounding face, so the cut is not empty.

    A half domain puts its symmetry plane at the edge of the mesh:
    ``cases/markelov1999`` models ``y >= 0`` and its centre plane is ``y = 0``,
    which is the boundary. VTK's cutter finds no cell interiors there and
    returns an **empty** slice -- so every field reads as "not present on this
    geometry" and no image is drawn, with nothing to suggest the plane was the
    problem.

    Moving the origin a ten-thousandth of the domain into the mesh cuts the
    first cell layer instead. At that distance the picture is the symmetry plane
    to far better than one cell.

    Only axis-aligned normals are adjusted. A tilted plane through a corner
    could be nudged in several directions and none of them is obviously right;
    it is returned untouched, and an empty cut is reported instead of guessed
    at.

    Args:
        origin: the requested plane origin, in metres.
        normal: the plane normal.
        bounds: ``(x_min, x_max, y_min, y_max, z_min, z_max)`` of the data.
        fraction: how far in to move, as a fraction of the extent along the
            normal.

    Returns:
        ``(origin, moved)`` -- the origin to use, and whether it changed.
    """
    axes = [index for index, component in enumerate(normal) if component != 0.0]
    if len(axes) != 1:
        return (float(origin[0]), float(origin[1]), float(origin[2])), False

    axis = axes[0]
    low, high = float(bounds[2 * axis]), float(bounds[2 * axis + 1])
    span = high - low
    if span <= 0.0:
        return (float(origin[0]), float(origin[1]), float(origin[2])), False

    delta = fraction * span
    value = float(origin[axis])
    if value <= low + delta:
        moved_to = low + delta
    elif value >= high - delta:
        moved_to = high - delta
    else:
        return (float(origin[0]), float(origin[1]), float(origin[2])), False

    adjusted = [float(component) for component in origin]
    adjusted[axis] = moved_to
    return (adjusted[0], adjusted[1], adjusted[2]), True


#: Fraction of extra framing left around a 3-D view's geometry.
VIEW_PADDING = 1.06


def _normalise(vector: Sequence[float]) -> tuple[float, float, float]:
    length = math.sqrt(sum(component * component for component in vector))
    if length == 0.0:
        return (0.0, 0.0, 0.0)
    return tuple(component / length for component in vector)  # type: ignore[return-value]


def view_extents(bounds: Sequence[float], normal: Sequence[float],
                 camera_up: Sequence[float]) -> tuple[float, float]:
    """Half-extents of *bounds* projected onto the camera's own axes.

    The camera sees along ``normal`` with ``camera_up`` up, so what ends up on
    screen is the box measured along ``normal x camera_up`` (right) and
    ``camera_up``. Those two numbers -- not any Cartesian extent -- are what
    both the viewport aspect ratio and the parallel scale have to come from.

    Args:
        bounds: ``(x_min, x_max, y_min, y_max, z_min, z_max)``.
        normal: the direction the camera looks along.
        camera_up: the camera's up vector.

    Returns:
        ``(half_width, half_height)`` in world units.
    """
    up = _normalise(camera_up)
    right = _normalise((
        normal[1] * camera_up[2] - normal[2] * camera_up[1],
        normal[2] * camera_up[0] - normal[0] * camera_up[2],
        normal[0] * camera_up[1] - normal[1] * camera_up[0],
    ))

    corners = [(bounds[0 + i], bounds[2 + j], bounds[4 + k])
               for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    along_up = [sum(c * u for c, u in zip(corner, up)) for corner in corners]
    along_right = [sum(c * r for c, r in zip(corner, right)) for corner in corners]

    return (0.5 * (max(along_right) - min(along_right)),
            0.5 * (max(along_up) - min(along_up)))


def view_shape(bounds: Sequence[float], normal: Sequence[float],
               camera_up: Sequence[float], height: int,
               padding: float = VIEW_PADDING) -> tuple[int, float]:
    """Viewport width and parallel scale that show *bounds* undistorted.

    A 3-D view is not sized for you the way a ``Slice`` is: VifPara takes an
    explicit width, and anything fixed there -- ``height * 4 / 3`` was the
    obvious guess -- stretches the geometry by whatever the difference is. The
    cylinder of ``cases/markelov1999`` is half again as tall as it is wide, and
    a 4:3 frame drew it nearly twice too wide.

    Deriving the width from the projected extents instead makes one world unit
    the same number of pixels across and down, which is the whole of what
    "represents the geometry" means here.

    Args:
        bounds: ``(x_min, x_max, y_min, y_max, z_min, z_max)``.
        normal: the direction the camera looks along.
        camera_up: the camera's up vector.
        height: viewport height in pixels.
        padding: fraction of extra room to leave around the geometry.

    Returns:
        ``(width_px, camera_parallel_scale)``. The scale is the half-height in
        world units, which is what ParaView's ``CameraParallelScale`` means.
    """
    half_width, half_height = view_extents(bounds, normal, camera_up)
    if half_height <= 0.0 or half_width <= 0.0:
        # A flat or degenerate body: keep it square rather than divide by zero.
        return max(int(height), 1), (max(half_width, half_height) * padding) or 1.0

    width = int(round(height * half_width / half_height))
    return max(width, 1), half_height * padding


def view_half_height(bounds: Sequence[float], normal: Sequence[float],
                     camera_up: Sequence[float], aspect: float,
                     padding: float = VIEW_PADDING) -> float:
    """The ``CameraParallelScale`` that fits *bounds* into a FIXED aspect view.

    Kept for a viewport whose width is not free to follow the geometry. Where
    it is, :func:`view_shape` is the better answer: it removes the distortion
    rather than padding around it.

    Args:
        bounds: ``(x_min, x_max, y_min, y_max, z_min, z_max)``.
        normal: the direction the camera looks along.
        camera_up: the camera's up vector.
        aspect: viewport width divided by height.
        padding: fraction of extra room to leave around the geometry.

    Returns:
        The half-height, always strictly positive.
    """
    half_width, half_height = view_extents(bounds, normal, camera_up)
    fitted = max(half_height, half_width / aspect if aspect > 0 else half_width)
    return (fitted * padding) or 1.0


def expression_fields(expression: str) -> set[str]:
    """The field arrays a Calculator expression refers to.

    Bare identifiers, minus the Calculator's own function and constant names.
    Used to check a derived field's inputs *before* handing ParaView an
    expression it cannot evaluate, which otherwise fails deep inside VTK with a
    message that names neither the field nor the spec entry.

    Args:
        expression: a resolved Calculator expression.

    Returns:
        The identifiers that must exist as arrays.
    """
    identifiers = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression))
    return identifiers - CALCULATOR_FUNCTIONS
