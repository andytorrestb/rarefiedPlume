r"""ParaView, driven directly: a reader, a layout, the views in it, one PNG.

**This module imports ParaView.** Like :mod:`plumetools.viz.render`, which is
its only caller, it can be imported by ``pvpython`` and by nothing else.

It needs nothing *but* ParaView. Everything here is ``paraview.simple`` and the
standard library, so the renderer runs under whatever ``pvpython`` is on the
machine -- there is no launcher to install and no version of anything to match
beyond ParaView's own.

:mod:`~plumetools.viz.render` decides *what* to draw: which time, which array,
what range. This module is the part that then talks to ParaView, and it is kept
separate because the two change for different reasons -- that one when a study
wants a different picture, this one when ParaView changes underneath it.

Where the call sequences come from
----------------------------------
This package first drew its images through `VifPara
<https://github.com/virtual-vehicle/VifPara>`_ 1.3.4 (MIT, copyright (c) 2025
Virtual Vehicle Research GmbH), a wrapper over ``paraview.simple`` that
supports ParaView 5.11 to 5.13 only. The classes below are the handful of
things it was used for, rewritten against ``paraview.simple`` itself.

They make **the same ParaView calls in the same order**, deliberately, including
the ones that look redundant. ParaView's session is global state and the order
things are created, shown and activated in decides what a later call sees; a
tidier sequence that drew a subtly different picture would be worse than an
untidy one that does not. Measured rather than assumed: on ParaView 5.11.2 every
slice this module writes is byte-identical to the one VifPara wrote from the
same case and spec, and a surface image differs from VifPara's by no more than
two VifPara runs differed from each other -- a few dozen pixels, by one level.

One ParaView version to draw like
---------------------------------
The images are meant to look the same whichever ParaView draws them, so that a
series started on one machine can be finished on another. ParaView 5.11 is the
reference, because it drew every image this package made before it could run
on anything else. Four things hold a newer ParaView to it:

* ``paraview.compatibility`` is declared as 5.11, ParaView's own mechanism for
  running an older script. It is what keeps a spec's ``Viridis (matplotlib)``
  meaning something after 6.1 renamed the preset to ``Viridis``, and it prints
  one line saying so at start-up.
* the view background is put back to ParaView 5's blue-grey
  (:func:`_draw_like_paraview_5`), which compatibility mode does not do.
* the reader is asked for point data by name, which ParaView 6 no longer
  derives by default (:class:`Case`).
* legend number formats are written as ``printf`` or as ``std::format``,
  whichever the running version reads natively (:func:`legend_number_format`).

Each is either a no-op on 5.11 or restates its default, so the call sequence
there is unchanged. Checked on 5.11.2 and 6.2.0.

Three consequences of the call sequence are worth knowing before changing
anything here.

1. Every exported image is 800 pixels tall
------------------------------------------
``image.height`` sizes the *view*. It fixes the aspect ratio, and how the
height is shared between the picture and its colour bar, but the layout as a
whole is then scaled to :data:`~plumetools.viz.resolve.EXPORT_HEIGHT`.

2. A slice never draws its own legend
-------------------------------------
:class:`SliceView` applies its colour map with the legend hidden, and
:class:`ColourBarView` shows the same cut with the *geometry* hidden. So a slice
with a legend is two views stacked in one layout -- the legend sits under the
picture instead of over part of it. :class:`SurfaceView` is the opposite and
draws its legend inline.

3. Some things are left behind until the session is reset
---------------------------------------------------------
A :class:`Case` shows its reader in a view of its own that nothing ever looks
at, and a :class:`SliceView` drawn without a legend leaves its colour-bar view.
Neither reaches an image, and :meth:`plumetools.viz.render.Renderer.close`
clears both between cases.
"""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from datetime import datetime
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Any, Sequence

import paraview

# Before paraview.simple is imported, which is why this module has to be the
# first of this package to import it. See "One ParaView version to draw like".
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

import paraview.simple as pv                               # noqa: E402
from paraview import servermanager                         # noqa: E402

from plumetools.viz.resolve import (                       # noqa: E402
    layout_shape,
    number_format,
    slice_frame,
)

__all__ = [
    "Case",
    "ColourBarView",
    "ColourMap",
    "Layout",
    "SessionLog",
    "SliceView",
    "SurfaceView",
    "VERSION",
    "legend_number_format",
    "log",
]

#: The ParaView that is actually running, ``(major, minor)`` -- not the
#: compatibility version declared above.
VERSION = (servermanager.vtkSMProxyManager.GetVersionMajor(),
           servermanager.vtkSMProxyManager.GetVersionMinor())


def legend_number_format(digits_after: int, conversion: str,
                         flags: str = "") -> str:
    """A six-wide legend number format in the dialect this ParaView reads.

    ``printf`` below 6.1 and ``std::format`` from it; see
    :func:`plumetools.viz.resolve.number_format`.
    """
    return number_format(ColourMap.DIGITS_BEFORE, digits_after, conversion,
                         flags=flags, std_format=VERSION >= (6, 1))


def _draw_like_paraview_5() -> None:
    """Put back the view background that ParaView 6 changed.

    ParaView 6 replaced the blue-grey view background with a warm grey, and
    declaring compatibility with 5.11 does not bring the old one back. A view
    whose spec sets no ``image.background`` takes the palette's, so the strip
    behind every colour bar -- and everything around a body in a surface view
    -- would change colour with the ParaView version, part way through a study
    if that is when the machine changed.

    Does nothing on ParaView 5, where this palette is already the one loaded.
    The settings are session state, so this is repeated for every case.
    """
    if VERSION >= (6, 0):
        pv.LoadPalette("BlueGrayBackground")


# --------------------------------------------------------------------------- #
# the log
# --------------------------------------------------------------------------- #

class SessionLog:
    """A timestamped record of what ParaView was asked to do.

    Every line goes to stdout, between the renderer's own, and to
    ``<output.log_directory>/<start time>.log`` once :meth:`set_path` has named
    a directory. The file is what is left when a batch job's stdout is not.

    The handlers are rebuilt lazily, on the first message after the path
    changes, so one process rendering several cases writes one file per case.
    All of them carry the same name -- the time the *process* started.
    """

    FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

    def __init__(self) -> None:
        self._path: str | None = None
        self._started: str | None = None
        self._ready = False
        self._logger: logging.Logger | None = None

    def set_path(self, directory: str | Path | None) -> None:
        """Write to *directory* from the next message on; ``None`` for stdout only."""
        self._path = None if directory is None else str(directory)
        self._ready = False

    def info(self, message: str) -> None:
        self._prepare().info(message)

    def warning(self, message: str) -> None:
        self._prepare().warning(message)

    def error(self, message: str) -> None:
        self._prepare().error(message)

    def _prepare(self) -> logging.Logger:
        if self._ready and self._logger is not None:
            return self._logger

        if self._logger is not None:
            for handler in self._logger.handlers[:]:
                self._logger.removeHandler(handler)
                handler.close()

        self._ready = True
        if self._started is None:
            self._started = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        self._logger = logging.getLogger("plumetools.viz")
        self._logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(self.FORMAT)

        stream = logging.StreamHandler(sys.stdout)
        stream.setFormatter(formatter)
        self._logger.addHandler(stream)

        if self._path is not None:
            directory = Path(self._path)
            directory.mkdir(parents=True, exist_ok=True)
            file = RotatingFileHandler(str(directory / f"{self._started}.log"),
                                       maxBytes=5 * 1024 * 1024, backupCount=10)
            file.setFormatter(formatter)
            self._logger.addHandler(file)

        self._logger.info(f"Logging initialized for {sys.argv[0]}.")
        return self._logger


#: The one log of this process.
log = SessionLog()


# --------------------------------------------------------------------------- #
# the case
# --------------------------------------------------------------------------- #

class Case:
    """A reconstructed OpenFOAM case, open in ParaView.

    Attributes:
        reader: the ``OpenFOAMReader`` proxy. Left with every array selected;
            see :mod:`plumetools.viz.render` for why it is never restricted.
    """

    def __init__(self, foam_file: str | Path) -> None:
        log.info("Creating new case.")
        _draw_like_paraview_5()
        # Point data is asked for by name because ParaView 6 stopped deriving
        # it by default, and `field_type: POINTS` has nothing to colour by
        # without it. On ParaView 5 this is the default, restated.
        self.reader = pv.OpenFOAMReader(FileName=str(foam_file), CaseType=1,
                                        Createcelltopointfiltereddata=1)
        self.reader.UpdatePipelineInformation()
        steps = len(self.reader.TimestepValues)
        # Not drawn into any image. See "left behind" in the module docstring.
        self._view = pv.CreateRenderView()
        self._display = pv.Show(self.reader, self._view)
        log.info(f"Case contains {steps} time steps.")

    def set_mesh_regions(self, regions: Sequence[str]) -> None:
        """Choose which of the reader's blocks are loaded, and reload."""
        try:
            self.reader.MeshRegions = list(regions)
            self.reader.UpdatePipeline()
        except AttributeError:
            log.error(f"One of the following mesh regions do not exist in the "
                      f"loaded case: {list(regions)}.")

    def region_names(self) -> list[str]:
        """Every block the reader offers: ``internalMesh``, ``patch/<name>``, ..."""
        info = self.reader.GetProperty("PatchArrayInfo")
        if info is None:
            return []
        # A flat list of (name, selected) pairs.
        return [info[index] for index in range(0, len(info), 2)]


# --------------------------------------------------------------------------- #
# colour
# --------------------------------------------------------------------------- #

@dataclass
class ColourMap:
    """One array, one preset, one range, and the legend that describes them.

    Attributes:
        field: the array to colour by.
        title: the legend's title.
        component: ``""`` for a scalar, ``"Magnitude"``, or ``"X"``/``"Y"``/``"Z"``.
        preset: a ParaView colour preset name.
        low, high: the range. Always explicit -- see
            :mod:`plumetools.viz.render` for why nothing here autoscales.
        log_scale: colour logarithmically.
        association: ``"CELLS"`` or ``"POINTS"``.
        orientation, location, size: where the legend sits in its view.
        format_type, digits_after: how the range labels are written; see
            :func:`legend_number_format`.
    """

    field: str
    title: str
    component: str
    preset: str
    low: float
    high: float
    log_scale: bool
    association: str
    orientation: str
    location: str
    size: float
    format_type: str
    digits_after: int

    #: Legend typeface and sizes. Not configurable; every image shares them.
    FONT = "Times"
    FONT_SIZE = 16
    #: Minimum width of a legend number, in characters.
    DIGITS_BEFORE = 6

    def apply(self, display: Any, render_view: Any, show_bar: bool) -> None:
        """Colour *display* by this map and show or hide its legend.

        ``ColorBy`` runs twice. The first call, by name alone, is what creates
        the *separate* transfer function -- one private to this representation,
        so two images of one array can carry two ranges; the second names the
        association and the component. Both autoscale as they go, which on a
        logarithmic map of a field that touches zero prints ParaView's "Ranges
        not valid for log-space" warning. It is harmless: the explicit range
        set straight afterwards is the one that is drawn.

        The legend is looked up in the **active** view, so the caller has to
        have made *render_view* active first.
        """
        if self.digits_after < 0:
            log.warning(f"Set digits after of color bar are negative "
                        f"({self.digits_after}). Please set a non-negative "
                        f"value. Correcting to 0.")
            self.digits_after = 0

        try:
            pv.ColorBy(display, self.field, separate=True)
        except RuntimeError:
            # The display has no geometry, or no such array. The legend is
            # still wanted -- this is how a colour-bar view gets one.
            pass
        lut = pv.GetColorTransferFunction(self.field, display, separate=True)
        lut.ScalarRangeInitialized = 1.0
        lut.UseLogScale = 1 if self.log_scale else 0
        lut.ApplyPreset(self.preset, True)

        pv.ColorBy(display, (self.association, self.field, self.component))
        lut.RescaleTransferFunction(self.low, self.high)
        display.LookupTable = lut

        bar = pv.GetScalarBar(lut)
        bar.SetPropertyWithName("Title", self.title)
        bar.SetPropertyWithName("ComponentTitle", self.component)
        bar.SetPropertyWithName("TitleFontFamily", self.FONT)
        bar.SetPropertyWithName("LabelFontFamily", self.FONT)
        bar.SetPropertyWithName("TitleFontSize", self.FONT_SIZE)
        bar.SetPropertyWithName("LabelFontSize", self.FONT_SIZE)
        bar.SetPropertyWithName("Orientation", self.orientation)
        bar.SetPropertyWithName("WindowLocation", self.location)
        bar.SetPropertyWithName(
            "RangeLabelFormat",
            legend_number_format(self.digits_after, self.format_type))
        bar.SetPropertyWithName("ScalarBarLength", self.size)

        display.SetScalarBarVisibility(render_view, bool(show_bar))


# --------------------------------------------------------------------------- #
# the layout
# --------------------------------------------------------------------------- #

class Layout:
    """A column of views, exported as one image.

    Args:
        rows: how many views, stacked top to bottom. Every one has to be filled
            by :meth:`place` before :meth:`save`.
    """

    def __init__(self, rows: int = 1) -> None:
        self.proxy = pv.CreateLayout(name="Layout Slice Matrix")
        self._views: dict[int, _View] = {}
        # The layout is a binary tree with a view at each leaf. Splitting the
        # last leaf again each time gives a column, and the cell that was split
        # is the one whose fraction later sets the height of the row above it.
        self._cells = [0]
        self._split_at: dict[int, int] = {}
        for _row in range(rows - 1):
            parent = self._cells[-1]
            first = self.proxy.SplitVertical(parent, 0.5)
            self._cells[-1:] = [first, first + 1]
            self._split_at[first] = self._split_at[first + 1] = parent
        log.info(f"Prepared layout with dimensions: {[1] * rows}.")

    def place(self, row: int, view: "_View") -> None:
        """Put *view* in *row*."""
        cell = self._cells[row]
        self._views[cell] = view
        self.proxy.AssignView(cell, view.render_view)

    def save(self, path: str) -> None:
        """Size the layout from the views in it, and write it to *path*."""
        log.info(f"Exporting snapshot of layout to {path}.")
        width, height, fractions = layout_shape(
            [(self._views[cell].width, self._views[cell].height)
             for cell in self._cells])
        log.info(f"Exported layout image has width {width}px and height "
                 f"{height}px.")

        for cell, fraction in zip(self._cells, fractions):
            self.proxy.SetSplitFraction(self._split_at[cell], fraction)
        self.proxy.SetSize(width, height)
        pv.SaveScreenshot(path, self.proxy)

    def delete(self) -> None:
        log.info(f"Deleting layout {self}. Layout should not be used anymore.")
        pv.Delete(self.proxy)


# --------------------------------------------------------------------------- #
# the views
# --------------------------------------------------------------------------- #

class _View:
    """A render view of a known pixel size, to go in one row of a layout."""

    def __init__(self, width: int, height: int,
                 show_orientation_axis: bool = False) -> None:
        self._width = width
        self._height = height
        self.render_view = pv.CreateRenderView()
        self._display: Any = None
        self._show_orientation_axis = show_orientation_axis

    @property
    def width(self) -> int:
        """Width in pixels. Never less than 1: a layout divides by it."""
        if self._width <= 0:
            log.warning(f"View {self} has width of {self._width}. Will set "
                        f"width to 10px pre rendering.")
            self._width = 10
        return self._width

    @property
    def height(self) -> int:
        """Height in pixels, with the same floor as :attr:`width`."""
        if self._height <= 0:
            log.warning(f"View {self} has height of {self._height}. Will set "
                        f"height to 10px pre rendering.")
            self._height = 10
        return self._height

    def render(self, layout: Layout, row: int = 0) -> None:
        """Draw this view into *row* of *layout*."""
        self._render_inside(layout, row)
        try:
            self.render_view.OrientationAxesVisibility = (
                1 if self._show_orientation_axis else 0)
        except AttributeError:                            # pragma: no cover
            pass

    def _render_inside(self, layout: Layout, row: int) -> None:
        raise NotImplementedError

    def _activate(self) -> None:
        """Make this the active view, which is where ParaView puts a legend."""
        pv.SetActiveView(None)
        pv.SetActiveView(self.render_view)

    def delete(self) -> None:
        log.info(f"Deleting view {self}. View should not be used anymore.")
        pv.Delete(self.render_view)


class ColourBarView(_View):
    """A view holding nothing but a legend, to sit under a :class:`SliceView`.

    It shows the same cut as the slice it belongs to and then hides it again:
    a scalar bar belongs to a representation, so there has to be one, but it
    does not have to be visible.
    """

    def __init__(self, cut: Any, colour_map: ColourMap,
                 origin: Sequence[float], normal: Sequence[float]) -> None:
        super().__init__(100, 100)
        self._cut = cut
        self._colour_map = colour_map
        self._origin = list(origin)
        self._normal = list(normal)

    def resize(self, width: int, height: int) -> None:
        self._width = width
        self._height = height

    def _render_inside(self, layout: Layout, row: int) -> None:
        log.info("Rendering slice color bar.")
        layout.place(row, self)
        self._activate()

        self._cut.SliceType = "Plane"
        self._cut.SliceType.Origin = self._origin
        self._cut.SliceType.Normal = self._normal

        display = pv.Show(self._cut, self.render_view)
        pv.Hide(self._cut, self.render_view)
        self._colour_map.apply(display, self.render_view, True)


class SliceView(_View):
    """A planar cut through *source*, seen square-on and sized from the cut.

    Args:
        source: any ParaView pipeline object -- the reader, or the tail of a
            chain of Calculators, which is what lets a derived field be cut
            like a written one.
        origin, normal: the cutting plane. The camera looks along *normal*.
        camera_up: the camera's up vector.
        height: view height in pixels; the width follows from the cut.
        colour_map: what to colour by.
        zoom, margin_x, margin_y, offset_x, offset_y: framing; see
            :func:`plumetools.viz.resolve.slice_frame`.
        show_orientation_axis: draw ParaView's XYZ widget.

    Attributes:
        colour_bar: the matching :class:`ColourBarView`. Created whether or not
            it is ever rendered.
    """

    def __init__(self, source: Any, origin: Sequence[float],
                 normal: Sequence[float], camera_up: Sequence[float],
                 height: int, colour_map: ColourMap, *, zoom: float = 1.0,
                 margin_x: float = 0, margin_y: float = 0,
                 offset_x: float = 0, offset_y: float = 0,
                 show_orientation_axis: bool = False) -> None:
        # The width is a placeholder until the cut has been measured.
        super().__init__(263, height, show_orientation_axis)
        self._origin = [float(component) for component in origin]
        self._normal = [float(component) for component in normal]
        self._camera_up = [float(component) for component in camera_up]
        self._colour_map = colour_map
        self._cut = pv.Slice(Input=source)
        self.colour_bar = ColourBarView(self._cut, colour_map,
                                        self._origin, self._normal)
        self._frame(zoom, margin_x, margin_y, offset_x, offset_y)

    def _frame(self, zoom: float, margin_x: float, margin_y: float,
               offset_x: float, offset_y: float) -> None:
        """Cut the plane, then point the camera at what came out."""
        self._activate()

        self._cut.SliceType = "Plane"
        self._cut.SliceType.Origin = self._origin
        self._cut.SliceType.Normal = self._normal

        self._display = pv.Show(self._cut, self.render_view)
        self._display.SetRepresentationType("Surface")
        self._display.PointSize = 1.0

        frame = slice_frame(
            self._cut.GetDataInformation().GetBounds(), self._normal,
            self._camera_up, self._height, zoom=zoom, margin_x=margin_x,
            margin_y=margin_y, offset_x=offset_x, offset_y=offset_y)

        view = self.render_view
        view.CameraFocalPoint = list(frame.focal_point)
        view.CameraPosition = list(frame.camera_position)
        view.CameraViewUp = self._camera_up
        view.CameraParallelProjection = 1
        self._width = frame.width
        self._height = frame.height
        view.CameraParallelScale = frame.parallel_scale

    def _render_inside(self, layout: Layout, row: int) -> None:
        log.info(f"Rendering slice to <{row}|0>.")
        layout.place(row, self)
        self._activate()
        self._colour_map.apply(self._display, self.render_view, False)

    def render_colour_bar(self, layout: Layout, row: int, height: int) -> None:
        """Draw the legend into *row*, as wide as the slice and *height* tall."""
        self.colour_bar.resize(self._width, height)
        self.colour_bar.render(layout, row)


class SurfaceView(_View):
    """The reader's own geometry in a parallel 3-D view, for a boundary field.

    Unlike :class:`SliceView` this is not sized for you: the width is explicit,
    and :func:`plumetools.viz.resolve.view_shape` is what derives one that does
    not stretch the geometry.

    Args:
        camera_position, camera_up, focal_point: the camera.
        width, height: view size in pixels.
        zoom: sets the initial parallel scale to ``1 / zoom``, which frames a
            two-unit box whatever the geometry; the caller overrides it.
        show_orientation_axis: draw ParaView's XYZ widget.
        show_colour_bar: draw the legend, inline.
    """

    def __init__(self, camera_position: Sequence[float],
                 camera_up: Sequence[float], focal_point: Sequence[float],
                 width: int, height: int, *, zoom: float = 1.0,
                 show_orientation_axis: bool = False,
                 show_colour_bar: bool = True) -> None:
        super().__init__(width, height, show_orientation_axis)
        self._show_colour_bar = show_colour_bar

        view = self.render_view
        view.CameraFocalPoint = [float(component) for component in focal_point]
        view.CameraPosition = [float(component) for component in camera_position]
        view.CameraViewUp = [float(component) for component in camera_up]
        view.CameraParallelProjection = 1
        view.CameraParallelScale = 1.0 / zoom
        self._activate()

    def add_case(self, case: Case, colour_map: ColourMap,
                 cell_data: bool) -> None:
        """Show whatever blocks *case* currently has loaded, coloured."""
        log.info(f"Adding case {case} to 3D Visualization.")
        self._activate()

        self._display = pv.Show(case.reader, self.render_view)
        colour_map.apply(self._display, self.render_view, self._show_colour_bar)

        self._display.Opacity = 1.0
        self._display.SetRepresentationType("Surface")
        self._display.PointSize = 1.0
        if cell_data:
            self._display.ColorArrayName = ["CELLS", colour_map.field]

    def _render_inside(self, layout: Layout, row: int) -> None:
        log.info(f"Rendering Visualization 3D to <{row}|0>.")
        layout.place(row, self)
        self._activate()
