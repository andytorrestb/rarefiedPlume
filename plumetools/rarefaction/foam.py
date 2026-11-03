r"""OpenFOAM side of the diagnostics: read a sampled time, write fields back into it.

:func:`load` turns one reconstructed time directory of a ``dsmcFoam`` case
into the inputs :mod:`~plumetools.rarefaction.core` works on;
:func:`write_fields` writes the results into the **same** time directory, as
ordinary ``volScalarField`` files. ParaView's OpenFOAM reader and OpenFOAM's
own utilities then find them beside ``rhoNMean`` -- on the mesh, and at the
time, they were computed from -- with nothing to convert.

What is read, and what is checked
---------------------------------
* The averaged fields named by ``settings.field_names`` (``rhoNMean`` ...), each
  checked to have exactly one value per cell of ``constant/polyMesh``.
* The averaging window, from ``<time>/uniform/functionObjects/
  functionObjectProperties``, where ``fieldAverage`` keeps ``totalIter`` and
  ``totalTime`` for every field it averages. Fields averaged over **different**
  windows -- a field added to ``fieldAverage`` on a restart, say -- cannot be
  combined, and are refused rather than mixed.
* The particle weight implied by the fields themselves,
  ``rhoNMean * V / dsmcRhoNMean``. It is one number in every occupied cell when
  the count and the density are what they are taken to be, and it should be
  ``nEquivalentParticles``. A spread, or a mismatch, means a field mapping is
  wrong, and the summary says so.
* ``C`` and ``V`` come from ``postProcess -func writeCellCentres`` and
  ``writeCellVolumes`` -- OpenFOAM's own centroids and volumes, which
  ``./Allpost`` writes.
* Cell bounding boxes and face neighbours come from ``constant/polyMesh``,
  read whole with the array readers in :mod:`plumetools.mesh.polymesh`.

Only a reconstructed case is read; ``processor*`` directories are not.

What is written
---------------
One file per output field, ``internalField nonuniform List<scalar>``, with
``zeroGradient`` on every ordinary patch and the patch's own type on a
constraint patch. A cell with nothing to compute from holds
``settings.invalid_value``; a low-confidence cell keeps its value, and
``diagValid`` marks it. A file of that
name that this module did not write is never overwritten.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from plumetools.foamio.errors import FoamDataError
from plumetools.foamio.patchfields import patch_field_type
from plumetools.foamio.primitives import BANNER, SEPARATOR
from plumetools.foamio.times import latest_time
from plumetools.foamio.volfield import read_internal
from plumetools.mesh.boundary import read_boundary
from plumetools.mesh.patches import face_geometry
from plumetools.mesh.polymesh import read_faces_flat, read_label_list, read_points_array
from plumetools.rarefaction.core import CellGraph

GENERATOR = "plumetools.rarefaction"

#: The shared readers' error, under the name this module has always exported.
PostError = FoamDataError

#: Patch types across which a cell sees its own mirror image.
MIRROR_PATCH_TYPES = ("symmetry", "symmetryPlane", "empty")

#: Relative spread of the implied particle weight that is still "one number".
WEIGHT_TOLERANCE = 1.0e-6

__all__ = ["CaseData", "PostError", "SamplingWindow", "latest_time", "load",
           "write_fields"]


@dataclass(frozen=True)
class SamplingWindow:
    """The ``fieldAverage`` window the averaged fields were accumulated over.

    Attributes:
        end_s: the time read -- the end of the window [s].
        duration_s: ``totalTime`` [s].
        n_steps: ``totalIter``, the time steps averaged.
        per_field: ``{field: (totalIter, totalTime)}`` as recorded.
    """

    end_s: float
    duration_s: float
    n_steps: int
    per_field: dict

    @property
    def start_s(self) -> float:
        return self.end_s - self.duration_s

    def as_dict(self) -> dict:
        return {"start_s": self.start_s, "end_s": self.end_s,
                "duration_s": self.duration_s, "averaged_steps": self.n_steps,
                "fields": sorted(self.per_field)}


@dataclass
class CaseData:
    """Everything read from one time directory.

    Attributes:
        case_dir, time: where it came from.
        n_cells: cells in ``constant/polyMesh``.
        inputs: ``{input name: array or object}`` for :class:`core.Context`.
        unavailable: ``{input name: why it could not be read}``.
        window: the sampling window, or ``None`` if none is recorded.
        patches: the mesh's patches, as :func:`read_boundary` gives them.
        checks: consistency checks, for the summary.
    """

    case_dir: Path
    time: str
    n_cells: int
    inputs: dict
    unavailable: dict
    window: SamplingWindow | None
    patches: dict
    checks: dict = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# reading
# --------------------------------------------------------------------------- #

def _cell_values(path: Path, n_cells: int) -> np.ndarray:
    return read_internal(path).over(n_cells, "cells in constant/polyMesh")


def read_window(case_dir: Path, time: str, base_names) -> tuple:
    """The averaging window of ``base_names``, from the function-object state.

    Returns:
        ``(window or None, reason)`` -- the reason being why there is none.

    Raises:
        PostError: if the fields were averaged over different windows.
    """
    path = (Path(case_dir) / time / "uniform" / "functionObjects"
            / "functionObjectProperties")
    if not path.is_file():
        return None, (f"no {path.relative_to(case_dir)}, where fieldAverage "
                      f"records how many steps it averaged")
    text = path.read_text(encoding="utf-8", errors="replace")
    found = {}
    for name in base_names:
        match = re.search(
            rf"^\s*{re.escape(name)}\s*\{{[^{{}}]*?totalIter\s+(\d+)\s*;"
            rf"[^{{}}]*?totalTime\s+([^;\s]+)\s*;", text, re.M | re.S)
        if match:
            found[name] = (int(match.group(1)), float(match.group(2)))
    absent = [name for name in base_names if name not in found]
    if absent:
        return None, (f"{path.relative_to(case_dir)} records no averaging "
                      f"window for {absent}")
    steps = {iters for iters, _ in found.values()}
    durations = [duration for _, duration in found.values()]
    if len(steps) > 1 or max(durations) - min(durations) > 1e-9 * max(durations):
        raise PostError(
            f"{path}: the fields were averaged over different windows -- "
            + ", ".join(f"{n}: {i} steps / {t:.6g} s" for n, (i, t) in found.items())
            + ". Means from different samples cannot be combined into one "
              "temperature or mean free path.")
    iters, duration = next(iter(found.values()))
    return SamplingWindow(end_s=float(time), duration_s=duration, n_steps=iters,
                          per_field=found), ""


def _particle_weight(case_dir: Path) -> float | None:
    path = Path(case_dir) / "constant" / "dsmcProperties"
    if not path.is_file():
        return None
    match = re.search(r"^\s*nEquivalentParticles\s+([^;\s]+)\s*;",
                      path.read_text(encoding="utf-8", errors="replace"), re.M)
    return float(match.group(1)) if match else None


def weight_check(number_density, volume, parcels_per_step,
                 configured: float | None) -> dict:
    """Whether ``n V / parcels`` is one number, and the configured one.

    That identity is what makes ``dsmcRhoN`` a parcel COUNT and ``rhoN`` a
    molecule density; see :mod:`plumetools.rarefaction.core`.
    """
    occupied = (parcels_per_step > 0.0) & (number_density > 0.0)
    if not occupied.any():
        return {"checked": False, "reason": "no occupied cell"}
    implied = number_density[occupied] * volume[occupied] / parcels_per_step[occupied]
    median = float(np.median(implied))
    spread = float((implied.max() - implied.min()) / median)
    document = {"checked": True, "implied_particle_weight": median,
                "relative_spread": spread,
                "configured_particle_weight": configured}
    consistent = spread <= WEIGHT_TOLERANCE
    if configured:
        document["relative_difference"] = abs(median - configured) / configured
        consistent &= document["relative_difference"] <= WEIGHT_TOLERANCE
    document["consistent"] = bool(consistent)
    return document


def _mesh_topology(case_dir: Path) -> tuple:
    poly = Path(case_dir) / "constant" / "polyMesh"
    for name in ("owner", "neighbour", "faces", "points", "boundary"):
        if not (poly / name).is_file():
            raise PostError(f"no {poly / name}; the case has no mesh to read")
    owner = read_label_list(poly / "owner")
    neighbour = read_label_list(poly / "neighbour")
    n_cells = int(max(owner.max(initial=-1), neighbour.max(initial=-1))) + 1
    return owner, neighbour, n_cells


def cell_extents(points, offsets, labels, owner, neighbour, n_cells) -> np.ndarray:
    """``(n_cells, 3)`` bounding-box dimensions of every cell [m].

    Per face first -- ``reduceat`` over the CSR face list -- then per cell over
    the faces it owns or neighbours.
    """
    starts = offsets[:-1]
    n_internal = neighbour.size
    extent = np.empty((n_cells, 3))
    for axis in range(3):
        coords = points[labels, axis]
        face_lo = np.minimum.reduceat(coords, starts)
        face_hi = np.maximum.reduceat(coords, starts)
        lo = np.full(n_cells, np.inf)
        hi = np.full(n_cells, -np.inf)
        np.minimum.at(lo, owner, face_lo)
        np.minimum.at(lo, neighbour, face_lo[:n_internal])
        np.maximum.at(hi, owner, face_hi)
        np.maximum.at(hi, neighbour, face_hi[:n_internal])
        extent[:, axis] = hi - lo
    return extent


def cell_graph(centres, points, offsets, labels, owner, neighbour,
               patches) -> CellGraph:
    """Face neighbours, plus a mirror image across every symmetry/empty face."""
    mirror_cells, mirror_offsets = [], []
    for info in patches.values():
        if info.type not in MIRROR_PATCH_TYPES or info.n_faces == 0:
            continue
        index = range(info.start_face, info.start_face + info.n_faces)
        faces = [labels[offsets[i]:offsets[i + 1]] for i in index]
        face_centres, area_vectors = face_geometry(faces, points)
        normals = area_vectors / np.linalg.norm(area_vectors, axis=1)[:, None]
        cells = owner[info.start_face:info.start_face + info.n_faces]
        height = np.einsum("ij,ij->i", face_centres - centres[cells], normals)
        mirror_cells.append(cells)
        mirror_offsets.append(2.0 * height[:, None] * normals)
    n_internal = neighbour.size
    return CellGraph(
        centres=centres, owner=owner[:n_internal], neighbour=neighbour,
        mirror_cells=(np.concatenate(mirror_cells) if mirror_cells
                      else np.empty(0, dtype=np.int64)),
        mirror_offsets=(np.concatenate(mirror_offsets) if mirror_offsets
                        else np.empty((0, 3))))


def load(case_dir: Path, settings, *, time: str | None = None) -> CaseData:
    """Read one time directory into diagnostic inputs.

    A missing averaged field or cell-geometry file is not an error here: the
    input is listed as unavailable, with the command that would supply it, and
    the diagnostics that do not need it still run.

    Raises:
        PostError: for no time directory, no mesh, a field from another mesh,
            or fields averaged over different windows.
    """
    case_dir = Path(case_dir)
    time = time or latest_time(case_dir)
    directory = case_dir / time
    if not directory.is_dir():
        raise PostError(f"no time directory {directory}")

    owner, neighbour, n_cells = _mesh_topology(case_dir)
    patches = read_boundary(case_dir)
    inputs: dict = {}
    unavailable: dict = {}

    def averaged(role: str):
        name = settings.input_field(role)
        path = directory / name
        if not path.is_file():
            unavailable[role] = (
                f"no {name} at t = {time}. fieldAverage writes it once the run "
                f"passes its timeStart (dsmc.average_start_s); an earlier time "
                f"has only one step's samples")
            return
        inputs[role] = _cell_values(path, n_cells)

    for role in ("number_density", "mass_density", "momentum_density",
                 "kinetic_energy_density", "parcels_per_step"):
        averaged(role)

    for role, function in (("cell_centres", "writeCellCentres"),
                           ("cell_volumes", "writeCellVolumes")):
        name = settings.input_field(role)
        path = directory / name
        if path.is_file():
            inputs[role] = _cell_values(path, n_cells)
        else:
            unavailable[role] = (f"no {name} at t = {time}; postProcess -func "
                                 f"{function} -time {time} writes it, and "
                                 f"./Allpost runs that")
    if "cell_volumes" in inputs:
        inputs["cell_volume"] = inputs.pop("cell_volumes")
    else:
        unavailable["cell_volume"] = unavailable["cell_volumes"]

    # The window: only fieldAverage outputs have one.
    base_names = []
    for role in ("number_density", "mass_density", "momentum_density",
                 "kinetic_energy_density", "parcels_per_step"):
        name = settings.input_field(role)
        if role in inputs and name.endswith("Mean"):
            base_names.append(name[:-len("Mean")])
    window, reason = (read_window(case_dir, time, base_names) if base_names
                      else (None, "no averaged input was read"))
    if window is not None:
        inputs["averaged_steps"] = window.n_steps
    else:
        unavailable["averaged_steps"] = (
            f"{reason}. particle_count_basis: accumulated needs it; per_step "
            f"does not")

    # The mesh, read whole once and shared by the extents and the graph.
    points = read_points_array(case_dir)
    offsets, labels = read_faces_flat(case_dir)
    inputs["cell_extent"] = cell_extents(points, offsets, labels, owner,
                                         neighbour, n_cells)
    if "cell_centres" in inputs:
        inputs["cell_graph"] = cell_graph(inputs["cell_centres"], points, offsets,
                                          labels, owner, neighbour, patches)
    else:
        unavailable["cell_graph"] = unavailable["cell_centres"]

    checks = {}
    if {"number_density", "parcels_per_step", "cell_volume"} <= set(inputs):
        checks["particle_weight"] = weight_check(
            inputs["number_density"], inputs["cell_volume"],
            inputs["parcels_per_step"], _particle_weight(case_dir))

    return CaseData(case_dir=case_dir, time=time, n_cells=n_cells, inputs=inputs,
                    unavailable=unavailable, window=window, patches=patches,
                    checks=checks)


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #

#: Values formatted and written per chunk: bounds the strings alive at once.
WRITE_CHUNK = 65536


def _write_streamed(path: Path, head: list[str], values: np.ndarray,
                    tail: list[str]) -> Path:
    """``head``, one line per value, ``tail`` -- streamed, LF line endings.

    Byte for byte what ``write_lines(path, head + [one line per value] +
    tail)`` writes, without holding a formatted string for every cell: a field
    of a few million cells is otherwise a few hundred megabytes of Python
    strings.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(head) + "\n")
        for start in range(0, values.size, WRITE_CHUNK):
            chunk = values[start:start + WRITE_CHUNK].tolist()
            handle.write("\n".join(f"{v:.10g}" for v in chunk) + "\n")
        handle.write("\n".join(tail) + "\n")
    return path


def _render(name: str, time: str, output, n_values: int, patches,
            invalid_value: float) -> tuple[list[str], list[str]]:
    """The lines before and after the values of one output field."""
    note = (f"Cells with no value hold {invalid_value:g}; diagValid = 0 marks"
            f" the rest that are not trusted, and diagStatus says why."
            if output.kind == "value" else "Defined in every cell.")
    head = [
        BANNER,
        f"// GENERATED by {GENERATOR} -- do not edit. Rewritten by each run of",
        "// the diagnostics. Units: " + output.units + ".",
        f"// {output.description}. {note}",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       volScalarField;",
        f'    location    "{time}";',
        f"    object      {name};",
        "}",
        SEPARATOR,
        "",
        f"dimensions      {output.dimensions};",
        "",
        "internalField   nonuniform List<scalar>",
        str(n_values),
        "(",
    ]
    tail = [")", ";", "", "boundaryField", "{"]
    for patch, info in patches.items():
        tail += [f"    {patch}", "    {",
                 f"        type            {patch_field_type(info.type)};",
                 "    }"]
    tail += ["}", "", "// " + "*" * 73 + " //"]
    return head, tail


def write_fields(data: CaseData, outputs: dict, registry: dict,
                 invalid_value: float) -> list[Path]:
    """Write each computed output into the time directory it was read from.

    Raises:
        PostError: rather than overwrite a file this module did not write.
    """
    directory = data.case_dir / data.time
    written = []
    for name, result in outputs.items():
        output = registry[name].output
        path = directory / name
        if path.is_file():
            with open(path, encoding="utf-8", errors="replace") as handle:
                head = handle.read(2048)
            if f"GENERATED by {GENERATOR}" not in head:
                raise PostError(
                    f"{path} exists and was not written by {GENERATOR}; it is "
                    f"not overwritten. Remove it, or name this diagnostic "
                    f"differently.")
        values = result.values
        if output.kind == "value":
            values = np.where(np.isfinite(values), values, float(invalid_value))
        values = np.asarray(values, dtype=np.float64)
        head, tail = _render(name, data.time, output, values.size, data.patches,
                             invalid_value)
        written.append(_write_streamed(path, head, values, tail))
    return written
