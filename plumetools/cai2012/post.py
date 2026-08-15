r"""Extract the four validation outputs, and score them against Cai.

§9 of the case specification asks for exactly four things, and this module
produces those and nothing else:

======  ===================================  =============================
A       centreline number density            ``n/n0`` vs ``X/D``   (Figs. 19-21)
B       centreline axial velocity            ``U1 sqrt(beta0)`` vs ``X/D``
C       centreline temperature               ``T/T0`` vs ``X/D``
D       number-density contour               ``n/n0`` in a centre plane (Figs. 6-8)
======  ===================================  =============================

Reading the solver's output
---------------------------
Everything comes from the **time-averaged** fields ``fieldAverage`` writes:
``rhoNMean``, ``momentumMean``, ``linearKEMean``, ``rhoMMean``. A single-timestep
``rhoN`` is one step's worth of parcels -- ``DSMCCloud::resetFields`` zeroes it
every step -- and comparing that with an analytical curve measures shot noise.

Cell centres and volumes are not in the mesh files in a form worth re-deriving,
so ``./Allpost`` writes them with OpenFOAM's own function objects::

    postProcess -func writeCellCentres    -> C, Cx, Cy, Cz
    postProcess -func writeCellVolumes    -> V

and this module reads them back. Volumes matter: the mesh is graded, so an
unweighted mean over a centreline tube would let one huge outer cell outvote a
hundred core cells.

Deriving U and T from what dsmcFoam writes
------------------------------------------
``dsmcFoam`` writes the moments, not the derived quantities:

.. code-block:: text

    rhoN       number density                  n
    rhoM       mass density                    rho = m n
    momentum   momentum density                rho u
    linearKE   translational kinetic energy density   (1/2) rho <v^2>

so

.. code-block:: text

    u    = momentum / rhoM
    <v^2> = 2 linearKE / rhoM
    3 R T = <v^2> - |u|^2

which is the same definition :mod:`plumetools.cai2012.analytical` uses, so the
DSMC and analytical temperatures are the same quantity rather than two things
with the same name.

Statistics, not pointwise equality
----------------------------------
Cai puts the worst-case centreline scatter at roughly 1%. DSMC results are
random variables; the metrics here are maximum, mean and RMS **relative** error
over the sampled range, and the comparison in ``docs/cai2012-case.md`` is against
Cai's reported error *trend*, not against his numbers as tolerances.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.cai2012 import analytical

#: Averaged fields this module reads, and whether each is a vector.
REQUIRED_FIELDS = {
    "rhoNMean": False,
    "rhoMMean": False,
    "momentumMean": True,
    "linearKEMean": False,
}


class PostError(ValueError):
    """A required field or time directory is missing, or is unreadable."""


# --------------------------------------------------------------------------- #
# reading OpenFOAM fields
# --------------------------------------------------------------------------- #

_LIST_HEAD = re.compile(r"internalField\s+(uniform|nonuniform)", re.S)


def read_internal_field(path: Path) -> np.ndarray:
    """Read the ``internalField`` of an OpenFOAM volume field.

    Args:
        path: the field file, e.g. ``<case>/0.05/rhoNMean``.

    Returns:
        ``(n_cells,)`` for a scalar field or ``(n_cells, 3)`` for a vector one.
        A ``uniform`` entry becomes a one-element array; the caller knows the
        cell count and can broadcast.

    Raises:
        PostError: if the file is missing or has no readable internal field.

    A purpose-built parser rather than a general dictionary reader: it only has
    to find one list, and it has to survive the ``boundaryField`` block that
    follows -- which is where the legacy line-offset parsers in this repository
    went wrong (findings RB-03, RB-04).
    """
    path = Path(path)
    if not path.is_file():
        raise PostError(f"no field file at {path}")
    text = path.read_text(encoding="utf-8", errors="replace")

    head = _LIST_HEAD.search(text)
    if not head:
        raise PostError(f"{path}: no internalField entry")

    if head.group(1) == "uniform":
        token = re.search(r"internalField\s+uniform\s+(\([^)]*\)|[-+0-9.eE]+)\s*;",
                          text)
        if not token:
            raise PostError(f"{path}: malformed uniform internalField")
        value = token.group(1)
        if value.startswith("("):
            return np.asarray([[float(v) for v in re.findall(r"[-+0-9.eE]+", value)]],
                              dtype=np.float64)
        return np.asarray([float(value)], dtype=np.float64)

    kind = re.search(r"internalField\s+nonuniform\s+List<(\w+)>", text)
    if not kind:
        raise PostError(f"{path}: nonuniform internalField with no List<> type")

    start = text.index("(", head.start())
    depth, i = 0, start
    while i < len(text):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                break
        i += 1
    payload = text[start + 1:i]
    values = np.fromstring(payload.replace("(", " ").replace(")", " "), sep=" ")
    return values.reshape(-1, 3) if kind.group(1) == "vector" else values


def latest_time(case_dir: Path) -> str:
    """The largest numeric time directory in a case.

    Raises:
        PostError: if the case has no time directory past ``0`` -- which means
            the solver never wrote anything, and there is nothing to
            post-process.
    """
    case_dir = Path(case_dir)
    times = []
    for entry in case_dir.iterdir():
        if not entry.is_dir():
            continue
        try:
            value = float(entry.name)
        except ValueError:
            continue
        if value > 0.0:
            times.append((value, entry.name))
    if not times:
        raise PostError(
            f"{case_dir} has no time directory past 0; the solver wrote nothing. "
            f"Run ./Allrun first.")
    return max(times)[1]


@dataclass(frozen=True)
class SampledField:
    """The averaged DSMC fields, on cell centres.

    Attributes:
        time: the time directory read.
        centres: ``(n, 3)`` cell centres [m].
        volumes: ``(n,)`` cell volumes [m^3].
        number_density: ``(n,)`` ``n`` [1/m^3].
        velocity: ``(n, 3)`` bulk velocity [m/s].
        temperature: ``(n,)`` translational temperature [K].
    """

    time: str
    centres: np.ndarray
    volumes: np.ndarray
    number_density: np.ndarray
    velocity: np.ndarray
    temperature: np.ndarray

    @property
    def n_cells(self) -> int:
        return int(self.centres.shape[0])


def read_case(case_dir: Path, species_mass_kg: float, *,
              time: str | None = None,
              mesh_time: str | None = None) -> SampledField:
    """Read the averaged fields and derive velocity and temperature.

    Args:
        case_dir: the case directory.
        species_mass_kg: molecular mass, for ``R = k_B/m``.
        time: the time directory, or ``None`` for the latest.
        mesh_time: where to read ``C`` and ``V`` from, if not from *time*.
            ``None`` reads them alongside the fields, which is what a single
            post-processing pass wants.

    Returns:
        A :class:`SampledField`.

    Raises:
        PostError: for a missing time directory, a missing field, or missing
            cell centres. The last one names the command that writes them,
            because it is a post-processing step people forget.

    Reading the geometry from another time
    --------------------------------------
    ``C`` and ``V`` are properties of the mesh, and the mesh does not move: the
    same 42 MB of cell centres and 18 MB of volumes are identical in every time
    directory. ``postProcess -func writeCellCentres`` writes them wherever it is
    pointed, so a study that post-processes *every* written frame -- which is
    what ``cases/cai2012-health`` does to get an error-against-time curve --
    would otherwise carry 60 MB of duplicated geometry per frame, several
    gigabytes over a matrix, to say the same thing each time.

    ``mesh_time`` points the geometry at one directory and the fields at
    another. It is checked against the field arrays' length like any other
    read, so pointing it at a *different mesh* fails on the cell count rather
    than silently pairing one case's densities with another's cell centres.

    Cells with no molecules have ``rhoM = 0``; their velocity and temperature are
    set to zero rather than to ``0/0``. That is not the same as "the gas there is
    at rest and at 0 K" -- it means *there was no gas sampled there* -- and the
    centreline extraction weights by number density so those cells contribute
    nothing rather than dragging an average down.
    """
    case_dir = Path(case_dir)
    time = time or latest_time(case_dir)
    directory = case_dir / time
    if not directory.is_dir():
        raise PostError(f"no time directory {directory}")

    geometry_time = mesh_time or time
    geometry_dir = case_dir / geometry_time
    if not geometry_dir.is_dir():
        raise PostError(
            f"no time directory {geometry_dir} to read the cell centres from")

    centres_path = geometry_dir / "C"
    if not centres_path.is_file():
        raise PostError(
            f"no cell centres at {centres_path}. They are written by\n"
            f"    postProcess -func writeCellCentres -time {geometry_time}\n"
            f"which ./Allpost runs. Without them there is no way to say where a "
            f"sampled value is.")
    centres = read_internal_field(centres_path)

    volumes_path = geometry_dir / "V"
    if volumes_path.is_file():
        volumes = read_internal_field(volumes_path)
    else:
        # Unweighted rather than absent: a graded mesh makes this a real
        # approximation, so it is reported by the caller, not hidden.
        volumes = np.ones(centres.shape[0], dtype=np.float64)
    if volumes.size == 1:
        volumes = np.full(centres.shape[0], float(volumes[0]))

    values = {}
    for name in REQUIRED_FIELDS:
        path = directory / name
        if not path.is_file():
            raise PostError(
                f"no {name} at {path}. It is written by the fieldAverage "
                f"function object in system/controlDict; a run that ended before "
                f"dsmc.average_start_s produces none.")
        array = read_internal_field(path)
        if array.shape[0] == 1 and centres.shape[0] > 1:
            array = np.broadcast_to(array, (centres.shape[0],) + array.shape[1:])
        elif array.shape[0] != centres.shape[0]:
            # Only reachable when the geometry came from somewhere else, and
            # the one failure mode that matters: a cell count that disagrees
            # means these densities belong to a different mesh, and pairing
            # them positionally would put every sample in the wrong place
            # while producing a plot that looks entirely ordinary.
            raise PostError(
                f"{path} has {array.shape[0]:,} values but the cell centres in "
                f"{centres_path} describe {centres.shape[0]:,} cells. These are "
                f"not the same mesh.")
        values[name] = np.asarray(array, dtype=np.float64)

    rho_m = values["rhoMMean"]
    momentum = values["momentumMean"]
    kinetic = values["linearKEMean"]

    occupied = rho_m > 0.0
    velocity = np.zeros((centres.shape[0], 3), dtype=np.float64)
    velocity[occupied] = momentum[occupied] / rho_m[occupied, None]

    mean_square = np.zeros(centres.shape[0], dtype=np.float64)
    mean_square[occupied] = 2.0 * kinetic[occupied] / rho_m[occupied]

    gas_constant = 1.380649e-23 / float(species_mass_kg)
    temperature = np.zeros(centres.shape[0], dtype=np.float64)
    temperature[occupied] = np.maximum(
        0.0,
        (mean_square[occupied] - (velocity[occupied] ** 2).sum(axis=1))
        / (3.0 * gas_constant))

    return SampledField(
        time=time,
        centres=np.asarray(centres, dtype=np.float64),
        volumes=np.asarray(volumes, dtype=np.float64),
        number_density=values["rhoNMean"],
        velocity=velocity,
        temperature=temperature,
    )


# --------------------------------------------------------------------------- #
# A / B / C -- the centreline
# --------------------------------------------------------------------------- #

def centerline(sampled: SampledField, cfg, geom, exit_state) -> dict:
    """Extract the DSMC centreline, normalised, next to the analytical one.

    Args:
        sampled: the read case.
        cfg: the case config.
        geom: the geometry.
        exit_state: the exit state, for ``n0``, ``beta0`` and ``T0``.

    Returns:
        A dict of equal-length arrays::

            x_over_D  n_over_n0  U_sqrt_beta0  T_over_T0  n_cells
            n_over_n0_analytical       U_sqrt_beta0_analytical
            T_over_T0_analytical
            n_over_n0_analytical_axis  U_sqrt_beta0_analytical_axis
            T_over_T0_analytical_axis

        Empty bins are dropped, so the arrays may be shorter than
        ``post.centerline_points``.

    Cells within ``post.centerline_radius_over_D`` of the axis are binned in
    ``x``. The density is a **volume-weighted** mean -- total molecules over
    total volume -- and the velocity and temperature are **molecule-weighted**,
    which is what the sampling supports: a cell with no molecules carries no
    information about the velocity there, and on a graded mesh an unweighted mean
    lets one huge outer cell outvote a hundred core cells.

    Averaging over a tube, not a line
    ---------------------------------
    A DSMC "centreline" is unavoidably an average over a tube of finite radius,
    and the plume profile is peaked on the axis, so that average sits **below**
    the on-axis analytical value. The bias is not small: at a tube radius of
    0.25 D it reaches **6.8%** in the density at ``X/D = 1``, which is the same
    size as the collisional departure this study exists to measure, and it looks
    exactly like physics.

    So the analytical solution is averaged over **the same cells with the same
    weights**, using :func:`plumetools.cai2012.analytical.moments`. That is what
    ``*_analytical`` holds, and it is what :func:`centerline_metrics` scores
    against. The on-axis curves -- the quantity Cai actually plots -- are also
    returned, as ``*_analytical_axis``, for the figures.
    """
    d = geom.diameter_m
    x = sampled.centres[:, 0]
    radial = np.hypot(sampled.centres[:, 1], sampled.centres[:, 2])

    tube = radial <= float(cfg.post.centerline_radius_over_D) * d
    x_max = float(cfg.post.centerline_x_over_D_max) * d
    inside = tube & (x >= 0.0) & (x <= x_max)

    n_bins = int(cfg.post.centerline_points)
    edges = np.linspace(0.0, x_max, n_bins + 1)
    index = np.clip(np.digitize(x[inside], edges) - 1, 0, n_bins - 1)

    density = sampled.number_density[inside]
    weight = sampled.volumes[inside] * density
    axial = sampled.velocity[inside, 0]
    temperature = sampled.temperature[inside]
    volume = sampled.volumes[inside]

    counts = np.bincount(index, minlength=n_bins)
    volume_sum = np.bincount(index, weights=volume, minlength=n_bins)
    molecules = np.bincount(index, weights=weight, minlength=n_bins)
    axial_sum = np.bincount(index, weights=weight * axial, minlength=n_bins)
    temperature_sum = np.bincount(index, weights=weight * temperature,
                                  minlength=n_bins)

    valid = (counts > 0) & (volume_sum > 0.0) & (molecules > 0.0)
    centres = 0.5 * (edges[:-1] + edges[1:])[valid]

    n_over_n0 = (molecules[valid] / volume_sum[valid]
                 / exit_state.number_density_per_m3)
    u_norm = (axial_sum[valid] / molecules[valid]) * math.sqrt(exit_state.beta0)
    t_norm = (temperature_sum[valid] / molecules[valid]) / exit_state.T0_K

    # The analytical solution through the SAME estimator: evaluated at the same
    # cell centres, binned with the same weights. Without this the tube-averaging
    # bias -- 6.8% in the density at X/D = 1 for a 0.25 D tube -- is charged to
    # the DSMC as if it were a physical departure.
    tube_x = x[inside]
    positive = tube_x > 0.0
    exact_cells = analytical.moments(
        np.where(positive, tube_x, 1.0), sampled.centres[inside, 1],
        sampled.centres[inside, 2], geom.radius_m, exit_state.speed_ratio)

    # A cell straddling the exit plane has no analytical value (the disk integral
    # is singular there); fall back to the closed-form axis value, which is the
    # exact limit and differs from the tube average by nothing measurable that
    # close to a uniform exit.
    exact_n = np.where(
        positive, exact_cells["n_over_n0"],
        analytical.centerline_density_ratio(
            np.zeros_like(tube_x), geom.radius_m, exit_state.speed_ratio))
    exact_u = np.where(
        positive, exact_cells["U_sqrt_beta0"][:, 0],
        analytical.centerline_axial_speed_ratio(
            np.zeros_like(tube_x), geom.radius_m, exit_state.speed_ratio))
    exact_t = np.where(
        positive, exact_cells["T_over_T0"],
        analytical.centerline_temperature_ratio(
            np.zeros_like(tube_x), geom.radius_m, exit_state.speed_ratio))

    exact_weight = volume * exact_n
    exact_volume_sum = np.bincount(index, weights=volume * exact_n,
                                   minlength=n_bins)
    exact_molecules = np.bincount(index, weights=exact_weight, minlength=n_bins)
    exact_axial = np.bincount(index, weights=exact_weight * exact_u,
                              minlength=n_bins)
    exact_temperature = np.bincount(index, weights=exact_weight * exact_t,
                                    minlength=n_bins)

    axis = analytical.centerline_profile(
        centres, geom.radius_m, exit_state.speed_ratio)

    return {
        "x_over_D": centres / d,
        "n_over_n0": n_over_n0,
        "U_sqrt_beta0": u_norm,
        "T_over_T0": t_norm,
        "n_cells": counts[valid].astype(np.float64),
        "n_over_n0_analytical": exact_volume_sum[valid] / volume_sum[valid],
        "U_sqrt_beta0_analytical": exact_axial[valid] / exact_molecules[valid],
        "T_over_T0_analytical": exact_temperature[valid] / exact_molecules[valid],
        "n_over_n0_analytical_axis": axis["n_over_n0"],
        "U_sqrt_beta0_analytical_axis": axis["U_sqrt_beta0"],
        "T_over_T0_analytical_axis": axis["T_over_T0"],
    }


# --------------------------------------------------------------------------- #
# D -- the density plane
# --------------------------------------------------------------------------- #

def density_plane(sampled: SampledField, cfg, geom, exit_state) -> dict:
    """Extract ``n/n0`` in a centre plane, as scattered points.

    Args:
        sampled: the read case.
        cfg: the case config.
        geom: the geometry.
        exit_state: the exit state, for ``n0``.

    Returns:
        ``{"x_over_D", "lateral_over_D", "n_over_n0", "plane", "lateral_axis"}``.

    The two cell layers straddling the plane are taken, not an interpolation onto
    a regular grid. A graded Cartesian mesh has no natural grid to interpolate
    onto, and a resampled field would hide where the mesh actually is -- which,
    on a contour of a quantity spanning three decades, is the thing most worth
    being able to see. ``matplotlib.tricontour`` draws it directly.
    """
    d = geom.diameter_m
    x = sampled.centres[:, 0]
    y = sampled.centres[:, 1]
    z = sampled.centres[:, 2]

    if cfg.post.plane == "xz":
        normal, lateral, lateral_axis = y, z, "Z"
    else:
        normal, lateral, lateral_axis = z, y, "Y"

    # The layer thickness adapts to the mesh: take everything within one core
    # cell of the plane, found from the smallest |normal| present.
    finest = float(np.min(np.abs(normal))) if normal.size else 0.0
    tolerance = max(finest * 1.5, 1.0e-12)

    half = float(cfg.post.plane_half_over_D) * d
    x_max = float(cfg.post.centerline_x_over_D_max) * d
    keep = ((np.abs(normal) <= tolerance) & (x >= 0.0) & (x <= x_max)
            & (np.abs(lateral) <= half))

    return {
        "x_over_D": x[keep] / d,
        "lateral_over_D": lateral[keep] / d,
        "n_over_n0": sampled.number_density[keep] / exit_state.number_density_per_m3,
        "plane": cfg.post.plane,
        "lateral_axis": lateral_axis,
    }


def analytical_plane(cfg, geom, exit_state) -> dict:
    """The analytical ``n/n0`` on a regular grid in the same centre plane.

    Drawn under the DSMC contours for comparison. Regular here because the
    analytical solution has no mesh, so there is nothing to be faithful to.
    ``x`` starts at half a grid step: the disk integral is singular in the exit
    plane itself.
    """
    d = geom.diameter_m
    x_max = float(cfg.post.centerline_x_over_D_max) * d
    half = float(cfg.post.plane_half_over_D) * d

    nx = int(cfg.post.plane_points_x)
    # An odd lateral count, so the axis is a grid line and the two halves are
    # exact mirrors of each other by construction rather than by floating-point
    # luck.
    n_half = max(2, (int(cfg.post.plane_points_lateral) + 1) // 2)
    x = np.linspace(0.0, x_max, nx + 1)[1:] - 0.5 * x_max / nx
    half_lateral = np.linspace(0.0, half, n_half)
    lateral = np.concatenate([-half_lateral[:0:-1], half_lateral])

    # Only the non-negative half is integrated; the rest is mirrored. The
    # solution is exactly axisymmetric (asserted in the tests), and the disk
    # quadrature refines itself near the exit plane, so halving the point count
    # halves what is otherwise the slowest part of post-processing.
    xh, lh = np.meshgrid(x, half_lateral, indexing="ij")
    zeros = np.zeros_like(xh)
    if cfg.post.plane == "xz":
        half_ratio = analytical.density_ratio(xh, zeros, lh, geom.radius_m,
                                              exit_state.speed_ratio)
    else:
        half_ratio = analytical.density_ratio(xh, lh, zeros, geom.radius_m,
                                              exit_state.speed_ratio)

    ratio = np.concatenate([half_ratio[:, :0:-1], half_ratio], axis=1)
    xx, ll = np.meshgrid(x, lateral, indexing="ij")

    return {
        "x_over_D": xx / d,
        "lateral_over_D": ll / d,
        "n_over_n0": ratio,
        "plane": cfg.post.plane,
    }


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #

def error_metrics(measured, exact, *, name: str) -> dict:
    """Maximum, mean and RMS **relative** error, plus where the worst one is.

    Args:
        measured: the DSMC values.
        exact: the analytical values at the same stations.
        name: what is being compared, for the keys.

    Returns:
        ``{f"{name}_max_rel_error", ..._mean_rel_error, ..._rms_rel_error,
        ..._max_rel_error_at_index, ..._n_points}``, all fractions not percents.

    Points where the analytical value is zero are dropped rather than producing
    an infinite relative error. Far downstream the density is small but never
    zero, so in practice nothing is dropped -- the guard is for the temperature,
    which approaches a small constant, and for a truncated run.
    """
    measured = np.asarray(measured, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    if measured.shape != exact.shape:
        raise ValueError(
            f"{name}: {measured.shape} measured values against {exact.shape} "
            f"analytical ones")

    usable = np.isfinite(measured) & np.isfinite(exact) & (np.abs(exact) > 0.0)
    if not np.any(usable):
        return {
            f"{name}_max_rel_error": float("nan"),
            f"{name}_mean_rel_error": float("nan"),
            f"{name}_rms_rel_error": float("nan"),
            f"{name}_max_rel_error_at_index": -1,
            f"{name}_n_points": 0,
        }

    relative = np.abs(measured[usable] - exact[usable]) / np.abs(exact[usable])
    worst = int(np.argmax(relative))
    return {
        f"{name}_max_rel_error": float(relative.max()),
        f"{name}_mean_rel_error": float(relative.mean()),
        f"{name}_rms_rel_error": float(math.sqrt(float((relative ** 2).mean()))),
        f"{name}_max_rel_error_at_index": int(np.flatnonzero(usable)[worst]),
        f"{name}_n_points": int(relative.size),
    }


def centerline_metrics(profile: dict) -> dict:
    """Error metrics for all three centreline quantities."""
    metrics = {}
    for key, label in (("n_over_n0", "density"),
                       ("U_sqrt_beta0", "velocity"),
                       ("T_over_T0", "temperature")):
        metrics.update(error_metrics(
            profile[key], profile[f"{key}_analytical"], name=label))
    worst = metrics.get("density_max_rel_error_at_index", -1)
    if worst >= 0:
        metrics["density_max_rel_error_at_x_over_D"] = float(
            profile["x_over_D"][worst])
    return metrics


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #

def write_csv(path: Path, columns: dict) -> Path:
    """Write a dict of equal-length arrays as a CSV, in key order."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    arrays = [np.asarray(columns[n], dtype=np.float64).ravel() for n in names]
    lengths = {a.size for a in arrays}
    if len(lengths) != 1:
        raise ValueError(
            f"{path}: columns have different lengths "
            f"{ {n: a.size for n, a in zip(names, arrays)} }")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(names) + "\n")
        for row in zip(*arrays):
            f.write(",".join(f"{v:.10g}" for v in row) + "\n")
    return path


def write_metrics(path: Path, document: dict) -> Path:
    """Write ``metrics.yaml`` for one case."""
    import yaml

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# GENERATED by plumetools.cai2012.post -- do not edit.\n"
        "# Relative errors are fractions, not percents.\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


def write_results(results_dir: Path, profile: dict, plane: dict,
                  metrics: dict) -> list[Path]:
    """Write ``centerline.csv``, ``density_plane.csv`` and ``metrics.yaml``."""
    results_dir = Path(results_dir)
    written = [
        write_csv(results_dir / "centerline.csv", {
            "x_over_D": profile["x_over_D"],
            "n_over_n0": profile["n_over_n0"],
            "n_over_n0_analytical": profile["n_over_n0_analytical"],
            "n_over_n0_analytical_axis": profile["n_over_n0_analytical_axis"],
            "U_sqrt_beta0": profile["U_sqrt_beta0"],
            "U_sqrt_beta0_analytical": profile["U_sqrt_beta0_analytical"],
            "U_sqrt_beta0_analytical_axis":
                profile["U_sqrt_beta0_analytical_axis"],
            "T_over_T0": profile["T_over_T0"],
            "T_over_T0_analytical": profile["T_over_T0_analytical"],
            "T_over_T0_analytical_axis": profile["T_over_T0_analytical_axis"],
            "n_cells": profile["n_cells"],
        }),
        write_csv(results_dir / "density_plane.csv", {
            "x_over_D": plane["x_over_D"],
            f"{plane['lateral_axis'].lower()}_over_D": plane["lateral_over_D"],
            "n_over_n0": plane["n_over_n0"],
        }),
        write_metrics(results_dir / "metrics.yaml", metrics),
    ]
    return written


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #

def plot_case(results_dir: Path, profile: dict, plane: dict,
              exact_plane: dict, cfg, label: str) -> list[Path]:
    """Draw the four validation figures. Returns ``[]`` without matplotlib.

    matplotlib is an optional dependency (``pip install -e ".[plots]"``), and a
    missing plotting library must not fail a validation run whose numbers are
    already in the CSVs.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    written = []

    panels = (
        ("n_over_n0", "$n/n_0$", "centerline-density.png", True),
        ("U_sqrt_beta0", r"$U_1\sqrt{\beta_0}$", "centerline-velocity.png", False),
        ("T_over_T0", "$T/T_0$", "centerline-temperature.png", False),
    )
    for key, ylabel, filename, log in panels:
        figure, axes = plt.subplots(figsize=(6.0, 4.2))
        axes.plot(profile["x_over_D"], profile[f"{key}_analytical_axis"],
                  "-", color="black", linewidth=1.6,
                  label="Cai collisionless, on axis")
        axes.plot(profile["x_over_D"], profile[key], "o", markersize=3.0,
                  color="tab:red", label=f"dsmcFoam, {label}")
        axes.set_xlabel("$X/D$")
        axes.set_ylabel(ylabel)
        if log:
            axes.set_yscale("log")
        axes.grid(True, which="both", alpha=0.3)
        axes.legend()
        axes.set_title(f"Cai 2012 centreline -- {label}")
        figure.tight_layout()
        path = results_dir / filename
        figure.savefig(path, dpi=150)
        plt.close(figure)
        written.append(path)

    levels = sorted(float(v) for v in cfg.post.contour_levels)
    figure, axes = plt.subplots(figsize=(7.0, 5.0))
    if plane["n_over_n0"].size >= 3:
        # tricontour on the raw cell centres: no resampling, so the contour
        # shows what the mesh actually resolved.
        positive = plane["n_over_n0"] > 0.0
        axes.tricontour(plane["x_over_D"][positive],
                        plane["lateral_over_D"][positive],
                        plane["n_over_n0"][positive],
                        levels=levels, colors="tab:red", linewidths=1.4)
    axes.contour(exact_plane["x_over_D"], exact_plane["lateral_over_D"],
                 exact_plane["n_over_n0"], levels=levels, colors="black",
                 linewidths=1.0, linestyles="dashed")
    axes.set_xlabel("$X/D$")
    axes.set_ylabel(f"${plane['lateral_axis']}/D$")
    axes.set_title(f"$n/n_0$ contours {levels} -- {label}\n"
                   f"solid: dsmcFoam    dashed: Cai collisionless")
    axes.set_aspect("equal", adjustable="box")
    figure.tight_layout()
    path = results_dir / "density-contour.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    written.append(path)

    return written


# --------------------------------------------------------------------------- #
# the study table
# --------------------------------------------------------------------------- #

#: Columns of ``results/study-table.csv``, in order.
STUDY_COLUMNS = (
    "case", "Kn", "n0_per_m3", "cell_over_mfp", "n_cells", "time",
    "density_max_rel_error_percent", "density_mean_rel_error_percent",
    "density_rms_rel_error_percent", "cai_reported_max_percent",
    "velocity_max_rel_error_percent", "temperature_max_rel_error_percent",
)


def study_row(metrics: dict) -> dict:
    """One row of the combined study table, from a case's ``metrics.yaml``."""
    case = metrics.get("case", {})
    errors = metrics.get("centerline", {})
    return {
        "case": case.get("name", "?"),
        "Kn": case.get("Kn"),
        "n0_per_m3": case.get("n0_per_m3"),
        "cell_over_mfp": case.get("cell_over_mean_free_path"),
        "n_cells": case.get("n_cells"),
        "time": case.get("time"),
        "density_max_rel_error_percent":
            _percent(errors.get("density_max_rel_error")),
        "density_mean_rel_error_percent":
            _percent(errors.get("density_mean_rel_error")),
        "density_rms_rel_error_percent":
            _percent(errors.get("density_rms_rel_error")),
        "cai_reported_max_percent": case.get("cai_reported_max_percent"),
        "velocity_max_rel_error_percent":
            _percent(errors.get("velocity_max_rel_error")),
        "temperature_max_rel_error_percent":
            _percent(errors.get("temperature_max_rel_error")),
    }


def _percent(value):
    return None if value is None else 100.0 * float(value)


def write_study_table(path: Path, rows: list) -> Path:
    """Write the combined study table, one row per case."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(STUDY_COLUMNS) + "\n")
        for row in rows:
            cells = []
            for column in STUDY_COLUMNS:
                value = row.get(column)
                if value is None:
                    cells.append("")
                elif isinstance(value, float):
                    cells.append(f"{value:.6g}")
                else:
                    cells.append(str(value))
            f.write(",".join(cells) + "\n")
    return path
