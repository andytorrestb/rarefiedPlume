r"""Post-processing: the loads on the plate, and the plume along its axis.

Plate loads
-----------
Standard ``dsmcFoam`` has no ready-made wall-pressure field. What it has is
``fD`` and ``q``, accumulated face by face in ``DSMCParcel::hitWallPatch``::

    deltaFD = nParticle*(preIMom - postIMom)/(deltaT*fA);
    deltaQ  = nParticle*(preIE   - postIE  )/(deltaT*fA);

-- the momentum and the energy the gas hands the wall, per unit area per unit
time. ``fD`` is a **force density**, ``[1 -1 -2]`` = Pa; ``q`` is a heat flux,
W/m^2, positive into the wall. With :math:`\hat n` the face normal pointing out
of the fluid, into the plate:

.. math::

    p = f_D \cdot \hat n, \qquad
    \tau = f_D - p\,\hat n, \qquad
    F = \sum f_D\, A, \qquad
    \dot Q = \sum q\, A.

:math:`p` is the surface pressure -- positive for a gas pushing on the plate --
and :math:`\tau` the shear. This is **not** the gas static pressure: ``n k T``
from the cell next to the wall omits the directed momentum of the impinging
plume, which at a speed ratio of five is nearly all of the load.

Both fields are reset every step, so a written snapshot is one step's tally.
``fDMean`` and ``qMean``, from the ``fieldAverage`` function object, are what is
read here.

Two numbers put the totals in proportion: the jet's **thrust** and **power** at
the exit. A plate that stopped the whole jet and re-emitted nothing would feel
exactly the thrust; a diffuse plate re-emits, which adds to it, and a plate the
plume partly misses feels less. Neither bound is a validation -- they are what
makes a wrong sign or a factor of the symmetry multiplicity visible at a glance.

The centreline
--------------
``n``, the axial velocity and the translational temperature along the plume
axis, from the averaged moments ``dsmcFoam`` writes (``rhoNMean``,
``rhoMMean``, ``momentumMean``, ``linearKEMean``), binned over a tube of radius
``post.centreline_radius_over_D`` -- volume-weighted for the density and
mass-weighted for the rest. It needs cell centres and volumes, which
``./Allpost`` writes with OpenFOAM's own function objects.

On an octree mesh a tube of fixed radius is not enough by itself: downstream
the cells are larger than the nozzle, and a tube a quarter of a diameter wide
passes between their centres and samples nothing. Where that happens the tube
is widened to the four cells round the axis -- see :func:`centreline`.

That has a visible consequence, and it is not physics. The plume is peaked on
its axis, so cells twice the size average over more of the peak and read lower:
**the curve steps down wherever the mesh coarsens**. Measured on the R-4D case
at ``x = 3 m``, where the cells beside the axis go from 62.5 mm to 125 mm: the
density reads 0.173 of the exit value just upstream and 0.131 just downstream,
while the sixteen fine cells under the footprint of the four coarse ones average
0.138 and the mass flux through the interface is continuous to 0.7%. Nothing
jumped; the sample got wider. A fixed wide tube would remove the steps and blur
the near field instead -- it would report 0.45 at the exit itself -- so the
sharpest sample available at each station is kept, ``centreline.csv`` carries
the cell size it was taken at, and the figure marks every place it changes.

Under a symmetry plane only part of the plate is in the mesh. Per-face
quantities are unaffected; the **totals** are reported both for the modelled
part and for the whole plate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools.foamio.errors import FoamDataError
from plumetools.foamio.times import latest_time
from plumetools.foamio.volfield import read_internal_field, read_patch_field
from plumetools.nozzleplume.constants import BOLTZMANN_J_PER_K
from plumetools.nozzleplume.meshread import read_patches

#: Averaged cell fields the centreline needs, and whether each is a vector.
CENTRELINE_FIELDS = {
    "rhoNMean": False,
    "rhoMMean": False,
    "momentumMean": True,
    "linearKEMean": False,
}

#: The error this module raises: the shared reader's, under the name callers
#: of this module have always caught. :func:`latest_time`,
#: :func:`read_internal_field` and :func:`read_patch_field` are the shared
#: readers in :mod:`plumetools.foamio`.
PostError = FoamDataError


# --------------------------------------------------------------------------- #
# plate loads
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class PlateLoads:
    """Per-face loads on the plate patch, and their totals.

    Attributes:
        time: the time directory read.
        centres: ``(n, 3)`` face centres [m].
        areas: ``(n,)`` face areas [m^2].
        normals: ``(n, 3)`` unit normals, out of the fluid.
        force_density: ``(n, 3)`` ``fDMean`` [Pa].
        heat_flux: ``(n,)`` ``qMean`` [W/m^2], positive into the plate.
        multiplicity: copies of the modelled domain that make the whole.
        symmetry_planes: the axes with a symmetry plane.
    """

    time: str
    centres: np.ndarray
    areas: np.ndarray
    normals: np.ndarray
    force_density: np.ndarray
    heat_flux: np.ndarray
    multiplicity: int = 1
    symmetry_planes: tuple = ()

    @property
    def pressure(self) -> np.ndarray:
        """``(n,)`` surface pressure ``fD . n`` [Pa]."""
        return np.einsum("ij,ij->i", self.force_density, self.normals)

    @property
    def shear(self) -> np.ndarray:
        """``(n,)`` shear stress magnitude [Pa]."""
        tangential = self.force_density - self.pressure[:, None] * self.normals
        return np.linalg.norm(tangential, axis=1)

    @property
    def side(self) -> np.ndarray:
        """``(n,)`` +1 on the upstream face, -1 on the lee face, 0 on an edge."""
        nx = self.normals[:, 0]
        return np.where(nx > 0.999, 1, np.where(nx < -0.999, -1, 0))

    @property
    def modelled_force_n(self) -> np.ndarray:
        """``(3,)`` force on the part of the plate that is in the mesh [N]."""
        return (self.force_density * self.areas[:, None]).sum(axis=0)

    @property
    def force_n(self) -> np.ndarray:
        """``(3,)`` force on the whole plate [N].

        The axial component adds across mirror images. A component normal to a
        symmetry plane cancels against its image, so it is zero for the whole
        plate however large it is on the modelled half.
        """
        force = self.modelled_force_n * self.multiplicity
        for axis in self.symmetry_planes:
            force["xyz".index(axis)] = 0.0
        return force

    @property
    def heat_load_w(self) -> float:
        """Heat delivered to the whole plate [W]."""
        return float((self.heat_flux * self.areas).sum()) * self.multiplicity

    def summary(self, exit_state) -> dict:
        """The totals, with the jet's thrust and power beside them."""
        front = self.side == 1
        pressure = self.pressure
        document = {
            "time": self.time,
            "n_faces": int(self.areas.size),
            "symmetry_multiplicity": int(self.multiplicity),
            "force_on_plate_n": [float(v) for v in self.force_n],
            "force_on_modelled_part_n": [float(v) for v in self.modelled_force_n],
            "heat_load_w": self.heat_load_w,
            "jet_thrust_n": float(exit_state.thrust_n()),
            "jet_power_w": float(exit_state.power_w()),
            "axial_force_over_jet_thrust":
                float(self.force_n[0] / exit_state.thrust_n()),
            "heat_load_over_jet_power":
                float(self.heat_load_w / exit_state.power_w()),
        }
        if front.any():
            area = float(self.areas[front].sum())
            peak = int(np.argmax(np.where(front, pressure, -np.inf)))
            document["upstream_face"] = {
                "n_faces": int(front.sum()),
                "area_m2": area * self.multiplicity,
                "mean_pressure_pa": float(
                    (pressure[front] * self.areas[front]).sum() / area),
                "peak_pressure_pa": float(pressure[peak]),
                "peak_pressure_at_m": [float(v) for v in self.centres[peak]],
                "mean_shear_pa": float(
                    (self.shear[front] * self.areas[front]).sum() / area),
                "peak_heat_flux_w_per_m2": float(self.heat_flux[front].max()),
                "exit_dynamic_pressure_pa": float(exit_state.dynamic_pressure_pa),
            }
        return document

    def report(self, exit_state) -> list[str]:
        document = self.summary(exit_state)
        fx, fy, fz = document["force_on_plate_n"]
        lines = [
            f"Plate loads (fDMean and qMean at t = {self.time}, "
            f"{document['n_faces']} faces)",
            f"  force on the plate       ({fx:.5g} {fy:.5g} {fz:.5g}) N",
            f"    axial / jet thrust     {document['axial_force_over_jet_thrust']:.4f}"
            f"   (thrust {document['jet_thrust_n']:.5g} N)",
            f"  heat load                {document['heat_load_w']:.5g} W",
            f"    / jet power            {document['heat_load_over_jet_power']:.4f}"
            f"   (power {document['jet_power_w']:.5g} W)",
        ]
        if self.multiplicity > 1:
            lines.append(
                f"  (whole plate: the mesh holds 1/{self.multiplicity} of it; a "
                f"component normal to a symmetry plane cancels and is zero)")
        front = document.get("upstream_face")
        if front:
            x, y, z = front["peak_pressure_at_m"]
            lines += [
                f"  upstream face            {front['area_m2']:.6g} m^2",
                f"    mean pressure          {front['mean_pressure_pa']:.5g} Pa",
                f"    peak pressure          {front['peak_pressure_pa']:.5g} Pa"
                f"   at (y, z) = ({y:.4g}, {z:.4g}) m",
                f"    mean shear             {front['mean_shear_pa']:.5g} Pa",
                f"    peak heat flux         {front['peak_heat_flux_w_per_m2']:.5g}"
                f" W/m^2",
            ]
        return lines


def plate_loads(case_dir: Path, cfg, geom, *, time: str | None = None) -> PlateLoads:
    """Read the averaged plate loads at one time.

    Raises:
        PostError: with no plate in ``case.yaml``, no time directory, no
            ``fDMean`` -- which a run that ended before ``dsmc.average_start_s``
            does not have -- or a field whose length disagrees with the patch.
    """
    if geom.plate is None:
        raise PostError("case.yaml has no plate, so there are no plate loads")
    case_dir = Path(case_dir)
    time = time or latest_time(case_dir)
    directory = case_dir / time
    if not directory.is_dir():
        raise PostError(f"no time directory {directory}")

    name = cfg.mesh.patch_names["plate"]
    patch = read_patches(case_dir, [name])[name]

    fields = {}
    for field in ("fDMean", "qMean"):
        path = directory / field
        if not path.is_file():
            raise PostError(
                f"no {field} at {path}. It is written by the fieldAverage "
                f"function object once the run passes dsmc.average_start_s; a "
                f"run that ended before then has none.")
        values = read_patch_field(path, name)
        if values.shape[0] == 1 and patch.n_faces > 1:
            values = np.repeat(values, patch.n_faces, axis=0)
        if values.shape[0] != patch.n_faces:
            raise PostError(
                f"{path} has {values.shape[0]} values on patch {name!r} but the "
                f"mesh has {patch.n_faces} faces there. These are not the same "
                f"mesh.")
        fields[field] = values

    return PlateLoads(
        time=time,
        centres=patch.centres,
        areas=patch.areas,
        normals=patch.normals,
        force_density=fields["fDMean"],
        heat_flux=fields["qMean"],
        multiplicity=geom.symmetry_multiplicity,
        symmetry_planes=tuple(geom.symmetry_planes),
    )


# --------------------------------------------------------------------------- #
# the centreline
# --------------------------------------------------------------------------- #

def centreline_tube(centres: np.ndarray, volumes: np.ndarray | None, cfg,
                    geom) -> tuple:
    """Which cells the centreline samples, and into which station.

    Args:
        centres: ``(n, 3)`` cell centres [m].
        volumes: ``(n,)`` cell volumes [m^3], or ``None`` -- in which case the
            tube is not widened where the cells outgrow it.
        cfg, geom: the case config and geometry.

    Returns:
        ``(inside, index, edges)``: a mask over the cells, the station of each
        cell in it, and the ``post.centreline_points + 1`` station edges.

    Shared with :mod:`plumetools.nozzleplume.diagnostics`, so the diagnostic
    centreline samples the same cells at the same stations as this one.
    """
    x = centres[:, 0]
    radial = np.hypot(centres[:, 1], centres[:, 2])
    x_end = geom.flow_length_m

    # The tube, widened where the mesh is coarser than it. The plume axis is a
    # lattice line -- y = 0 and z = 0 are lattice planes -- so the four cells
    # round it have their centres 0.71 cell sizes from it and the next ring
    # 1.58. Three quarters of a cell therefore takes exactly the axis cells, and
    # nothing else, wherever the configured tube would have caught none.
    reach = np.maximum(
        float(cfg.post.centreline_radius_over_D) * geom.diameter_m,
        0.75 * np.cbrt(volumes) if volumes is not None else 0.0)
    inside = (radial <= reach) & (x >= 0.0) & (x <= x_end)

    n_bins = int(cfg.post.centreline_points)
    edges = np.linspace(0.0, x_end, n_bins + 1)
    index = np.clip(np.digitize(x[inside], edges) - 1, 0, n_bins - 1)
    return inside, index, edges


def centreline(case_dir: Path, cfg, geom, exit_state, *,
               time: str | None = None) -> dict:
    """The plume along its axis: ``n``, ``U_x`` and ``T``, normalised by the exit.

    Returns:
        A dict of equal-length arrays ``x_m``, ``x_over_D``, ``n_over_n_exit``,
        ``U_over_U_exit``, ``T_over_T_exit``, ``n_cells`` and ``cell_size_m``
        -- the last being the size of the cells each sample was taken from,
        which is what a step in the others at a refinement interface has to be
        read against. Empty bins are dropped. The run ends at the plate when
        there is one.

    Raises:
        PostError: for a missing field, or missing cell centres -- naming the
            command that writes them.

    The temperature is the **translational** one, ``(<v^2> - |u|^2) / 3R``. With
    internal modes the overall temperature also weighs the internal energy; that
    is not what this reports.
    """
    case_dir = Path(case_dir)
    time = time or latest_time(case_dir)
    directory = case_dir / time
    if not directory.is_dir():
        raise PostError(f"no time directory {directory}")

    centres_path = directory / "C"
    if not centres_path.is_file():
        raise PostError(
            f"no cell centres at {centres_path}. They are written by\n"
            f"    postProcess -func writeCellCentres -time {time}\n"
            f"which ./Allpost runs.")
    centres = read_internal_field(centres_path)
    n_cells = centres.shape[0]

    volumes_path = directory / "V"
    volumes = (read_internal_field(volumes_path) if volumes_path.is_file()
               else np.ones(n_cells))
    if volumes.size == 1:
        volumes = np.full(n_cells, float(volumes[0]))

    values = {}
    for name in CENTRELINE_FIELDS:
        array = read_internal_field(directory / name)
        if array.shape[0] == 1 and n_cells > 1:
            array = np.repeat(array, n_cells, axis=0)
        if array.shape[0] != n_cells:
            raise PostError(
                f"{directory / name} has {array.shape[0]:,} values but the cell "
                f"centres describe {n_cells:,} cells. These are not the same "
                f"mesh.")
        values[name] = array

    inside, index, edges = centreline_tube(
        centres, volumes if volumes_path.is_file() else None, cfg, geom)
    n_bins = edges.size - 1

    volume = volumes[inside]
    density = values["rhoNMean"][inside]
    rho = values["rhoMMean"][inside]
    momentum_x = values["momentumMean"][inside, 0]
    momentum_sq = (values["momentumMean"][inside] ** 2).sum(axis=1)
    kinetic = values["linearKEMean"][inside]

    def binned(weights):
        return np.bincount(index, weights=weights, minlength=n_bins)

    counts = np.bincount(index, minlength=n_bins)
    volume_sum = binned(volume)
    molecules = binned(volume * density)
    mass = binned(volume * rho)
    valid = (counts > 0) & (mass > 0.0) & (molecules > 0.0)

    # Mass-weighted bulk velocity; thermal energy = total KE - bulk KE, taken
    # cell by cell so a velocity gradient across the tube is not read as heat.
    occupied = rho > 0.0
    thermal = np.zeros_like(rho)
    thermal[occupied] = (kinetic[occupied]
                         - 0.5 * momentum_sq[occupied] / rho[occupied])
    velocity = binned(volume * momentum_x)[valid] / mass[valid]
    gas_constant = BOLTZMANN_J_PER_K / exit_state.species.mass_kg
    temperature = np.maximum(
        0.0, 2.0 * binned(volume * thermal)[valid] / (3.0 * gas_constant * mass[valid]))

    station = 0.5 * (edges[:-1] + edges[1:])[valid]
    return {
        "x_m": station,
        "x_over_D": station / geom.diameter_m,
        "n_over_n_exit": (molecules[valid] / volume_sum[valid]
                          / exit_state.number_density_per_m3),
        "U_over_U_exit": velocity / exit_state.speed_m_per_s,
        "T_over_T_exit": temperature / exit_state.temperature_K,
        "n_cells": counts[valid].astype(np.float64),
        "cell_size_m": binned(np.cbrt(volume))[valid] / counts[valid],
        "time": time,
    }


def sampling_changes(profile: dict, factor: float = 1.5) -> list:
    """Stations where the centreline's sampling cells change size.

    Args:
        profile: what :func:`centreline` returns.
        factor: how much the cell size has to change, either way, to count. An
            octree changes by two; a bin that straddles an interface holds both
            sizes and reads in between, which 1.5 lets through once, not twice.

    Returns:
        ``[x_m]``, the first station at each new cell size.

    Each is a place where a step in the profile is the sample getting wider or
    narrower, not the flow changing. See the module docstring.
    """
    sizes = np.asarray(profile.get("cell_size_m", ()), dtype=np.float64)
    marks = []
    if sizes.size:
        current = sizes[0]
        for x, size in zip(profile["x_m"][1:], sizes[1:]):
            if size >= factor * current or size <= current / factor:
                marks.append(float(x))
                current = size
    return marks


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #

def write_csv(path: Path, columns: dict) -> Path:
    """Write a dict of equal-length arrays as a CSV, in key order."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(columns)
    arrays = [np.asarray(columns[n], dtype=np.float64).ravel() for n in names]
    if len({a.size for a in arrays}) != 1:
        raise ValueError(
            f"{path}: columns have different lengths "
            f"{ {n: a.size for n, a in zip(names, arrays)} }")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(names) + "\n")
        for row in zip(*arrays):
            f.write(",".join(f"{v:.10g}" for v in row) + "\n")
    return path


def write_yaml(path: Path, document: dict, title: str) -> Path:
    """Write a generated YAML document with a provenance comment."""
    import yaml

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"# GENERATED by plumetools.nozzleplume.post -- do not edit.\n# {title}\n\n"
        + yaml.safe_dump(document, sort_keys=False, default_flow_style=False),
        encoding="utf-8", newline="\n")
    return path


def write_plate_results(results_dir: Path, loads: PlateLoads, exit_state) -> list[Path]:
    """Write ``plate_loads.csv`` (one row per face) and ``plate_summary.yaml``."""
    results_dir = Path(results_dir)
    return [
        write_csv(results_dir / "plate_loads.csv", {
            "x_m": loads.centres[:, 0],
            "y_m": loads.centres[:, 1],
            "z_m": loads.centres[:, 2],
            "area_m2": loads.areas,
            "side": loads.side,
            "pressure_pa": loads.pressure,
            "shear_pa": loads.shear,
            "heat_flux_w_per_m2": loads.heat_flux,
            "fD_x_pa": loads.force_density[:, 0],
            "fD_y_pa": loads.force_density[:, 1],
            "fD_z_pa": loads.force_density[:, 2],
        }),
        write_yaml(results_dir / "plate_summary.yaml", loads.summary(exit_state),
                   "Loads on the plate. 'side' in plate_loads.csv is +1 on the "
                   "upstream face, -1 on the lee face, 0 on an edge."),
    ]


def write_centreline_results(results_dir: Path, profile: dict) -> Path:
    """Write ``centreline.csv``."""
    return write_csv(Path(results_dir) / "centreline.csv", {
        key: profile[key] for key in (
            "x_m", "x_over_D", "n_over_n_exit", "U_over_U_exit",
            "T_over_T_exit", "n_cells", "cell_size_m")})


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #

def plot_results(results_dir: Path, loads: PlateLoads | None,
                 profile: dict | None, label: str) -> list[Path]:
    """Draw the figures. Returns ``[]`` without matplotlib.

    matplotlib is an optional dependency (``pip install -e ".[plots]"``), and a
    missing plotting library must not fail a post-processing run whose numbers
    are already in the CSVs.
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

    if loads is not None and int((loads.side == 1).sum()) >= 3:
        front = loads.side == 1
        y, z = loads.centres[front, 1], loads.centres[front, 2]
        for values, title, unit, filename in (
                (loads.pressure[front], "surface pressure", "Pa",
                 "plate-pressure.png"),
                (loads.heat_flux[front], "heat flux", "W/m$^2$",
                 "plate-heat-flux.png")):
            span_y, span_z = max(np.ptp(y), 1e-9), max(np.ptp(z), 1e-9)
            width = 8.0
            figure, axes = plt.subplots(
                figsize=(width, max(2.5, min(8.0, width * span_z / span_y)) + 1.0))
            # tripcolor on the raw face centres: no resampling, so the map shows
            # what the mesh resolved and nothing smoother.
            image = axes.tripcolor(y, z, values, shading="gouraud", cmap="viridis")
            figure.colorbar(image, ax=axes, label=f"{title} [{unit}]")
            axes.set_xlabel("y [m]")
            axes.set_ylabel("z [m]")
            axes.set_aspect("equal", adjustable="box")
            axes.set_title(f"Plate upstream face: {title}\n{label}, "
                           f"t = {loads.time} s", fontsize=10)
            figure.tight_layout()
            path = results_dir / filename
            figure.savefig(path, dpi=150)
            plt.close(figure)
            written.append(path)

    if profile is not None and profile["x_m"].size:
        marks = sampling_changes(profile)
        figure, panels = plt.subplots(3, 1, figsize=(6.5, 8.0), sharex=True)
        for axes, key, ylabel, log in (
                (panels[0], "n_over_n_exit", "$n / n_{exit}$", True),
                (panels[1], "U_over_U_exit", "$U_x / U_{exit}$", False),
                (panels[2], "T_over_T_exit", "$T_{tr} / T_{exit}$", True)):
            values = profile[key]
            shown = values > 0.0 if log else np.ones_like(values, dtype=bool)
            axes.plot(profile["x_m"][shown], values[shown], "o-", markersize=2.5,
                      linewidth=0.8, color="tab:red")
            # Where the sampling cells change size. A step in the curve at one
            # of these is the sample getting wider, not the flow changing.
            for index, x in enumerate(marks):
                axes.axvline(x, color="0.35", linestyle=":", linewidth=0.9,
                             label="sampling cells change size: a step here\n"
                                   "is the estimator, not the flow"
                             if index == 0 and axes is panels[0] else None)
            if log:
                axes.set_yscale("log")
            axes.set_ylabel(ylabel)
            axes.grid(True, which="both", alpha=0.3)
        if marks:
            panels[0].legend(fontsize=7, loc="upper right")
        panels[2].set_xlabel("x [m] from the exit plane")
        panels[0].set_title(f"Plume centreline -- {label}, t = {profile['time']} s",
                            fontsize=10)
        figure.tight_layout()
        path = results_dir / "centreline.png"
        figure.savefig(path, dpi=150)
        plt.close(figure)
        written.append(path)

    return written
