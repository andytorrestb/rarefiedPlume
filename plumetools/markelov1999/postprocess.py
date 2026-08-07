r"""Post-processing: mesh statistics, particle statistics, and wall pressures.

Generic across every pressure case -- nothing here is keyed to a particular
reservoir pressure or a particular mesh.

Surface pressure
----------------
There is no ready-made wall-pressure field in standard ``dsmcFoam``. What it
does provide is ``fD``, accumulated in ``DSMCParcel::hitWallPatch``::

    deltaFD = cloud.nParticle()*(preIMom - postIMom)/(deltaT*fA);
    cloud.fDBF()[wppIndex][wppLocalFace] += deltaFD;

That is the momentum the gas delivers to the face, per unit area per unit time --
a **force density**, dimensions ``[1 -1 -2 0 0 0 0]`` = kg m^-1 s^-2 = Pa. The
surface pressure is its wall-normal component::

    p = fD . n_hat

with ``n_hat`` the outward face normal, so a gas pushing on the wall gives a
positive pressure. The tangential part of ``fD`` is the shear stress and is
reported separately rather than folded in.

This is **not** the gas static pressure. Substituting ``n k T`` from the cell
adjacent to the wall would be a different quantity -- it omits the directed
momentum of the impinging plume, which at these speed ratios is most of the load.

``fD`` is reset every timestep by ``DSMCCloud::resetFields``, so a written
snapshot is a one-step momentum tally and is far too noisy to read. The
``fieldAverage`` function object in the generated ``controlDict`` produces
``fDMean``, which is what :func:`wall_pressure` reads.

Averaging windows
-----------------
Pressure is averaged over a small angular and axial window, area-weighted, never
from a single face: one face's value depends entirely on where snappyHexMesh
happened to put it, so it would change with the refinement level rather than with
the physics.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools import geometry as facegeom
from plumetools.inflow import read_patch_geometry
from plumetools.markelov1999.foamfields import read_patch_field
from plumetools.markelov1999.geometry import MarkelovGeometry, cylinder_surface_angles
from plumetools.mesh.boundary import read_boundary
from plumetools.mesh.polymesh import read_faces, read_points


# --------------------------------------------------------------------------- #
# 1. mesh statistics
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class MeshStats:
    """Mesh statistics, from ``checkMesh``'s log and from the mesh itself."""

    n_cells: int
    n_points: int
    n_faces: int
    patch_face_counts: dict
    cell_type_counts: dict
    max_non_orthogonality: float | None
    average_non_orthogonality: float | None
    max_skewness: float | None
    max_aspect_ratio: float | None
    min_cell_volume_m3: float | None
    max_cell_volume_m3: float | None
    check_mesh_ok: bool | None

    def as_dict(self) -> dict:
        return {
            "n_cells": self.n_cells,
            "n_points": self.n_points,
            "n_faces": self.n_faces,
            "patch_face_counts": self.patch_face_counts,
            "cell_type_counts": self.cell_type_counts,
            "max_non_orthogonality": self.max_non_orthogonality,
            "average_non_orthogonality": self.average_non_orthogonality,
            "max_skewness": self.max_skewness,
            "max_aspect_ratio": self.max_aspect_ratio,
            "min_cell_volume_m3": self.min_cell_volume_m3,
            "max_cell_volume_m3": self.max_cell_volume_m3,
            "check_mesh_ok": self.check_mesh_ok,
        }

    def report(self) -> list[str]:
        lines = [
            "Mesh statistics",
            f"  cells                {self.n_cells}",
            f"  points               {self.n_points}",
            f"  faces                {self.n_faces}",
        ]
        if self.cell_type_counts:
            mix = ", ".join(f"{k} {v}" for k, v in sorted(self.cell_type_counts.items()))
            lines.append(f"  cell types           {mix}")
        for label, value, fmt in (
            ("max non-orthogonality", self.max_non_orthogonality, "{:.2f}"),
            ("mean non-orthogonality", self.average_non_orthogonality, "{:.2f}"),
            ("max skewness", self.max_skewness, "{:.3f}"),
            ("max aspect ratio", self.max_aspect_ratio, "{:.3f}"),
            ("min cell volume [m^3]", self.min_cell_volume_m3, "{:.6e}"),
            ("max cell volume [m^3]", self.max_cell_volume_m3, "{:.6e}"),
        ):
            if value is not None:
                lines.append(f"  {label:<20} " + fmt.format(value))
        if self.check_mesh_ok is not None:
            lines.append(f"  checkMesh            "
                         f"{'OK' if self.check_mesh_ok else 'reported failures'}")
        lines.append("  patch faces          " + ", ".join(
            f"{k} {v}" for k, v in self.patch_face_counts.items()))
        return lines


#: A float, not "any run of digits, dots and e/E".
#:
#: checkMesh writes ``Max volume = 1.647404671e-05.  Cell volumes OK.`` -- with a
#: sentence-ending period immediately after the number. A loose character class
#: swallows it and ``float()`` then fails on ``'1.647404671e-05.'``. This pattern
#: requires the exponent to be followed by digits, so the period cannot be part
#: of the match.
_FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

_CHECKMESH_PATTERNS = {
    "n_cells": (r"^\s*cells:\s+(\d+)", int),
    "n_points": (r"^\s*points:\s+(\d+)", int),
    "n_faces": (r"^\s*faces:\s+(\d+)", int),
    "max_skewness": (rf"Max skewness\s*=\s*({_FLOAT})", float),
    "max_aspect_ratio": (rf"Max aspect ratio\s*=\s*({_FLOAT})", float),
}

_CELL_TYPES = ("hexahedra", "prisms", "wedges", "pyramids", "tet wedges",
               "tetrahedra", "polyhedra")


def parse_check_mesh_log(path: Path) -> dict:
    """Extract statistics from a ``log.checkMesh``.

    Returns an empty-ish dict rather than raising when the log is absent: the
    mesh-derived counts below do not need it, and a case that was meshed
    elsewhere should still produce a summary.
    """
    path = Path(path)
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")

    out: dict = {}
    for key, (pattern, cast) in _CHECKMESH_PATTERNS.items():
        m = re.search(pattern, text, re.M)
        if m:
            out[key] = cast(m.group(1))

    ortho = re.search(
        rf"non-orthogonality Max:\s*({_FLOAT})\s*average:\s*({_FLOAT})", text)
    if ortho:
        out["max_non_orthogonality"] = float(ortho.group(1))
        out["average_non_orthogonality"] = float(ortho.group(2))

    volumes = re.search(
        rf"Min volume\s*=\s*({_FLOAT})\.?\s*Max volume\s*=\s*({_FLOAT})", text)
    if volumes:
        out["min_cell_volume_m3"] = float(volumes.group(1))
        out["max_cell_volume_m3"] = float(volumes.group(2))

    types = {}
    for name in _CELL_TYPES:
        m = re.search(rf"^\s*{re.escape(name)}:\s+(\d+)", text, re.M)
        if m and int(m.group(1)) > 0:
            types[name] = int(m.group(1))
    if types:
        out["cell_type_counts"] = types

    out["check_mesh_ok"] = "Mesh OK." in text
    return out


def mesh_statistics(case_dir: Path) -> MeshStats:
    """Mesh statistics for a case.

    Patch face counts come from ``constant/polyMesh/boundary`` -- always
    available. Everything else prefers ``log.checkMesh``, since recomputing
    non-orthogonality and skewness here would be reimplementing ``checkMesh``
    badly; the point counts fall back to the mesh files when the log is absent.
    """
    case_dir = Path(case_dir)
    patches = read_boundary(case_dir)
    parsed = parse_check_mesh_log(case_dir / "log.checkMesh")

    n_points = parsed.get("n_points")
    n_faces = parsed.get("n_faces")
    if n_points is None:
        n_points = len(read_points(case_dir))
    if n_faces is None:
        n_faces = len(read_faces(case_dir))

    return MeshStats(
        n_cells=int(parsed.get("n_cells", 0)),
        n_points=int(n_points),
        n_faces=int(n_faces),
        patch_face_counts={name: p.n_faces for name, p in patches.items()},
        cell_type_counts=parsed.get("cell_type_counts", {}),
        max_non_orthogonality=parsed.get("max_non_orthogonality"),
        average_non_orthogonality=parsed.get("average_non_orthogonality"),
        max_skewness=parsed.get("max_skewness"),
        max_aspect_ratio=parsed.get("max_aspect_ratio"),
        min_cell_volume_m3=parsed.get("min_cell_volume_m3"),
        max_cell_volume_m3=parsed.get("max_cell_volume_m3"),
        check_mesh_ok=parsed.get("check_mesh_ok"),
    )


# --------------------------------------------------------------------------- #
# 2. particle statistics
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class ParticleStats:
    """Particle counts and timings, parsed from the solver log."""

    initial_particles: int | None
    final_particles: int | None
    mean_particles: float | None
    particles_inserted_total: int | None
    collisions_total: int | None
    runtime_s: float | None
    n_time_steps: int

    def as_dict(self) -> dict:
        return {
            "initial_particles": self.initial_particles,
            "final_particles": self.final_particles,
            "mean_particles": self.mean_particles,
            "particles_inserted_total": self.particles_inserted_total,
            "collisions_total": self.collisions_total,
            "runtime_s": self.runtime_s,
            "n_time_steps": self.n_time_steps,
        }

    def report(self) -> list[str]:
        def show(label, value, fmt="{}"):
            return (f"  {label:<24} " + (fmt.format(value) if value is not None
                                         else "not recorded"))
        return [
            "Particle statistics",
            show("time steps", self.n_time_steps),
            show("initial particles", self.initial_particles),
            show("steady/final particles", self.final_particles),
            show("mean particles", self.mean_particles, "{:.4e}"),
            show("particles inserted", self.particles_inserted_total),
            show("collisions", self.collisions_total),
            show("runtime [s]", self.runtime_s, "{:.1f}"),
        ]


def parse_solver_log(path: Path, steady_fraction: float = 0.5) -> ParticleStats:
    """Parse ``log.dsmcFoam`` for particle counts, collisions and runtime.

    Args:
        path: the solver log.
        steady_fraction: the trailing fraction of time steps treated as steady
            when averaging the particle count. The default halves the run, which
            matches the ``dsmc.average_start_s`` default.

    Returns:
        A :class:`ParticleStats`. Fields the log does not contain come back
        ``None`` rather than zero -- "the solver never printed this" and "the
        solver printed zero" are different facts, and a study table that silently
        showed 0 collisions for an unparsed log would be worse than one showing a
        gap.
    """
    path = Path(path)
    if not path.is_file():
        return ParticleStats(None, None, None, None, None, None, 0)

    text = path.read_text(encoding="utf-8", errors="replace")

    counts = [int(m) for m in re.findall(
        r"Number of dsmc particles\s*=\s*(\d+)", text)]
    inserted = [int(m) for m in re.findall(
        r"Particles inserted\s*=\s*(\d+)", text)]
    collisions = [int(m) for m in re.findall(r"Collisions\s*=\s*(\d+)", text)]
    exec_times = [float(m) for m in re.findall(
        r"ExecutionTime = ([\d.eE+-]+) s", text)]

    steady = counts[int(len(counts) * (1.0 - steady_fraction)):] if counts else []

    return ParticleStats(
        initial_particles=counts[0] if counts else None,
        final_particles=counts[-1] if counts else None,
        mean_particles=float(np.mean(steady)) if steady else None,
        particles_inserted_total=int(sum(inserted)) if inserted else None,
        collisions_total=int(sum(collisions)) if collisions else None,
        runtime_s=exec_times[-1] if exec_times else None,
        n_time_steps=len(counts),
    )


# --------------------------------------------------------------------------- #
# 3-5. wall pressure
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class PressureWindow:
    """An area-weighted surface average over one window."""

    name: str
    azimuth_deg: float
    half_angle_deg: float
    axial_half_height_m: float
    n_faces: int
    area_m2: float
    pressure_pa: float
    shear_pa: float

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "azimuth_deg": self.azimuth_deg,
            "half_angle_deg": self.half_angle_deg,
            "axial_half_height_m": self.axial_half_height_m,
            "n_faces": self.n_faces,
            "area_m2": self.area_m2,
            "pressure_pa": self.pressure_pa,
            "shear_pa": self.shear_pa,
        }


def latest_time_dir(case_dir: Path) -> Path | None:
    """The numerically largest time directory, or ``None``.

    Time ``0`` is excluded: it holds the initial fields, and an average field is
    not written there.
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
            times.append((value, entry))
    return max(times)[1] if times else None


def angular_difference_deg(a, b: float) -> np.ndarray:
    """Smallest absolute angle between azimuths [deg], wrapping at +-180.

    Needed because the windward window is centred on 180 deg, where the azimuth
    wraps: a face at -179 deg is 2 deg from the window centre, not 359.
    """
    d = (np.asarray(a, dtype=np.float64) - b + 180.0) % 360.0 - 180.0
    return np.abs(d)


def wall_pressure(case_dir: Path, cfg, geom: MarkelovGeometry, *,
                  time_dir: Path | None = None,
                  field: str | None = None) -> tuple[list, dict]:
    """Area-weighted surface pressure in each configured window.

    Args:
        case_dir: the case directory.
        cfg: the case config.
        geom: the geometry.
        time_dir: which time to read; the latest if omitted.
        field: the field to read; ``post.surface_pressure_field`` if omitted.

    Returns:
        ``(windows, extras)`` where ``windows`` is a list of
        :class:`PressureWindow` and ``extras`` carries the pressure ratio and the
        field provenance.

    Raises:
        FileNotFoundError: if there is no time directory or no such field.
        ValueError: if a window captures no faces -- which means the window is
            configured for a geometry the mesh does not have, and returning a
            NaN would let that reach the study table.
    """
    case_dir = Path(case_dir)
    field = field or cfg.post.surface_pressure_field
    time_dir = time_dir or latest_time_dir(case_dir)
    if time_dir is None:
        raise FileNotFoundError(
            f"{case_dir} has no time directory, so there is nothing to "
            f"post-process. Has the solver run?")

    patch = cfg.post.cylinder_patch
    values = read_patch_field(Path(time_dir) / field, patch)

    _labels, faces, points = read_patch_geometry(case_dir, patch)
    centres = facegeom.centroids(faces, points)
    areas = facegeom.areas(faces, points)
    normals = facegeom.normals(faces, points)

    if len(values) == 1 and len(faces) > 1:
        values = np.repeat(values, len(faces), axis=0)
    if len(values) != len(faces):
        raise ValueError(
            f"{field} on patch {patch!r} has {len(values)} values but the mesh "
            f"has {len(faces)} faces")

    # OpenFOAM winds a boundary face so its normal points OUT of the fluid, i.e.
    # into the wall. hitWallPatch adds the momentum the gas delivers to the wall,
    # which is along that same direction -- so p = fD . n is positive for a gas
    # pushing on the surface, with no sign flip.
    outward = normals * np.sign(np.einsum(
        "ij,ij->i",
        normals,
        centres - np.array([geom.cylinder_centre_x_m, 0.0, 0.0]),
    ))[:, None]

    normal_component = np.einsum("ij,ij->i", values, outward)
    tangential = values - normal_component[:, None] * outward
    shear_magnitude = np.linalg.norm(tangential, axis=1)

    psi_deg = np.degrees(cylinder_surface_angles(centres, geom))
    z = centres[:, 2] - geom.cylinder_centre_z_m

    windows = []
    for spec in cfg.post.windows:
        name = spec["name"]
        azimuth = float(spec["azimuth_deg"])
        half_angle = float(spec["half_angle_deg"])
        half_height = float(spec["axial_half_height_m"])

        selected = ((angular_difference_deg(psi_deg, azimuth) <= half_angle)
                    & (np.abs(z) <= half_height))
        n_selected = int(np.count_nonzero(selected))
        if n_selected == 0:
            raise ValueError(
                f"pressure window {name!r} (azimuth {azimuth} +- {half_angle} deg, "
                f"|z| <= {half_height} m) captured no faces on patch {patch!r}. "
                f"Widen it, or check that the cylinder is where the config says.")

        area = float(areas[selected].sum())
        windows.append(PressureWindow(
            name=name,
            azimuth_deg=azimuth,
            half_angle_deg=half_angle,
            axial_half_height_m=half_height,
            n_faces=n_selected,
            area_m2=area,
            pressure_pa=float((normal_component[selected] * areas[selected]).sum()
                              / area),
            shear_pa=float((shear_magnitude[selected] * areas[selected]).sum() / area),
        ))

    by_name = {w.name: w for w in windows}
    ratio, ratio_note = None, None
    if "windward" in by_name and "leeward" in by_name:
        windward = by_name["windward"].pressure_pa
        leeward = by_name["leeward"].pressure_pa
        if leeward > 0.0 and windward > 0.0:
            ratio = windward / leeward
        else:
            # A non-positive averaged surface pressure is not a physical result;
            # it is what an under-sampled wall looks like. A handful of
            # MaxwellianThermal reflections can leave with more outward normal
            # momentum than they arrived with, so the one-sided average goes
            # negative before enough hits accumulate. Dividing by it would put a
            # meaningless number -- possibly a large one -- into the study table.
            ratio_note = (
                f"not computed: windward {windward:.4e} Pa, leeward "
                f"{leeward:.4e} Pa. A non-positive averaged wall pressure means "
                f"too few particles have struck that window to average, not a "
                f"suction. Run longer, or start fieldAverage later.")

    extras = {
        "field": field,
        "time": Path(time_dir).name,
        "pressure_ratio": ratio,
        "pressure_ratio_note": ratio_note,
        "patch": patch,
    }
    return windows, extras


def pressure_report(windows, extras) -> list[str]:
    """Report lines for the wall-pressure results."""
    lines = [
        f"Cylinder surface pressure (from {extras['field']} at t = {extras['time']}, "
        f"patch {extras['patch']!r})",
        "  the wall-normal component of the time-averaged force density "
        "[1 -1 -2 0 0 0 0] = Pa",
        f"  {'window':<14} {'azimuth':>8} {'faces':>6} {'area [m^2]':>12} "
        f"{'p [Pa]':>13} {'shear [Pa]':>12}",
    ]
    for w in windows:
        lines.append(
            f"  {w.name:<14} {w.azimuth_deg:>7.0f}d {w.n_faces:>6} {w.area_m2:>12.6e} "
            f"{w.pressure_pa:>13.6e} {w.shear_pa:>12.6e}")
    if extras.get("pressure_ratio") is not None:
        lines.append(f"  windward / leeward pressure ratio: "
                     f"{extras['pressure_ratio']:.4f}")
    else:
        lines.append("  windward / leeward pressure ratio: "
                     + (extras.get("pressure_ratio_note")
                        or "not computed (a window is missing)"))
    return lines


# --------------------------------------------------------------------------- #
# study table
# --------------------------------------------------------------------------- #

#: Columns of the combined study table, in order.
STUDY_COLUMNS = (
    "case_name",
    "pressure_psi",
    "pressure_pa",
    "gap_in",
    "windward_pressure_pa",
    "leeward_pressure_pa",
    "pressure_ratio",
    "n_cells",
    "n_particles_mean",
    "particles_per_cell_cylinder",
    "particles_per_cell_plate",
    "runtime_s",
)


def summary_row(summary: dict) -> dict:
    """One study-table row from a case summary.

    Missing values become the empty string rather than 0 or NaN: a blank cell
    reads as "not measured", which is the truth when a case has not been run,
    while a zero would read as a result.
    """
    def get(*path, default=""):
        node = summary
        for key in path:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return default if node is None else node

    windows = {w["name"]: w for w in get("pressure", "windows", default=[])} \
        if isinstance(get("pressure", "windows", default=[]), list) else {}

    return {
        "case_name": get("case", "name"),
        "pressure_psi": get("case", "pressure_psi"),
        "pressure_pa": get("case", "pressure_pa"),
        "gap_in": get("case", "gap_in"),
        "windward_pressure_pa": windows.get("windward", {}).get("pressure_pa", ""),
        "leeward_pressure_pa": windows.get("leeward", {}).get("pressure_pa", ""),
        "pressure_ratio": get("pressure", "pressure_ratio"),
        "n_cells": get("mesh", "n_cells"),
        "n_particles_mean": get("particles", "mean_particles"),
        "particles_per_cell_cylinder": get("resolution_audit", "region_occupancy",
                                           "cylinder"),
        "particles_per_cell_plate": get("resolution_audit", "region_occupancy",
                                        "plate"),
        "runtime_s": get("particles", "runtime_s"),
    }


def write_study_table(path: Path, rows: list) -> Path:
    """Write the combined study table as CSV.

    Rows keep the order they are given, which ``AllpostCases`` makes the
    deterministic case order, so the table is reproducible byte for byte.
    """
    import csv

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(STUDY_COLUMNS))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in STUDY_COLUMNS})
    return path


def plot_study(table_path: Path, out_dir: Path) -> list[Path]:
    """Plot windward, leeward and ratio against reservoir pressure.

    Args:
        table_path: the CSV written by :func:`write_study_table`.
        out_dir: where to write the figures.

    Returns:
        The figure paths, or an empty list if matplotlib is absent.

    matplotlib is an optional dependency (``pip install -e ".[plots]"``), so its
    absence is reported and skipped rather than made a hard failure of the whole
    post-processing run.
    """
    import csv

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []

    with open(table_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    def column(name):
        out = []
        for row in rows:
            try:
                out.append(float(row[name]))
            except (TypeError, ValueError):
                out.append(float("nan"))
        return np.asarray(out)

    psi = column("pressure_psi")
    order = np.argsort(psi)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    for column_name, ylabel, filename, logy in (
        ("windward_pressure_pa", "windward pressure [Pa]",
         "windward_pressure.png", True),
        ("leeward_pressure_pa", "leeward pressure [Pa]",
         "leeward_pressure.png", True),
        ("pressure_ratio", "windward / leeward pressure ratio",
         "pressure_ratio.png", False),
    ):
        values = column(column_name)
        if np.all(np.isnan(values)):
            continue
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
        ax.plot(psi[order], values[order], "o-")
        ax.set_xlabel("reservoir pressure [psi]")
        ax.set_ylabel(ylabel)
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.set_title("AIAA 99-3455, 6 in gap  --  NOT validated against "
                     "DAC or experiment", fontsize=8)
        fig.tight_layout()
        path = out_dir / filename
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)

    return written
