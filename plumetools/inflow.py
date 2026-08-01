"""Evaluate the source-flow model over a case's inflow patch.

This is the orchestration the per-case ``runInflow.py`` calls. It reads the mesh,
derives face geometry, evaluates :mod:`plumetools.sourceflow`, and returns the
result; writing is :func:`plumetools.foamio.fields.write_inflow_fields`.

Nothing here spawns a subprocess. The original ``readMeshStats`` shelled out to
``checkMesh`` purely to learn the point and face counts, then parsed its log --
which meant a missing OpenFOAM installation left the script in an unbounded
parse loop (findings RB-01, RB-05). Those counts come from the mesh files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools import geometry, sourceflow
from plumetools.config import CaseConfig
from plumetools.constants import AMU_TO_KG
from plumetools.mesh.boundary import read_boundary, require_patch
from plumetools.mesh.polymesh import read_faces, read_points
from plumetools.mesh.sets import read_face_set


@dataclass(frozen=True)
class InflowResult:
    """Per-face inflow conditions on the patch, all arrays in face order.

    Face order is the contract with dsmcFoam+: values are consumed in patch-face
    order, so ``labels`` is carried explicitly rather than reconstructed.

    Attributes:
        labels: global face labels, as strings, in patch order.
        faces: point-label lists per face, in winding order.
        vertices: ``(n_unique_points, 3)`` coordinates of the patch's points [m].
        centroids: ``(n_faces, 3)`` face centroids [m].
        normals: ``(n_faces, 3)`` unit face normals. Not used by the model
            (finding GP-02); retained for mesh verification.
        theta: ``(n_faces,)`` polar angle from +z [rad].
        phi: ``(n_faces,)`` azimuth from +x in the x-y plane [rad].
        areas: ``(n_faces,)`` face areas [m^2].
        rhoN: ``(n_faces,)`` number density [1/m^3].
        U: ``(n_faces, 3)`` velocity [m/s].
        T: ``(n_faces,)`` translational temperature [K].
    """

    labels: list
    faces: list
    vertices: np.ndarray
    centroids: np.ndarray
    normals: np.ndarray
    areas: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    rhoN: np.ndarray
    U: np.ndarray
    T: np.ndarray


def read_patch_geometry(case_dir: Path, patch_name: str):
    """Read the faces and points of one patch.

    Args:
        case_dir: an OpenFOAM case directory.
        patch_name: patch to read, e.g. ``"inflow"``.

    Returns:
        ``(labels, faces, points)`` where ``labels`` are global face labels in
        patch order, ``faces`` their point-label lists, and ``points`` the full
        mesh point array.

    Face order comes from ``constant/polyMesh/sets/<patch>`` when that file
    exists, matching what the original did. When it does not -- for instance
    immediately after ``blockMesh``, before ``topoSet`` has run -- the patch's
    contiguous ``[startFace, startFace + nFaces)`` range is used instead, which is
    the same ordering.
    """
    case_dir = Path(case_dir)
    patch = require_patch(read_boundary(case_dir), patch_name)
    all_faces = read_faces(case_dir)
    points = read_points(case_dir)

    set_path = case_dir / "constant" / "polyMesh" / "sets" / patch_name
    if set_path.is_file():
        labels = read_face_set(case_dir, patch_name)
        if len(labels) != patch.n_faces:
            raise ValueError(
                f"faceSet {patch_name!r} has {len(labels)} faces but the boundary "
                f"declares {patch.n_faces}; re-run topoSet"
            )
    else:
        labels = list(range(patch.start_face, patch.start_face + patch.n_faces))

    return labels, [all_faces[i] for i in labels], points


def compute_inflow(case_dir: Path, cfg: CaseConfig) -> InflowResult:
    """Evaluate the source-flow model over the case's inflow patch.

    Args:
        case_dir: an OpenFOAM case directory containing ``constant/polyMesh``.
        cfg: the case configuration.

    Returns:
        An :class:`InflowResult`.

    Raises:
        NotImplementedError: for a ``model`` or ``angular.form`` other than the
            ones extracted here.

    Two frozen inconsistencies are applied deliberately, both configurable:

    * The density is evaluated at ``legacy.rhoN_T0_K`` while velocity and
      temperature use ``stagnation.T0_K``. In the resolution-study cases those
      are 300 K and 3000 K, so the inflow's mass and energy flux disagree by a
      factor of ten (finding SM-02). Set ``legacy.rhoN_T0_K: null`` to use one
      temperature throughout.
    * ``geometry.sphere_radius_m`` is used as a fixed radius rather than each
      face's own ``|r|``. That is correct here and deliberate: the inflow surface
      is a sphere of radius exactly 0.5 m and the centroids fall inside it purely
      from faceting (finding SM-09).
    """
    if cfg.model != "source_flow":
        raise NotImplementedError(f"model {cfg.model!r} is not implemented")
    if cfg.angular.form != "legacy_separable":
        raise NotImplementedError(
            f"angular.form {cfg.angular.form!r} is not implemented; only "
            f"'legacy_separable' has been extracted"
        )

    labels, faces, points = read_patch_geometry(case_dir, cfg.geometry.patch)

    centroids = geometry.centroids(faces, points)
    normals = geometry.normals(faces, points)
    theta, phi = geometry.spherical(centroids, radius=cfg.geometry.sphere_radius_m)

    mass_kg = cfg.gas.molar_mass_g_per_mol * AMU_TO_KG
    rho_T0 = cfg.legacy.rhoN_T0_K if cfg.legacy.rhoN_T0_K is not None else cfg.stagnation.T0_K

    rhoN = sourceflow.number_density(
        theta, phi,
        gamma=cfg.gas.gamma,
        T0_K=rho_T0,
        mass_kg=mass_kg,
        molar_mass_g_per_mol=cfg.gas.molar_mass_g_per_mol,
        p0_pa=cfg.stagnation.p0_pa,
        throat_radius_m=cfg.stagnation.throat_radius_m,
        sphere_radius_m=cfg.geometry.sphere_radius_m,
        exponent_offset=cfg.angular.exponent_offset,
        normalization_gamma=cfg.angular.normalization_gamma,
        quadrature_points=cfg.angular.quadrature_points,
    )

    if cfg.legacy.clip_beyond_limiting_angle:
        theta_l = sourceflow.limiting_angle(cfg.gas.gamma)
        rhoN = np.where(np.abs(phi) > theta_l, 0.0, rhoN)

    v_limit = sourceflow.limiting_velocity(cfg.gas.gamma, cfg.stagnation.T0_K, mass_kg)
    U = sourceflow.velocity(v_limit, theta, phi)
    T = np.full(len(labels), cfg.stagnation.T0_K, dtype=np.float64)

    unique = sorted({p for f in faces for p in f})
    return InflowResult(
        labels=[str(x) for x in labels],
        faces=faces,
        vertices=points[unique],
        centroids=centroids,
        normals=normals,
        areas=geometry.areas(faces, points),
        theta=theta,
        phi=phi,
        rhoN=np.asarray(rhoN, dtype=np.float64),
        U=U,
        T=T,
    )


def area_weighted_number_density(inflow: InflowResult) -> float:
    """The single uniform number density standard dsmcFoam would use [1/m^3].

    Args:
        inflow: an :class:`InflowResult`.

    Returns:
        Area-weighted mean of the per-face number density, i.e. the value that
        preserves the total inflow particle flux if the angular profile is
        collapsed to a constant.

    OpenFOAM's own ``dsmcFoam`` has no per-face number density: its ``FreeStream``
    inflow model reads ONE scalar per species from ``constant/dsmcProperties``
    (verified in v2512, ``FreeStream.C:93``). Per-face velocity and temperature
    *are* supported, so the plume's angular *direction* survives -- but its
    angular *density* profile does not.

    This is the value to put in ``FreeStreamCoeffs.numberDensities``. It is a
    lossy collapse, not a translation: on the reference case the per-face density
    spans a factor of ~200, and replacing it with any constant discards the
    structure the source-flow model exists to produce. Use the MNF fork's
    ``dsmcFreeStreamInflowFieldPatch``, or a custom InflowBoundaryModel, to keep
    it. See docs/solver-compatibility.md.
    """
    total = float(inflow.areas.sum())
    if total <= 0.0:
        raise ValueError("inflow patch has zero area")
    return float((inflow.rhoN * inflow.areas).sum() / total)
