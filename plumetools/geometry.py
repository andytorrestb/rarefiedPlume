"""Face geometry: centroids, normals, and spherical coordinates.

All lengths are metres and all angles radians.
"""

from __future__ import annotations

import numpy as np


def centroid(vertices: np.ndarray, legacy_triangle: bool = False) -> np.ndarray:
    """Vertex-average centroid of one face.

    Args:
        vertices: ``(n, 3)`` face vertex coordinates [m], in winding order.
        legacy_triangle: divide by the literal ``3.0`` regardless of ``n``,
            reproducing the pre-refactor behaviour. See below.

    Returns:
        ``(3,)`` centroid [m].

    The default divides by ``len(vertices)``. On a triangle that is arithmetically
    the original expression with ``3.0`` replaced by ``len()``, so it reproduces
    the golden exactly; on a quadrilateral it is 4/3 times smaller than the legacy
    result, which is the correction.

    ``legacy_triangle=True`` exists only to document finding AR-01: the archived
    wake-cylinder cases run this model on quadrilateral inflow faces, where the
    fixed ``3.0`` put every centroid at 4/3 of its true radius (measured mean
    \\|r\\| 0.3946 m against a true 0.2960 m). No active case sets it.

    Note this is the vertex average, not the area centroid. The two coincide for
    triangles and for parallelograms, and differ for general polygons -- preserved
    as-is because the golden depends on it.
    """
    v = np.asarray(vertices, dtype=np.float64)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"expected (n, 3) vertices, got {v.shape}")
    if len(v) < 3:
        raise ValueError(f"a face needs at least 3 vertices, got {len(v)}")
    return v.sum(axis=0) / (3.0 if legacy_triangle else len(v))


def centroids(faces: list[list[int]], points: np.ndarray,
              legacy_triangle: bool = False) -> np.ndarray:
    """Centroids of many faces. See :func:`centroid`.

    Args:
        faces: point-label lists, one per face.
        points: ``(n_points, 3)`` mesh coordinates [m].
        legacy_triangle: forwarded to :func:`centroid`.

    Returns:
        ``(n_faces, 3)`` centroids [m], in the order of ``faces``.
    """
    return np.array(
        [centroid(points[f], legacy_triangle=legacy_triangle) for f in faces],
        dtype=np.float64,
    )


def normal(vertices: np.ndarray) -> np.ndarray:
    """Unit normal of one planar face.

    Args:
        vertices: ``(n, 3)`` face vertex coordinates [m], in winding order.

    Returns:
        ``(3,)`` unit normal.

    Raises:
        ValueError: if the first three vertices are collinear (degenerate face).

    Computed as ``(p0 - p2) x (p0 - p1)``, the orientation settled by commit
    b22f573 ("fix cross product calc so that normal vectors point out"). Only the
    first three vertices are used, so the result is meaningful only for planar
    faces.

    This is **not** used by the source-flow model -- inflow velocity direction
    comes from the centroid's radial direction, not the face normal (finding
    GP-02). It is used to verify generated meshes: see
    :func:`plumetools.foamio.blockmesh` and the mesh checks in the test suite.
    """
    v = np.asarray(vertices, dtype=np.float64)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"expected (n, 3) vertices, got {v.shape}")
    if len(v) < 3:
        raise ValueError(f"a face needs at least 3 vertices, got {len(v)}")

    n = np.cross(v[0] - v[2], v[0] - v[1])
    magnitude = np.linalg.norm(n)
    if magnitude == 0.0:
        raise ValueError("degenerate face: first three vertices are collinear")
    return n / magnitude


def normals(faces: list[list[int]], points: np.ndarray) -> np.ndarray:
    """Unit normals of many faces. See :func:`normal`."""
    return np.array([normal(points[f]) for f in faces], dtype=np.float64)


def area(vertices: np.ndarray) -> float:
    """Area of one planar polygon [m^2].

    Args:
        vertices: ``(n, 3)`` face vertex coordinates [m], in winding order.

    Returns:
        Area [m^2], always positive.

    Triangle-fan sum from the first vertex, so it is exact for any planar polygon
    -- triangles, the quads blockMesh produces, and the higher-order faces
    snappyHexMesh leaves behind.
    """
    v = np.asarray(vertices, dtype=np.float64)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"expected (n, 3) vertices, got {v.shape}")
    if len(v) < 3:
        raise ValueError(f"a face needs at least 3 vertices, got {len(v)}")
    fan = np.cross(v[1:-1] - v[0], v[2:] - v[0])
    return float(np.linalg.norm(fan.sum(axis=0)) / 2.0)


def areas(faces: list[list[int]], points: np.ndarray) -> np.ndarray:
    """Areas of many faces [m^2]. See :func:`area`."""
    return np.array([area(points[f]) for f in faces], dtype=np.float64)


def spherical(points: np.ndarray, radius: float | None = None,
              polar_axis: str = "z") -> tuple[np.ndarray, np.ndarray]:
    """Convert Cartesian positions to the (theta, phi) the source-flow model uses.

    Args:
        points: ``(n, 3)`` positions [m].
        radius: sphere radius [m] to divide by. ``None`` uses each point's own
            ``|r|``. A fixed value reproduces the legacy hard-coded ``r = 0.5``.
        polar_axis: ``"z"`` (the only value the legacy code implements).

    Returns:
        ``(theta, phi)``, each ``(n,)`` [rad], with
        ``theta = arccos(z / radius)`` in ``[0, pi]`` and
        ``phi = arctan2(y, x)`` in ``(-pi, pi]``.

    Two conventions worth stating explicitly, because neither is written down
    anywhere in the original code:

    * **theta is measured from +z, which is not the plume axis.** For
      ``cases/3d-inflow`` the plume axis is +x -- the inflow hemisphere spans
      ``x in [0, 0.5]``, one-sided. The angular dependence compensates with an
      ``abs(theta - pi/2)`` shift, making the model a separable product of two
      angles in different planes rather than a function of the single off-axis
      angle ``arccos(x/r)`` (finding SM-04).
    * **A fixed ``radius`` is deliberate, not an approximation.** The committed
      inflow surface is a sphere of radius exactly 0.5 m; face centroids fall at
      0.4990-0.4997 m purely from faceting. Passing ``radius=None`` here would
      inject that faceting spread into the model (finding SM-09).
    """
    if polar_axis != "z":
        raise NotImplementedError(
            f"polar_axis={polar_axis!r} is not implemented; the source-flow model "
            f"as written measures theta from +z"
        )
    p = np.asarray(points, dtype=np.float64)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"expected (n, 3) points, got {p.shape}")

    r = np.full(len(p), float(radius)) if radius is not None else np.linalg.norm(p, axis=1)
    if np.any(r <= 0.0):
        raise ValueError("radius must be positive")

    ratio = p[:, 2] / r
    if np.any(np.abs(ratio) > 1.0):
        worst = np.abs(ratio).max()
        raise ValueError(
            f"z/radius = {worst:.6f} exceeds 1, so arccos is undefined; the fixed "
            f"radius is smaller than the points it is applied to"
        )
    return np.arccos(ratio), np.arctan2(p[:, 1], p[:, 0])
