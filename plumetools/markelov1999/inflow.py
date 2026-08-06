"""Evaluate the paper-faithful source-flow model over a case's inflow patch.

Reads the meshed inflow patch, derives each face's off-axis angle and radius from
its centroid, and evaluates :mod:`plumetools.markelov1999.sourceflow` there.

Deliberately **not** :mod:`plumetools.inflow`
---------------------------------------------
That module is the orchestration for the frozen legacy model. It hard-wires the
legacy conventions -- ``geometry.spherical`` with ``theta`` from ``+z``, a fixed
sphere radius substituted for each face's own ``|r|``, and the
``legacy.rhoN_T0_K`` split that evaluates density at one temperature and velocity
at another (SM-02). Threading the paper-faithful model through those switches
would mean adding a second meaning to every one of them.

Two differences of substance from the legacy orchestration:

**Each face's own radius.** The legacy code substitutes the nominal sphere radius
for every face, which is right for its case: the committed inflow surface is an
exact sphere and the centroid spread is pure faceting (SM-09). Here the inflow
patch comes out of snappyHexMesh as snapped polygons whose centroids sit up to a
cell inside the true sphere, and since density goes as ``1/r**2`` that faceting
is a real effect on the flux. Using the actual ``|r|`` makes the per-face density
consistent with the position the particle is actually injected at, which is what
the flux verification in :mod:`plumetools.markelov1999.flux` then checks.

**One temperature.** Density, velocity and temperature all use
``stagnation.T0_K``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from plumetools import geometry as facegeom
from plumetools.inflow import read_patch_geometry
from plumetools.markelov1999 import sourceflow as sf
from plumetools.markelov1999.constants import molecular_mass_kg


@dataclass(frozen=True)
class InflowResult:
    """Per-face inflow conditions, every array in patch-face order.

    Face order is the contract with the solver: ``0/boundaryU`` and friends are
    consumed in patch-face order, so it is carried explicitly rather than
    reconstructed.

    Attributes:
        labels: global face labels, as strings, in patch order.
        faces: point-label lists per face, in winding order.
        centroids: ``(n, 3)`` face centroids [m].
        normals: ``(n, 3)`` unit face normals, oriented as OpenFOAM winds them
            (out of the fluid, so on this cavity they point toward the origin).
        areas: ``(n,)`` face areas [m^2].
        radius: ``(n,)`` distance of each centroid from the orifice centre [m].
        theta: ``(n,)`` off-axis angle from ``+x`` [rad].
        rhoN: ``(n,)`` number density [1/m^3].
        rhoM: ``(n,)`` mass density [kg/m^3].
        U: ``(n, 3)`` velocity [m/s], radial with magnitude ``v_limit``.
        T: ``(n,)`` temperature [K].
        v_limit: the limiting speed [m/s].
        mass_kg: molecular mass used throughout [kg].
    """

    labels: list
    faces: list
    centroids: np.ndarray
    normals: np.ndarray
    areas: np.ndarray
    radius: np.ndarray
    theta: np.ndarray
    rhoN: np.ndarray
    rhoM: np.ndarray
    U: np.ndarray
    T: np.ndarray
    v_limit: float
    mass_kg: float

    @property
    def n_faces(self) -> int:
        return len(self.labels)

    def summary(self) -> list[str]:
        """Report lines describing what the model produced on this patch."""
        return [
            f"  {self.n_faces} inflow faces, total area "
            f"{self.areas.sum():.6e} m^2",
            f"  |r|    {self.radius.min():.6f} .. {self.radius.max():.6f} m",
            f"  theta  {np.degrees(self.theta.min()):.2f} .. "
            f"{np.degrees(self.theta.max()):.2f} deg from +x",
            f"  |U|    {np.linalg.norm(self.U, axis=1).mean():.4f} m/s "
            f"(v_limit {self.v_limit:.4f})",
            f"  T      {self.T[0]:g} K",
            f"  rhoN   {self.rhoN.min():.6e} .. {self.rhoN.max():.6e} 1/m^3 "
            f"(spread {self.rhoN.max() / max(self.rhoN.min(), 1e-300):.1f}x)",
        ]


def molecular_mass_from_config(cfg) -> float:
    """The molecular mass [kg] both the model and the solver use.

    ``dsmc.species.mass_kg`` when set, otherwise derived from the molar mass. One
    number, whichever way it is obtained -- the analytical model and
    ``constant/dsmcProperties`` read it from here, so they cannot disagree.
    """
    explicit = cfg.dsmc.species.mass_kg
    if explicit is not None:
        return float(explicit)
    return molecular_mass_kg(cfg.dsmc.species.molar_mass_g_per_mol)


def compute_inflow(case_dir: Path, cfg) -> InflowResult:
    """Evaluate the source-flow model over the case's inflow patch.

    Args:
        case_dir: an OpenFOAM case directory containing ``constant/polyMesh``.
        cfg: the case configuration.

    Returns:
        An :class:`InflowResult`.

    Raises:
        NotImplementedError: for a ``model`` other than
            ``markelov1999_axisymmetric``. Running the legacy model here would
            silently apply the ``+z`` angle convention to this geometry.
        ValueError: if any face lies outside the plume cone, which for a
            hemispherical source would mean the mesh is not the geometry the
            model assumes.
    """
    if cfg.model != sf.MODEL_NAME:
        raise NotImplementedError(
            f"plumetools.markelov1999.inflow implements model {sf.MODEL_NAME!r}, "
            f"but this case declares {cfg.model!r}"
        )

    labels, faces, points = read_patch_geometry(case_dir, cfg.geometry.patch)

    centroids = facegeom.centroids(faces, points)
    normals = facegeom.normals(faces, points)
    areas = facegeom.areas(faces, points)

    radius = np.linalg.norm(centroids, axis=1)
    theta = sf.off_axis_angle(centroids, plume_axis=cfg.geometry.plume_axis)

    mass_kg = molecular_mass_from_config(cfg)
    gamma = cfg.gas.gamma
    T0 = cfg.stagnation.T0_K

    rhoM = sf.mass_density(
        theta, radius,
        gamma=gamma,
        T0_K=T0,
        mass_kg=mass_kg,
        p0_pa=cfg.stagnation.p0_pa,
        orifice_radius_m=cfg.stagnation.throat_radius_m,
        exponent_offset=cfg.angular.exponent_offset,
        quadrature_points=cfg.angular.quadrature_points,
    )
    rhoN = rhoM / mass_kg

    if not np.all(np.isfinite(rhoN)):
        raise ValueError(
            f"{np.count_nonzero(~np.isfinite(rhoN))} of {len(rhoN)} inflow faces "
            f"produced a non-finite number density. The paper-faithful angular "
            f"function cannot do this by construction, so suspect the mesh: a face "
            f"centroid at the origin would."
        )

    v_limit = sf.limiting_velocity(gamma, T0, mass_kg)
    U = sf.radial_velocity(v_limit, centroids)
    T = np.full(len(labels), float(T0), dtype=np.float64)

    return InflowResult(
        labels=[str(x) for x in labels],
        faces=faces,
        centroids=centroids,
        normals=normals,
        areas=areas,
        radius=radius,
        theta=theta,
        rhoN=np.asarray(rhoN, dtype=np.float64),
        rhoM=np.asarray(rhoM, dtype=np.float64),
        U=U,
        T=T,
        v_limit=float(v_limit),
        mass_kg=float(mass_kg),
    )


def faces_outside_the_cone(inflow: InflowResult, gamma: float) -> np.ndarray:
    """Indices of inflow faces at or beyond the limiting angle.

    Those faces have exactly zero density, so they inject nothing. On a
    hemispherical source with a diatomic gas the set is always empty -- the cone
    reaches 130 deg and the hemisphere only 90 -- which is why it is *checked*
    rather than assumed: a non-empty result means the inflow surface is not the
    hemisphere the configuration describes.
    """
    return np.flatnonzero(inflow.theta >= sf.limiting_angle(gamma))


def area_weighted_number_density(inflow: InflowResult) -> float:
    """Area-weighted mean number density [1/m^3].

    What stock ``FreeStream`` would have to use, since it takes one scalar per
    species. Reported by ``./Allrun`` **for comparison only**: this case family
    runs ``plumeFieldInflow``, which reads the per-face field, so the collapse is
    never actually applied. Printing it makes the size of what would be lost
    visible -- on this geometry the per-face density spans a wide range and any
    single value discards the angular structure the model exists to produce.
    """
    total = float(inflow.areas.sum())
    if total <= 0.0:
        raise ValueError("inflow patch has zero area")
    return float((inflow.rhoN * inflow.areas).sum() / total)
