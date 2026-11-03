r"""CFD interface state -> DSMC inflow faces, face for face, with its flux error.

Canonical interface fields (backend-independent), one row per interface face,
in CFD face order, saved to ``results/interface_fields.npz``:

==========================  ======  ==========================================
``numberDensity``           1/m^3   ``rho / m``
``massDensity``             kg/m^3
``velocity``                m/s     Cartesian
``translationalTemperature`` K
``faceArea``                m^2
``faceNormal``              --      unit, CFD -> DSMC
``faceCentre``              m
``valid``                   bool    state positive and finite
``cfdFace``                 --      the CFD mesh face label (provenance)
``cfdCell``                 --      the continuum-side CFD cell
==========================  ======  ==========================================

Argon is single-species: ``n = rho/m`` and ``p = n k T`` are checked, not
assumed (:func:`consistency`).

Face correspondence
-------------------
v2512 ``subsetMesh`` writes no addressing maps (``faceMap``/``cellMap``; read
in ``subsetMesh.C``), so each DSMC ``interface`` patch face is matched to its
CFD face geometrically -- nearest face centre (``scipy.spatial.cKDTree``) --
and the match is then *verified*: the centres must coincide to round-off, the
DSMC outward area vector must equal ``-n A``, and the match must be one-to-one.
A face without a match, or with an invalid state, is refused -- nothing is
extrapolated across missing data.

Mapping error
-------------
The state is copied, not interpolated, so the only change is the face geometry
the DSMC mesh recomputes. :func:`flux_change` integrates the Maxwellian forward
fluxes on both geometries and reports the difference, so "conservative" is a
measured statement.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import halfrange
from .meshgeom import MeshGeometry
from .physics import K_BOLTZMANN

CANONICAL = ("numberDensity", "massDensity", "velocity", "translationalTemperature",
             "faceArea", "faceNormal", "faceCentre", "valid", "cfdFace", "cfdCell")


class MappingError(ValueError):
    """An interface face that cannot be supplied with a valid state."""


def canonical(surf, face_arrays: dict, gas) -> dict:
    rho = face_arrays["rho"]
    return {
        "numberDensity": rho / gas.mass_kg, "massDensity": rho,
        "velocity": face_arrays["U"], "translationalTemperature": face_arrays["T"],
        "faceArea": surf.area, "faceNormal": surf.normal, "faceCentre": surf.centre,
        "valid": face_arrays["valid"], "cfdFace": surf.faces, "cfdCell": surf.cfd_cell,
    }


def consistency(fields: dict, gas, pressure=None) -> dict:
    n = fields["numberDensity"]
    rho = fields["massDensity"]
    out = {"max_rel_rho_vs_m_n": float(np.max(np.abs(rho / (gas.mass_kg * n) - 1.0)))}
    if pressure is not None:
        p = n * K_BOLTZMANN * fields["translationalTemperature"]
        out["max_rel_p_vs_n_k_T"] = float(np.max(np.abs(p / pressure - 1.0)))
    return out


def save(path: Path, fields: dict) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **{k: np.asarray(v) for k, v in fields.items()})
    return path


def load(path: Path) -> dict:
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


def match_cells(dsmc_centres: np.ndarray, cfd_centres: np.ndarray, cfd_sizes: np.ndarray,
                *, rel_tol: float = 1e-6) -> np.ndarray:
    """CFD cell label of every DSMC cell, by centre; verified one-to-one."""
    from scipy.spatial import cKDTree
    dist, idx = cKDTree(cfd_centres).query(dsmc_centres)
    scale = cfd_sizes[idx]
    if np.any(dist > rel_tol * scale):
        raise MappingError(f"{int((dist > rel_tol * scale).sum())} DSMC cells have no "
                           "coincident CFD cell")
    if np.unique(idx).size != idx.size:
        raise MappingError("two DSMC cells matched one CFD cell")
    return idx


def match_interface(dsmc_mesh: MeshGeometry, fields: dict,
                    patch: str = "interface", *, rel_tol: float = 1e-6) -> dict:
    """For each DSMC ``patch`` face, the row of ``fields`` that supplies it.

    Raises:
        MappingError: an unmatched face, a face whose geometry disagrees, or a
            matched row with an invalid state.
    """
    from scipy.spatial import cKDTree
    sl = dsmc_mesh.patch_slice(patch)
    new_faces = np.arange(sl.start, sl.stop)
    centre = dsmc_mesh.face_centres[new_faces]
    dist, rows = cKDTree(fields["faceCentre"]).query(centre)
    rows = rows.astype(np.int64)
    if rows.size != fields["cfdFace"].size or np.unique(rows).size != rows.size:
        raise MappingError(f"interface face counts differ or repeat: DSMC {rows.size}, "
                           f"CFD {fields['cfdFace'].size}")
    scale = np.sqrt(fields["faceArea"][rows])
    dc = dist / scale
    area_vec = dsmc_mesh.face_areas[new_faces]
    expected = -fields["faceNormal"][rows] * fields["faceArea"][rows, None]
    da = np.linalg.norm(area_vec - expected, axis=1) / fields["faceArea"][rows]
    if np.max(dc) > rel_tol or np.max(da) > 1e-6:
        raise MappingError(f"interface geometry mismatch: max centre offset {np.max(dc):.3e} "
                           f"(per sqrt area), max area-vector error {np.max(da):.3e}")
    bad = ~fields["valid"][rows]
    if bad.any():
        raise MappingError(f"{int(bad.sum())} interface faces carry an invalid CFD state")
    return {"rows": rows, "dsmc_faces": new_faces, "max_centre_offset": float(np.max(dc)),
            "max_area_vector_error": float(np.max(da)),
            "dsmc_area_vectors": area_vec}


def flux_change(fields: dict, match: dict, gas) -> dict:
    """Maxwellian forward fluxes on the CFD vs the DSMC face geometry."""
    r = match["rows"]
    pre = halfrange.integrate(halfrange.incoming_flux(
        fields["massDensity"], fields["velocity"], fields["translationalTemperature"],
        fields["faceNormal"], gas.gas_constant), fields["faceArea"])
    av = match["dsmc_area_vectors"]
    a = np.linalg.norm(av, axis=1)
    n_in = -av / a[:, None]       # DSMC patch area points out of DSMC; inflow is -it
    post = halfrange.integrate(halfrange.incoming_flux(
        fields["massDensity"][r], fields["velocity"][r],
        fields["translationalTemperature"][r], n_in, gas.gas_constant), a)
    scale_p = abs(pre["mass"]) * np.linalg.norm(fields["velocity"], axis=1).max()
    return {
        "pre_mapping": pre, "post_mapping": post,
        "mass_rel": post["mass"] / pre["mass"] - 1.0,
        "energy_rel": post["energy"] / pre["energy"] - 1.0,
        "momentum_abs": (np.asarray(post["momentum"]) - np.asarray(pre["momentum"])).tolist(),
        "momentum_ref_scale": float(scale_p),
        "momentum_rel_to_scale": (np.abs(np.asarray(post["momentum"])
                                         - np.asarray(pre["momentum"])) / scale_p).tolist(),
    }


def inflow_from_fields(fields: dict, match: dict) -> dict:
    r = match["rows"]
    return {"n": fields["numberDensity"][r], "U": fields["velocity"][r],
            "T": fields["translationalTemperature"][r]}


def source_inflow(n_faces: int, source) -> dict:
    return {"n": np.full(n_faces, source.number_density),
            "U": np.tile(source.velocity_vector, (n_faces, 1)),
            "T": np.full(n_faces, source.temperature_K)}
