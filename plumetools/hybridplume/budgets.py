r"""Control-volume budgets: CFD precursor, DSMC region, and the hybrid interface.

For a control volume with every boundary normal outward::

    d/dt int rho dV   + sum_f Fm A = 0
    d/dt int rho U dV + sum_f Fp A = 0
    d/dt int rho E dV + sum_f FE A = 0

Residuals are reported **dimensionally** and **normalised by declared positive
reference scales** -- the source's mass, momentum and energy inflow -- never by
a quantity that is itself near zero (transverse momentum, a vacuum boundary).

Which fluxes
------------
* CFD mass: the solver's own face flux ``phi`` summed per patch by the
  ``surfaceFieldValue`` function objects every step. That is the numerical
  flux the continuity equation used.
* CFD momentum and energy: no numerical face flux is written by
  ``rhoCentralFoam``. They are **reconstructed** from the written boundary
  values as inviscid fluxes (``rho U Un + p n``, ``rho Un (cp T + |U|^2/2)``)
  and labelled ``reconstructed_inviscid``; viscous boundary terms are omitted.
  The reconstructed mass flux is reported beside ``phi`` so the reconstruction
  error is visible.
* DSMC: from particles (:mod:`.particles`), never from injected counts alone.

Accumulation terms are retained until stationarity is shown: the time
derivative is a least-squares slope over the last averaging window.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from plumetools.foamio.volfield import read_patch_field

from .mesh import PATCHES
from .meshgeom import MeshGeometry


def _slope(t: np.ndarray, y: np.ndarray) -> float:
    if t.size < 2:
        return float("nan")
    a = np.vstack([t - t.mean(), np.ones_like(t)]).T
    return float(np.linalg.lstsq(a, y, rcond=None)[0][0])


def reconstructed_patch_fluxes(case: Path, time: str, mesh: MeshGeometry, gas) -> dict:
    out = {}
    d = Path(case) / time
    for name, info in mesh.patches.items():
        if info.n_faces == 0:
            continue
        sl = mesh.patch_slice(name)
        s = mesh.face_areas[sl]
        a = np.linalg.norm(s, axis=1)
        n = s / a[:, None]
        p = read_patch_field(d / "p", name).reshape(-1)
        t = read_patch_field(d / "T", name).reshape(-1)
        u = read_patch_field(d / "U", name).reshape(-1, 3)
        p = np.broadcast_to(p, (info.n_faces,)) if p.size == 1 else p
        t = np.broadcast_to(t, (info.n_faces,)) if t.size == 1 else t
        u = np.broadcast_to(u, (info.n_faces, 3)) if u.shape[0] == 1 else u
        rho = p / (gas.gas_constant * t)
        un = np.einsum("ij,ij->i", u, n)
        mass = rho * un * a
        mom = (rho[:, None] * u * un[:, None] + p[:, None] * n) * a[:, None]
        energy = rho * un * (gas.cp * t + 0.5 * np.einsum("ij,ij->i", u, u)) * a
        out[name] = {"mass": float(mass.sum()), "momentum": mom.sum(0).tolist(),
                     "energy": float(energy.sum()), "area_m2": float(a.sum())}
    return out


def cfd_budgets(case: Path, time: str, mesh: MeshGeometry, source, series: dict, cfg) -> dict:
    gas = source.gas
    transit = float(mesh.points[:, 0].max()) / source.velocity
    a_src = float(np.linalg.norm(mesh.face_areas[mesh.patch_slice(PATCHES["source"])],
                                 axis=1).sum())
    ref_mass = source.mass_density * source.velocity * a_src
    ref_mom = (source.mass_density * source.velocity ** 2 + source.pressure) * a_src
    ref_energy = ref_mass * (gas.cp * source.temperature_K + 0.5 * source.velocity ** 2)
    doc = {"reference_scales": {"mass_kg_s": ref_mass, "momentum_N": ref_mom,
                                "energy_W": ref_energy,
                                "basis": "uniform source state on the staircase source area"}}

    recon = reconstructed_patch_fluxes(case, time, mesh, gas)
    doc["reconstructed_inviscid_patch_fluxes"] = recon

    dm = series.get("domainMass")
    phis = {k: series.get(k) for k in ("boundaryMassFlux", "exitPlaneMassFlux",
                                       "farFieldMassFlux")}
    if dm is not None and all(v is not None for v in phis.values()):
        t_end = dm[-1, 0]
        win = dm[:, 0] >= t_end - transit
        tw = dm[win, 0]
        mass_rate = _slope(tw, dm[win, 1])
        mom_rate = [_slope(tw, dm[win, 2 + k]) for k in range(3)]
        energy_rate = _slope(tw, dm[win, 5]) if dm.shape[1] > 5 else float("nan")

        def window_mean(arr):
            w = arr[:, 0] >= t_end - transit
            return float(arr[w, 1].mean())
        phi = {"source": window_mean(phis["boundaryMassFlux"]),
               "exitPlane": window_mean(phis["exitPlaneMassFlux"]),
               "farField": window_mean(phis["farFieldMassFlux"])}
        net = sum(phi.values())
        resid = mass_rate + net
        doc["mass"] = {
            "numerical_phi_per_patch_kg_s": phi, "net_outflow_kg_s": net,
            "accumulation_kg_s": mass_rate, "residual_kg_s": resid,
            "residual_over_reference": resid / ref_mass,
            "source_phi_over_rhoUA_staircase": -phi["source"] / ref_mass,
            "source_phi_over_rhoUA_exact_disk": -phi["source"] / (
                source.mass_density * source.velocity * source.area_m2),
            "window_s": [float(tw[0]), float(t_end)],
        }
        mom_net = np.sum([recon[p]["momentum"] for p in recon], axis=0)
        mom_res = np.asarray(mom_rate) + mom_net
        doc["momentum"] = {
            "reconstructed_net_outflow_N": mom_net.tolist(), "accumulation_N": mom_rate,
            "residual_N": mom_res.tolist(),
            "residual_over_reference": (mom_res / ref_mom).tolist(),
            "note": "reconstructed inviscid boundary fluxes; transverse components "
                    "are normalised by the AXIAL source scale, not by themselves"}
        e_net = sum(recon[p]["energy"] for p in recon)
        doc["energy"] = {
            "reconstructed_net_outflow_W": e_net, "accumulation_W": energy_rate,
            "residual_W": energy_rate + e_net,
            "residual_over_reference": (energy_rate + e_net) / ref_energy}
        m0 = dm[win, 1][0]
        doc["stationarity"] = {
            "window_transits": 1.0,
            "relative_mass_change_over_window": float((dm[win, 1][-1] - m0) / dm[win, 1][-1]),
            "relative_energy_change_over_window": float(
                (dm[win, 5][-1] - dm[win, 5][0]) / dm[win, 5][-1]) if dm.shape[1] > 5 else None,
            "accumulation_over_reference": mass_rate / ref_mass}
        doc["mass"]["reconstructed_vs_phi"] = {
            p: {"reconstructed_kg_s": recon[k]["mass"], "phi_kg_s": phi[p]}
            for p, k in (("source", PATCHES["source"]), ("exitPlane", PATCHES["exit_plane"]),
                         ("farField", PATCHES["far_field"])) if k in recon}
    else:
        doc["mass"] = {"available": False,
                       "reason": "function-object series missing (postProcessing/)"}
    return doc


def dsmc_budget(measured: dict, accumulation: dict, references: dict) -> dict:
    """Particle-based budget of the DSMC region from crossing measurements.

    ``measured``: per-patch ``{"in": {...}, "out": {...}}`` rates (mass kg/s,
    momentum N, energy W), with "out" the physical flux LEAVING the region.
    """
    net = {"mass": 0.0, "momentum": np.zeros(3), "energy": 0.0}
    for patch, d in measured.items():
        net["mass"] += d["out"]["mass"] - d["in"]["mass"]
        net["momentum"] += np.asarray(d["out"]["momentum"]) - np.asarray(d["in"]["momentum"])
        net["energy"] += d["out"]["energy"] - d["in"]["energy"]
    res_m = accumulation.get("mass", 0.0) + net["mass"]
    res_p = np.asarray(accumulation.get("momentum", np.zeros(3))) + net["momentum"]
    res_e = accumulation.get("energy", 0.0) + net["energy"]
    return {"net_outflow": {"mass": net["mass"], "momentum": net["momentum"].tolist(),
                            "energy": net["energy"]},
            "accumulation": accumulation,
            "residual": {"mass": res_m, "momentum": res_p.tolist(), "energy": res_e},
            "residual_over_reference": {
                "mass": res_m / references["mass"],
                "momentum": (res_p / references["momentum"]).tolist(),
                "energy": res_e / references["energy"]}}
