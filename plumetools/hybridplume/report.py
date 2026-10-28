"""Collect a case's step reports into one summary, and draw its figures.

Figures (PNG, ``results/``):

* ``centreline.png`` -- n/ne, T/Te, Ux/ae against x/D as three panels (never a
  second y-axis), one colour per source of data in a fixed order: DSMC
  fieldAverage, DSMC parcel snapshots, CFD precursor, collisionless reference;
* ``breakdown_slice.png`` -- Kn_GLL on the cells nearest the z = 0 plane, a
  single-hue sequential ramp on a log scale, with the retained continuum region
  outlined;
* ``interface_mach.png`` -- the area distribution of the interface normal Mach
  number, with Mn = 1 marked.

Every figure carries the case label, so a SMOKE result cannot be mistaken for a
converged one.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

#: Fixed categorical order (reference palette, light surface).
SERIES = {"dsmc_fa": "#2a78d6", "dsmc_lag": "#eb6834", "cfd": "#1baf7a", "ref": "#eda100",
          "hyb": "#e87ba4"}
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK = "#2b2b2a"
MUTED = "#8a8985"

STEPS = ("setup", "mesh_check_cfd", "cfd_report_cfd", "breakdown", "interface",
         "dsmc_setup", "dsmc_fields", "dsmc_post", "crossing_window", "reference")


def _load(res: Path, name: str):
    p = res / f"{name}.json"
    return json.loads(p.read_text()) if p.is_file() else None


def summarise(case_dir: Path) -> dict:
    res = Path(case_dir) / "results"
    docs = {k: _load(res, k) for k in STEPS}
    present = [k for k, v in docs.items() if v is not None]
    label = next((d.get("label") for d in docs.values() if d and d.get("label")), None)
    head = {"label": label, "smoke": bool(label and label.upper().startswith("SMOKE"))}
    if docs["setup"]:
        s = docs["setup"]
        head["source"] = {k: s["source"]["derived"][k] for k in (
            "number_density_per_m3", "velocity_m_s", "mach", "reynolds_D")}
        head["source_area_error"] = s["source_discretisation"]["relative_area_error"]
        head["mesh_cells"] = s["mesh"]["n_cells"]
    if docs["cfd_report_cfd"] and "residual_over_reference" in (docs["cfd_report_cfd"].get("mass") or {}):
        c = docs["cfd_report_cfd"]
        head["cfd_mass_residual_over_reference"] = c["mass"]["residual_over_reference"]
        head["cfd_exitplane_inflow_over_source"] = (
            c["mass"]["numerical_phi_per_patch_kg_s"]["exitPlane"]
            / c["mass"]["numerical_phi_per_patch_kg_s"]["source"])
        head["cfd_stationarity"] = c.get("stationarity")
    if docs["interface"]:
        i = docs["interface"]
        head["interface"] = {k: i.get(k) for k in (
            "mode", "status", "reasons", "return_mass_fraction", "interface_faces",
            "outer_boundary_contact_fraction")}
    if docs["dsmc_post"]:
        p = docs["dsmc_post"]
        head["dsmc"] = {k: p.get(k) for k in ("injection", "return_flux", "stationarity",
                                              "transverse_momentum_over_axial")}
        head["centreline_comparison"] = p.get("centreline_comparison")
    if docs["crossing_window"]:
        head["crossings"] = docs["crossing_window"].get("exit_rates")
    return {"label": label, "stages_present": present, "headline": head, "steps": docs}


def _style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.grid(True, color="#e6e5e0", linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED)


def plots(case_dir: Path) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, LogNorm

    res = Path(case_dir) / "results"
    label = (summarise(case_dir)["label"]) or ""
    out = []

    cl = res / "centreline.csv"
    if cl.is_file():
        data = np.genfromtxt(cl, delimiter=",", names=True)
        x = data["x_over_D"]
        fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
        panels = (("n", "n / n$_e$", "{}_over_ne", True), ("T", "T / T$_e$", "{}_over_Te", False),
                  ("U", "U$_x$ / a$_e$", "{}_over_ae_x", False))
        for ax, (q, ylabel, pat, logy) in zip(axes, panels):
            for key, name, colour, style in (
                    (f"{q}_fa", "DSMC (fieldAverage)", SERIES["dsmc_fa"], "-"),
                    (f"{q}_lag", "DSMC (parcel snapshots)", SERIES["dsmc_lag"], "--"),
                    (f"{q}_cfd", "CFD precursor", SERIES["cfd"], "-"),
                    (f"{q}_ref", "collisionless reference", SERIES["ref"], ":"),
                    (f"{q}_hyb", "hybrid (CFD retained + DSMC)", SERIES["hyb"], "-")):
                col = pat.format(key)
                if col in data.dtype.names and np.isfinite(data[col]).any():
                    ax.plot(x, data[col], style, color=colour, linewidth=2, label=name)
            if logy:
                ax.set_yscale("log")
            _style(ax, "x / D (axis tube)", ylabel)
        axes[0].legend(frameon=False, fontsize=8)
        fig.suptitle(label, color=INK, fontsize=10)
        path = res / "centreline.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        out.append(path)

    bdz = res / "breakdown.npz"
    cfd = Path(case_dir) / "cfd"
    if bdz.is_file() and (cfd / "constant" / "polyMesh" / "owner").is_file():
        from .meshgeom import read_mesh
        mesh = read_mesh(cfd)
        with np.load(bdz) as z:
            kn = z["Kn_GLL"]
            retained = z["retained"]
        c = mesh.cell_centres
        zmin = np.min(np.abs(c[:, 2]))
        sel = np.isclose(np.abs(c[:, 2]), zmin) & (c[:, 2] > 0)
        size = mesh.cell_size()[sel]
        fig, ax = plt.subplots(figsize=(7, 5.6), constrained_layout=True)
        vals = np.where(np.isfinite(kn[sel]) & (kn[sel] > 0), kn[sel], np.nan)
        lo = max(np.nanmin(vals), 1e-4)
        sc = ax.scatter(c[sel, 0], c[sel, 1], c=vals, s=(72 * size[:, 0] / 0.02) ** 2 * 0.25,
                        marker="s", cmap=ListedColormap(SEQ_BLUE),
                        norm=LogNorm(vmin=lo, vmax=np.nanmax(vals)), linewidths=0)
        r = retained[sel]
        ax.scatter(c[sel][r, 0], c[sel][r, 1], s=3, color=SERIES["cfd"], label="retained CFD region")
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("Kn_GLL (log)", color=INK)
        ax.set_aspect("equal")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        _style(ax, "x [m]", "y [m] (cells nearest z = 0)")
        ax.set_title(label, color=INK, fontsize=9)
        path = res / "breakdown_slice.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        out.append(path)

    itf = res / "interface_fields.npz"
    if itf.is_file():
        with np.load(itf) as z:
            n = z["faceNormal"]
            u = z["velocity"]
            t = z["translationalTemperature"]
            a = z["faceArea"]
        from .config import build_gas
        from .workflow import load_case
        gas = build_gas(load_case(case_dir))
        mn = np.einsum("ij,ij->i", u, n) / gas.sound_speed(t)
        fig, ax = plt.subplots(figsize=(6, 3.6), constrained_layout=True)
        ax.hist(mn, bins=40, weights=a / a.sum(), color=SERIES["dsmc_fa"], rwidth=0.9)
        ax.axvline(1.0, color=INK, linewidth=1)
        ax.text(1.02, ax.get_ylim()[1] * 0.92, "Mn = 1", color=INK, fontsize=8)
        _style(ax, "normal Mach number Mn = U.n / c (CFD -> DSMC)", "area fraction")
        ax.set_title(label, color=INK, fontsize=9)
        path = res / "interface_mach.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        out.append(path)
    return out
