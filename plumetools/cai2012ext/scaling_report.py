"""Numerical exports and paired Cartesian/source-centered scientific figures."""
from __future__ import annotations

import copy
from pathlib import Path
import time

import numpy as np

from .scaling import CaiPredictor, CoordinateMap
from . import scaling_fit as sf
from .scaling_workflow import (LABELS, _fit_kwargs, load_available, options, run_assessment,
                               synthetic_cases, write_csv, write_json)


def _mean_profile(coordinate, values, weights, edges):
    indices = np.clip(np.searchsorted(edges, coordinate, side="right")-1, 0, len(edges)-2)
    result = []
    for index in range(len(edges)-1):
        mask = indices == index
        if np.any(mask):
            result.append((np.average(coordinate[mask], weights=weights[mask]),
                           np.average(values[mask], weights=weights[mask])))
    return np.asarray(result).reshape(-1, 2)


def _plot_maps(plt, destination, case, prepared, predictions, title):
    from matplotlib.colors import LogNorm, Normalize
    points = prepared.points
    radius = np.linalg.norm(points, axis=1)
    angle = np.rad2deg(np.arctan2(points[:, 1], points[:, 0]))
    fields = {"DSMC" if not case.metadata.get("synthetic") else "SYNTHETIC reference": prepared.density,
              **{LABELS[k]: v for k, v in predictions.items()}}
    for coordinate_name, coords, labels in [
            ("cartesian", points, ("$x/D_e$", r"$\rho/D_e$")),
            ("source-centered", np.column_stack((radius, angle)),
             ("$r/D_e$", r"$\theta$ [deg]"))]:
        for residual in [False, True]:
            shown = dict(list(fields.items())[1:]) if residual else fields
            fig, axes = plt.subplots(1, len(shown), figsize=(3.0*len(shown), 3.4), squeeze=False)
            for ax, (name, values) in zip(axes[0], shown.items()):
                value = ((values-prepared.density)/np.maximum(prepared.density, 1e-10)
                         if residual else np.maximum(values, 1e-10))
                artist = ax.scatter(coords[:, 0], coords[:, 1], c=value, s=8, rasterized=True,
                                    cmap="RdBu_r" if residual else "viridis",
                                    norm=Normalize(-.5, .5) if residual else LogNorm(1e-4, 1.0))
                if not residual and len(coords) >= 3:
                    try:
                        ax.tricontour(coords[:, 0], coords[:, 1], values,
                                      levels=[.001, .01, .1], colors="black", linewidths=.6)
                    except (ValueError, RuntimeError):
                        pass
                ax.set(xlabel=labels[0], ylabel=labels[1], title=name)
                ax.grid(alpha=.2)
            fig.suptitle(title, fontsize=9)
            fig.subplots_adjust(left=.06, right=.88, bottom=.17, top=.78, wspace=.42)
            colorbar_axis = fig.add_axes([.92, .20, .009, .53])
            fig.colorbar(artist, cax=colorbar_axis,
                         label="(model-reference)/reference" if residual else "$n/n_e$")
            fig.savefig(destination / f"{coordinate_name}-{'residual' if residual else 'density'}.png", dpi=120)
            plt.close(fig)


def _plot_profiles(plt, destination, case, prepared, predictions, cfg, title):
    points = prepared.points
    r = np.linalg.norm(points, axis=1)
    theta = np.rad2deg(np.arctan2(points[:, 1], points[:, 0]))
    values = {"SYNTHETIC reference" if case.metadata.get("synthetic") else "DSMC": prepared.density,
              **{LABELS[k]: v for k, v in predictions.items()}}
    diagnostic_radii = cfg["evaluation"]["angular_radii_D"]
    fig, axes = plt.subplots(1, 1+len(diagnostic_radii), figsize=(5*(1+len(diagnostic_radii)), 3.8))
    axes = np.atleast_1d(axes)
    profiles = []
    masks = [(theta < 7.5, points[:, 0], np.linspace(0, max(points[:, 0])*1.001, 31),
              axes[0], "$x/D_e$", "Centerline tube: theta < 7.5 deg")]
    for radius, ax in zip(diagnostic_radii, axes[1:]):
        radius = float(radius)
        masks.append((np.abs(r-radius) <= .5, theta, np.linspace(0, 90, 37), ax,
                      r"$\theta$ [deg]", f"Angular shell: r/D_e = {radius:g} +/- 0.5"))
    for mask, coord, edges, ax, xlabel, subtitle in masks:
        for name, field in values.items():
            profile = _mean_profile(coord[mask], field[mask], prepared.volumes[mask], edges)
            if len(profile):
                ax.plot(profile[:, 0], profile[:, 1], ".-", label=name, markersize=3)
                profiles.extend(dict(region=subtitle, model=name, coordinate=p[0], density=p[1]) for p in profile)
        ax.set(xlabel=xlabel, ylabel="$n/n_e$", yscale="log", title=subtitle)
        ax.grid(True, which="both", alpha=.3)
        if np.any(mask):
            ax.legend(fontsize=7)
        else:
            ax.text(.5, .5, "No samples in this shell", transform=ax.transAxes, ha="center")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(destination / "density-profiles.png", dpi=120)
    plt.close(fig)
    write_csv(destination / "profiles.csv", profiles)


def _plot_jacobians(plt, destination, prepared, jacobians, log_bound, title):
    from matplotlib.colors import LogNorm
    points = prepared.points
    polar = np.column_stack((np.linalg.norm(points, axis=1),
                             np.rad2deg(np.arctan2(points[:, 1], points[:, 0]))))
    fig, axes = plt.subplots(2, len(jacobians), figsize=(3*len(jacobians), 6), squeeze=False)
    for column, (family, values) in zip(axes.T, jacobians.items()):
        for ax, coords, labels in zip(column, [points, polar],
                                      [("$x/D_e$", r"$\rho/D_e$"),
                                       ("$r/D_e$", r"$\theta$ [deg]")]):
            artist = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=8, rasterized=True,
                                cmap="viridis", norm=LogNorm(np.exp(-3*log_bound), np.exp(3*log_bound)))
            ax.set(xlabel=labels[0], ylabel=labels[1], title=LABELS[family])
            ax.grid(alpha=.2)
    fig.subplots_adjust(left=.07, right=.87, bottom=.1, top=.85, hspace=.45, wspace=.45)
    fig.colorbar(artist, cax=fig.add_axes([.92, .15, .012, .64]), label="3D Jacobian J")
    fig.suptitle(title + "; geometric Jacobian only, not a density multiplier", fontsize=9)
    fig.savefig(destination / "mapping-jacobians.png", dpi=120)
    plt.close(fig)


def _plot_contours(plt, destination, prepared, predictions, title):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    all_rows = []
    for family, field in predictions.items():
        rows = sf.contour_diagnostics(prepared, field)
        all_rows.extend(dict(model=family, **row) for row in rows)
        for q, ax in zip([.1, .01, .001], axes):
            subset = [row for row in rows if row["level"] == q and
                      row.get("model_radius") and row.get("dsmc_radius")]
            ax.plot([row["angle_deg"] for row in subset],
                    [row["dsmc_radius"]/row["model_radius"] for row in subset], ".-",
                    label=LABELS[family], markersize=3)
            ax.axhline(1, color="grey", linewidth=.7)
            ax.set(xlabel=r"$\theta$ [deg]", title=f"n/n_e = {q:g}")
            ax.grid(alpha=.3)
    axes[0].set_ylabel(r"$\xi_{reference}/\xi_{model}$")
    for q, ax in zip([.1, .01, .001], axes):
        if not any(row["level"] == q and row.get("model_radius") and row.get("dsmc_radius")
                   for row in all_rows):
            ax.text(.5, .5, "No valid crossings", transform=ax.transAxes, ha="center", fontsize=8)
    axes[-1].legend(fontsize=7)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(destination / "contour-displacement.png", dpi=120)
    plt.close(fig)
    write_csv(destination / "contour-displacement.csv", all_rows)
    # Ray-estimator overlays preserve absent crossings as gaps in each line.
    fig, axes = plt.subplots(2, 3, figsize=(11, 6))
    reference_rows = sf.contour_diagnostics(prepared, prepared.baseline)
    series = [("SYNTHETIC reference" if prepared.case.metadata.get("synthetic") else "DSMC",
               reference_rows, "dsmc_radius")]
    series.extend((LABELS[family], sf.contour_diagnostics(prepared, field), "model_radius")
                  for family, field in predictions.items())
    for q, column in zip([.1, .01, .001], axes.T):
        for name, rows, key in series:
            selected = [row for row in rows if row["level"] == q]
            angles = np.deg2rad([row["angle_deg"] for row in selected])
            radii = np.array([np.nan if row[key] is None else row[key] for row in selected])
            column[0].plot(radii*np.cos(angles), radii*np.sin(angles), ".-", label=name, markersize=3)
            column[1].plot(radii, np.rad2deg(angles), ".-", label=name, markersize=3)
        column[0].set(xlabel="$x/D_e$", ylabel=r"$\rho/D_e$", title=f"n/n_e = {q:g}")
        column[1].set(xlabel="$r/D_e$", ylabel=r"$\theta$ [deg]")
        for ax in column:
            ax.grid(alpha=.3)
            if not any(np.any(np.isfinite(line.get_xdata())) for line in ax.lines):
                ax.text(.5, .5, "No valid crossings", transform=ax.transAxes, ha="center", fontsize=8)
    axes[0, -1].legend(fontsize=6)
    fig.suptitle(title + "; fixed ray estimators; missing crossings remain gaps", fontsize=9)
    fig.tight_layout()
    fig.savefig(destination / "contour-overlays.png", dpi=120)
    plt.close(fig)


def render_report(output, report, *, cases=None):
    """Report uses stored predictors; no coefficient refitting occurs here."""
    from .figures import _pyplot
    output = Path(output)
    cfg = report["config"]
    smoke = bool(report["synthetic"])
    if cases is None:
        cases = synthetic_cases(cfg) if smoke else load_available(report["inventory"], cfg)[0]
    plt = _pyplot()
    began = time.perf_counter()
    for case in cases:
        prepared = sf.prepare_case(case, floor=float(cfg["fit"]["floor"]),
                                   angle_edges=np.deg2rad(cfg["fit"]["angle_edges_deg"]),
                                   radial_bins=_fit_kwargs(cfg, smoke)["radial_bins"],
                                   evaluator_options=options(cfg, smoke))
        if not len(prepared.points):
            write_json(output / "figures" / case.name / "empty-coverage.json", prepared.coverage)
            continue
        for subset in cfg["fit"]["subsets"]:
            if subset == "suitable" and not case.suitable:
                continue
            for track in cfg["fit"]["tracks"]:
                records = [f for f in report["fits"] if f.get("subset") == subset and
                           f.get("track") == track and f.get("predictor") and
                           f.get("stratum") == str(case.metadata.get("stratum", "synthetic"))]
                predictions = {}
                jacobians = {}
                destination = output / "figures" / case.name / subset / track
                destination.mkdir(parents=True, exist_ok=True)
                for record in records:
                    predictor = CaiPredictor.load(output / record["predictor"])
                    family = record["family"]
                    predictions[family] = predictor.density(prepared.points, case.kn)
                    temperature = (predictor.temperature(prepared.points, case.kn)
                                   if cfg["evaluation"]["temperature_transfer"]
                                   else np.full(len(prepared.points), np.nan))
                    jacobian = predictor.mapping.jacobian(prepared.points, case.kn)
                    jacobians[family] = jacobian
                    mapped = predictor.mapping.map_points(prepared.points, case.kn)
                    write_csv(destination / f"{family}-field.csv", [
                        dict(X=p[0], P=p[1], radius_D=np.linalg.norm(p),
                             theta_deg=np.rad2deg(np.arctan2(p[1], p[0])), reference=n,
                             predicted=pred, T_transfer=temp, mapped_X=m[0], mapped_P=m[1],
                             jacobian=j, volume=v)
                        for p, n, pred, temp, m, j, v in zip(prepared.points, prepared.density,
                            predictions[family], temperature, mapped, jacobian, prepared.volumes)])
                    probe = np.array([[0, .5], [.01, .49], [.01, .51], [.1, .5]])
                    numerical_probe = prepared.points[np.linspace(0, max(0, len(prepared.points)-1),
                                                                   min(12, len(prepared.points)), dtype=int)]
                    write_json(destination / f"{family}-source-consistency.json", {
                        "physical_X_P": probe, "mapped_X_P": predictor.mapping.map_points(probe, case.kn),
                        "lip_disk_radius": .5, "exit_plane_maps_into_itself": True,
                        "note": "A mapped lip radius differing from 0.5 changes the emitting disk boundary.",
                        "source": predictor.mapping.source_diagnostics(case.kn),
                        "evaluator": predictor.evaluation_diagnostics(prepared.points, case.kn),
                        "map_quadrature": predictor.mapping.quadrature_convergence(numerical_probe, case.kn),
                        "disk_quadrature": predictor.quadrature_convergence(numerical_probe, case.kn),
                        "origin_jacobian": "see mapping diagnostics; an angular map need not be differentiable at origin"})
                if predictions and plt:
                    title = (f"{'SYNTHETIC / ' if smoke else ''}{case.name}; {subset}; track {track}"
                             + ("" if case.suitable else "; QUESTIONABLE REFERENCE"))
                    _plot_maps(plt, destination, case, prepared, predictions, title)
                    _plot_jacobians(plt, destination, prepared, jacobians,
                                    float(cfg["fit"]["log_bound"]), title)
                    _plot_profiles(plt, destination, case, prepared, predictions, cfg, title)
                    _plot_contours(plt, destination, prepared, predictions, title)
    _mapping_curves(output, report, plt)
    write_csv(output / "subset-comparison.csv", [
        {"key": f.get("key"), "subset": f.get("subset"), "track": f.get("track"),
         "family": f.get("family"), "mapping": f.get("mapping"),
         "L_F": f.get("evaluation", {}).get("L_F"), "L_C": f.get("evaluation", {}).get("L_C"),
         "coefficients": f.get("mapping", {}).get("coefficients"),
         "objective": f.get("diagnostics", {}).get("objective"),
         "evaluation": f.get("evaluation"), "identifiable": f.get("identifiable")}
        for f in report["fits"]])
    metric_rows = []
    for fitted in report["fits"]:
        for role, key in [("training", "evaluation"), ("holdout", "holdout_evaluation")]:
            for case_metrics in (fitted.get(key) or {}).get("cases", []):
                for region, metrics in case_metrics["regions"].items():
                    metric_rows.append({"fit": fitted["key"], "role": role,
                                        "case": case_metrics["name"], "Kn": case_metrics["kn"],
                                        "region": region, **metrics,
                                        "temperature_transfer": case_metrics.get("temperature_transfer"),
                                        "coverage": case_metrics["coverage"],
                                        "parameter_uncertainty": fitted.get("parameter_uncertainty")})
    write_csv(output / "regional-metrics.csv", metric_rows)
    write_json(output / "report-timing.json", {"evaluation_exports_and_plots_seconds": time.perf_counter()-began})


def _mapping_curves(output, report, plt):
    rows = []
    for record in report["fits"]:
        if not record.get("predictor"):
            continue
        mapping = CaiPredictor.load(output / record["predictor"]).mapping
        for kn in np.geomspace(1e-4, 100, 31):
            points = np.array([[1, 0], [1, 1], [10, 0], [10, 10]])
            mapped = mapping.map_points(points, kn)
            scales = mapping.scale_factors(kn) if mapping.family != "radial_angular" else {}
            for p, m, j in zip(points, mapped, mapping.jacobian(points, kn)):
                rows.append(dict(key=record["key"], Kn=kn, X=p[0], P=p[1],
                                 mapped_X=m[0], mapped_P=m[1], jacobian=j,
                                 scale_axial=scales.get("axial"), scale_transverse=scales.get("transverse"),
                                 radius_ratio=np.linalg.norm(m)/np.linalg.norm(p)))
    write_csv(output / "mapping-versus-kn.csv", rows)
    if not rows or not plt:
        return
    groups = sorted({(record["key"].split("/")[0], record["track"])
                     for record in report["fits"] if record.get("predictor")})
    for group_index, (stratum, track) in enumerate(groups):
        subsets = report["config"]["fit"]["subsets"]
        fig, axes = plt.subplots(len(subsets), 2, figsize=(10, 3.2*len(subsets)), squeeze=False)
        for subset, axis_row in zip(subsets, axes):
            for record in report["fits"]:
                if (not record.get("key", "").startswith(stratum+"/") or
                        record.get("track") != track or record.get("subset") != subset):
                    continue
                selected = [r for r in rows if r["key"] == record["key"] and r["X"] == 1 and r["P"] == 0]
                for ax, key in zip(axis_row, ["radius_ratio", "jacobian"]):
                    ax.plot([r["Kn"] for r in selected], [r[key] for r in selected],
                            label=LABELS[record["family"]])
                    ax.set(xscale="log", xlabel="$Kn_e$", ylabel=key.replace("_", " "), title=subset)
                    ax.grid(alpha=.3)
                    ax.legend(fontsize=7)
        fig.suptitle(f"{stratum}; track {track}; (X,P)=(1,0); J is not a density multiplier", fontsize=9)
        fig.tight_layout()
        fig.savefig(output / f"mapping-versus-kn-{stratum}-{track}.png", dpi=120)
        if group_index == 0:
            fig.savefig(output / "mapping-versus-kn.png", dpi=120)
        plt.close(fig)


def run_sensitivities(cases, cfg, output, *, smoke=False):
    """Explicit opt-in refits; held-out Kn are never used to choose settings."""
    variants = []
    settings = cfg["evaluation"].get("sensitivity", {})
    for values_key, fit_key in [("floors", "floor"), ("quadrature_orders", "quadrature_order"),
                                ("log_bounds", "log_bound"), ("regularization", "regularization"),
                                ("degrees", "degree"), ("angular_orders", "angular_order")]:
        variants.extend((f"{fit_key}-{value}", fit_key, value) for value in settings.get(values_key, []))
    for factor in settings.get("angle_bin_factors", []):
        width = np.diff(cfg["fit"]["angle_edges_deg"])[0] * factor
        variants.append((f"angle-width-{width}", "angle_edges_deg", np.linspace(0, 90, round(90/width)+1).tolist()))
    for factor in settings.get("radial_bin_factors", []):
        variants.append((f"radial-bins-{factor}", "radial_bins", max(4, int(cfg["fit"]["radial_bins"]*factor))))
    summaries = []
    for name, key, value in variants:
        changed = copy.deepcopy(cfg)
        changed["fit"][key] = value
        result = run_assessment(cases, changed, Path(output) / "sensitivity" / name, smoke=smoke)
        summaries.append(dict(variant=name, setting=key, value=value, fits=result["fits"]))
    write_json(Path(output) / "sensitivities.json", summaries)
