"""Study-scoped orchestration; defaults to inventory and never launches DSMC.

The serial optimizer keeps numerical libraries to one thread. ``cores`` is an
upper resource declaration, not a promise to spawn that many workers.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import time

import numpy as np
import yaml

from .scaling import CaiPredictor, CoordinateMap
from . import scaling_fit as sf
from . import scaling_inventory as si

LABELS = {"baseline": "Cai baseline", "axial": "Axial scale",
          "axial_transverse": "Axial/transverse scale",
          "radial_angular": "Radial/angular map"}


def clean(value):
    """Strict JSON: unavailable/nonfinite diagnostics become null, never NaN."""
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(document), indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def write_csv(path, rows):
    rows = [clean(row) for row in rows]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows({k: json.dumps(v) if isinstance(v, (dict, list)) else v
                         for k, v in row.items()} for row in rows)


def load_config(study_root, path=None):
    path = Path(path or Path(study_root) / "scaling.yaml").resolve()
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    cores = int(cfg["resources"]["cores"])
    if not 1 <= cores <= 8:
        raise ValueError("local core cap must be between one and eight")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[key] = "1"
    cfg["_path"] = str(path)
    return cfg


def options(cfg, smoke=False):
    ev = cfg["evaluation"]
    return {"n_radial": 12 if smoke else int(ev["n_radial"]),
            "n_azimuthal": 24 if smoke else int(ev["n_azimuthal"]),
            "chunk_size": int(cfg["resources"]["chunk_size"])}


def mapping_from_config(cfg, family):
    settings = cfg["fit"]
    return CoordinateMap(family=family, degree=int(settings["degree"]),
                         angular_order=int(settings["angular_order"]),
                         feature=settings["feature"], log_bound=float(settings["log_bound"]),
                         quadrature_order=int(settings["quadrature_order"]))


def synthetic_cases(cfg):
    """Small, explicitly synthetic direct-Cai fields; these are not DSMC data."""
    radii = np.geomspace(.25, 60, 12)
    angles = np.deg2rad([0, 20, 40, 60])
    points = np.vstack([np.column_stack((radii * np.cos(a), radii * np.sin(a)))
                        for a in angles])
    limit = int(cfg["resources"]["smoke_max_points_per_case"])
    if limit < len(points):
        points = points[np.linspace(0, len(points)-1, limit, dtype=int)]
    truth = CaiPredictor(CoordinateMap("axial", [-.2]), **options(cfg, True))
    return [sf.ScalingCase(name=f"SYNTHETIC_Kn{kn:g}", kn=kn, points=points.copy(),
                           density=truth.density(points, kn), volumes=np.ones(len(points)),
                           temperature=truth.temperature(points, kn), suitable=True,
                           metadata={"synthetic": True, "truth": truth.mapping.to_dict(),
                                     "uncertainty": "unavailable: deterministic synthetic field"})
            for kn in [1.0, 10.0]]


def load_available(inv, cfg):
    cases, skipped = [], []
    total_bytes = 0
    for row in inv["cases"]:
        if not row.get("compatible") or not row.get("fit_ready"):
            skipped.append({"name": row["name"], "reason": row.get("warnings", []),
                            "result_status": row.get("result_status")})
            continue
        try:
            case = si.load_case(row, max_bytes=int(cfg["resources"]["max_input_bytes"]))
            if case is not None:
                size = sum(array.nbytes for array in (case.points, case.density, case.volumes,
                           case.occupied, case.temperature) if array is not None)
                size += sum(value.nbytes for value in case.metadata.values() if isinstance(value, np.ndarray))
                total_bytes += size
                if total_bytes > int(cfg["resources"].get("max_total_input_bytes", 1073741824)):
                    raise MemoryError("Total input memory cap exceeded; reduce registered input size or raise the explicit cap on HPC")
                case.metadata["stratum"] = row["stratum"]
                cases.append(case)
        except (ValueError, OSError) as exc:
            skipped.append({"name": row["name"], "reason": str(exc)})
    return cases, skipped


def _fit_kwargs(cfg, smoke):
    settings = cfg["fit"]
    return dict(max_nfev=int(cfg["resources"]["smoke_max_nfev"] if smoke else settings["max_nfev"]),
                regularization=float(settings["regularization"]),
                floor=float(settings["floor"]),
                angle_edges=np.deg2rad(settings["angle_edges_deg"]),
                radial_bins=min(12, int(settings["radial_bins"])) if smoke else int(settings["radial_bins"]),
                coefficient_bound=float(settings["coefficient_bound"]),
                evaluate_temperature=bool(cfg["evaluation"]["temperature_transfer"]),
                evaluator_options=options(cfg, smoke))


def run_assessment(cases, cfg, output, *, smoke=False, inv=None, skipped=()):
    """Separate fit tracks and suitability subsets, stratified by source physics."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    records = []
    strata = sorted({str(c.metadata.get("stratum", "synthetic")) for c in cases})
    for index, stratum in enumerate(strata):
        members = [c for c in cases if str(c.metadata.get("stratum", "synthetic")) == stratum]
        for subset in cfg["fit"]["subsets"]:
            selected = [c for c in members if subset != "suitable" or c.suitable]
            for track in cfg["fit"]["tracks"]:
                for family in cfg["fit"]["families"]:
                    key = f"stratum{index}/{subset}/{track}/{family}"
                    destination = output / key
                    mapping = mapping_from_config(cfg, family)
                    result = sf.fit(selected, mapping, track=track,
                                    holdout_kn=tuple(cfg["fit"]["holdout_kn"]),
                                    **_fit_kwargs(cfg, smoke))
                    document = result.to_dict()
                    document.update(key=key, subset=subset, stratum=stratum, family=family,
                                    synthetic=smoke, label=LABELS[family])
                    if int(cfg["fit"].get("block_resamples", 0)):
                        block_options = _fit_kwargs(cfg, smoke)
                        block_options.pop("evaluate_temperature", None)
                        document["block_refits"] = sf.block_refits(
                            selected, result.mapping, track=track,
                            n_resamples=int(cfg["fit"]["block_resamples"]),
                            holdout_kn=tuple(cfg["fit"]["holdout_kn"]), **block_options)
                    document["parameter_uncertainty"] = (
                        "descriptive correlation-qualified block ensemble; not a confidence interval"
                        if document.get("block_refits", {}).get("available") else
                        "statistical parameter uncertainty unavailable")
                    if selected and result.identifiable:
                        predictor = CaiPredictor(result.mapping, speed_ratio=selected[0].speed_ratio,
                                                 **options(cfg, smoke))
                        destination.mkdir(parents=True, exist_ok=True)
                        predictor.save(destination / "predictor.json")
                        train, held, _ = sf.select_cases(selected, holdout_kn=tuple(cfg["fit"]["holdout_kn"]))
                        eval_options = {k: v for k, v in _fit_kwargs(cfg, smoke).items()
                                        if k in ("floor", "angle_edges", "radial_bins", "evaluator_options")}
                        eval_options["temperature"] = bool(cfg["evaluation"]["temperature_transfer"])
                        document["evaluation"] = sf.evaluate(train, result.mapping, **eval_options)
                        document["evaluation_role"] = "training descriptive agreement"
                        document["holdout_evaluation"] = (sf.evaluate(held, result.mapping, **eval_options)
                                                             if held else None)
                        document["diagnostic_sensitivity"] = sf.sensitivity_report(
                            train, result.mapping, evaluator_options=options(cfg, smoke),
                            floor=float(cfg["fit"]["floor"]),
                            angle_widths_deg=(15., 30.),
                            radial_bin_counts=(6, 12) if smoke else (40, 80))
                        document["predictor"] = str((destination / "predictor.json").relative_to(output))
                    write_json(destination / "fit.json", document)
                    records.append(document)
    if not strata:
        records.append({"status": "no eligible whole-field inputs", "identifiable": False,
                        "subsets": {"exploratory": "unidentifiable", "suitable": "unidentifiable"}})
    report = {"schema_version": 1, "synthetic": smoke, "config": cfg,
              "inventory": inv, "skipped": list(skipped), "fits": records,
              "timing_seconds": {"assessment": time.perf_counter()-started},
              "interpretation": "Descriptive empirical composition; no conservation or unseen-Kn claim.",
              "statistical_uncertainty": "unavailable unless independent time-block fields are registered",
              "production_dsmc_executed": False}
    write_json(output / "assessment.json", report)
    write_csv(output / "coefficients.csv", [
        {"key": r.get("key"), "family": r.get("family"), "track": r.get("track"),
         "subset": r.get("subset"), "identifiable": r.get("identifiable"),
         "success": r.get("success"),
         "coefficients": r.get("mapping", {}).get("coefficients"),
         "objective": r.get("diagnostics", {}).get("objective"),
         "L_F": r.get("evaluation", {}).get("L_F"),
         "L_C": r.get("evaluation", {}).get("L_C"),
         "mapping": r.get("mapping"), "diagnostics": r.get("diagnostics")}
        for r in records])
    write_summary(output, report)
    return report


def write_summary(output, report):
    inv = report.get("inventory") or {}
    text = ["# Empirical Cai coordinate scaling", "",
            "SYNTHETIC NUMERICAL CHECK ONLY; no DSMC evidence." if report["synthetic"] else
            "Study-scoped descriptive assessment. No DSMC or scheduler jobs were launched.", "",
            "Density is Cai evaluated at mapped coordinates, without an amplitude or determinant factor.",
            "Temperature is a transfer test using density coefficients. Velocity is unchanged.", "",
            f"Eligible fit records: {sum(bool(r.get('identifiable')) for r in report['fits'])}.",
            f"Missing whole-field results at Kn: {inv.get('missing_results', [])}.", "",
            "| Case | Kn | Source | Fit input | Suitable | Warnings |",
            "|---|---:|---|---|---|---|"]
    for row in inv.get("cases", []):
        text.append(f"| {row['name']} | {row.get('Kn')} | {row.get('source_type')} | "
                    f"{row.get('result_status')} | {row.get('suitable')} | "
                    + "; ".join(row.get("warnings", [])).replace("|", "/") + " |")
    text += ["", "Separate F and C fits, their cross-diagnostics, regional scores, coverage and bound flags",
             "are stored in each fit.json. Exploratory and suitable fits are reported separately.", "",
             "A missing or inadequate suitable subset cannot establish model effectiveness.",
             "The h=1/(1+Kn) feature saturates at low Kn. Better descriptive fit alone does not",
             "establish physical correctness, predictive validity or a self-consistent kinetic solution.", ""]
    (Path(output) / "README.md").write_text("\n".join(text), encoding="utf-8")


def main(argv=None, *, study_root=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", nargs="?", default="inventory",
                        choices=["inventory", "prepare", "fit", "report", "evaluate", "smoke"])
    parser.add_argument("--config")
    parser.add_argument("--output")
    parser.add_argument("--full", action="store_true", help="explicitly request full available-field assessment")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--predictor", help="serialized predictor for evaluate")
    parser.add_argument("--points", help="CSV with X,P columns for evaluate")
    parser.add_argument("--kn", type=float)
    parser.add_argument("--sensitivities", action="store_true", help="refit configured numerical sensitivity cases")
    args = parser.parse_args(argv)
    root = Path(study_root or Path.cwd()).resolve()
    try:
        cfg = load_config(root, args.config)
        output = Path(args.output or root / cfg["output"]).resolve()
        if args.mode == "smoke" and not args.output:
            output = output / "smoke"
        registration = root / cfg["registration"] if cfg.get("registration") else None
        if args.mode == "evaluate":
            if not args.predictor or not args.points or args.kn is None:
                parser.error("evaluate requires --predictor, --points and --kn")
            predictor = CaiPredictor.load(args.predictor)
            points = np.genfromtxt(args.points, delimiter=",", names=True)
            points = np.column_stack((np.atleast_1d(points["X"]), np.atleast_1d(points["P"])))
            begin = time.perf_counter()
            density = predictor.density(points, args.kn)
            temperature = predictor.temperature(points, args.kn)
            write_csv(output / "predictions.csv", [dict(X=p[0], P=p[1], Kn=args.kn,
                      n_over_ne=n, T_over_Te=t) for p, n, t in zip(points, density, temperature)])
            write_json(output / "evaluation-timing.json", {"seconds": time.perf_counter()-begin})
        elif args.mode == "prepare":
            paths = si.prepare_missing_specs(root, output / "specifications", registration)
            print(f"Prepared {len(paths)} configuration-only YAML files; no meshes or DSMC runs.")
        elif args.mode == "report":
            report = json.loads((output / "assessment.json").read_text(encoding="utf-8"))
            write_summary(output, report)
            if not args.no_plots:
                from .scaling_report import render_report
                render_report(output, report)
        else:
            begin = time.perf_counter()
            inv = si.inventory(root, registration)
            write_json(output / "inventory.json", inv)
            write_csv(output / "inventory.csv", inv["cases"])
            if args.mode != "inventory":
                if args.mode == "fit" and not args.full:
                    parser.error("fit requires --full; use smoke for a bounded synthetic numerical check")
                cases, skipped = (synthetic_cases(cfg), []) if args.mode == "smoke" else load_available(inv, cfg)
                preprocessing = time.perf_counter()-begin
                report = run_assessment(cases, cfg, output, smoke=args.mode == "smoke", inv=inv, skipped=skipped)
                report["timing_seconds"]["inventory_and_loading"] = preprocessing
                write_json(output / "assessment.json", report)
                if not args.no_plots:
                    from .scaling_report import render_report
                    render_report(output, report, cases=cases)
                if args.sensitivities:
                    from .scaling_report import run_sensitivities
                    run_sensitivities(cases, cfg, output, smoke=args.mode == "smoke")
        print(f"Scaling {args.mode}: {output}")
        return 0
    except (ValueError, OSError, KeyError, MemoryError) as exc:
        parser.exit(2, f"Allscale: {exc}\n")
