"""Reports of an interface case: the report file, the manifest row, the
console lines, the NON-PRODUCTION text, and the evidence ``collect`` brings
back from a cluster.

These format what :func:`~plumetools.cai2012ext.interface_case.derive` and
:func:`~plumetools.cai2012ext.interface_assessment.assess` produced; the
comparisons with the 1999 terminal model and the decision table are the only
quantities computed here, and only for the report. Split out of
:mod:`plumetools.cai2012ext.interface_case`.
"""

from __future__ import annotations

import math
from collections import deque
from pathlib import Path

import numpy as np
import yaml

from plumetools.cai2012ext import interface
from plumetools.cai2012ext.interface_assessment import STATUS_CANDIDATE, assess

#: The per-case report, written by ``Allmesh`` and ``runInflow.py``.
REPORT_FILE = Path("results") / "interface-report.yaml"


def report(derived, assessment: dict | None = None, *,
           decision_table: bool = False) -> dict:
    """Everything about the case, for ``results/interface-report.yaml``."""
    assessment = assessment or assess(derived)
    source, selection, settings = derived.source, derived.selection, derived.settings
    d = source.diameter_m
    out = {
        "model": interface.MODEL_NAME,
        "purpose": settings.purpose,
        "non_production": assessment["status"] != STATUS_CANDIDATE,
        "settings": settings.as_dict(),
        "source": source.as_dict(),
        "selection": selection.as_dict(d, source.r_a),
        "diagnostic_radius_over_D": settings.diagnostic_radius_over_D,
        "diagnostic_radius_reproduces_cai": settings.diagnostic_radius_over_D == 10.0,
        "assessment": assessment,
        "terminal_model_comparison": interface.compare_with_terminal_model(
            source, sorted({1.0, 2.0, 5.0, 10.0} | (
                {selection.radius_m / d} if selection.admissible else set()))),
    }
    if decision_table:
        out["angular_shape_comparison"] = interface.angular_shape_comparison(source)
    if derived.meshable:
        out.update({
            "mesh": derived.planned.mesh_plan.as_dict(),
            "run": derived.planned.run.as_dict(),
            "reference_cell": derived.planned.reference,
            "estimated_parcels": derived.planned.parcels,
            "ideal_sphere_flux_balance": derived.planned.ideal_flux,
            "resolution_estimate": derived.planned.resolution,
        })
    if decision_table:
        r_max, _ = interface.domain_radius_limit(derived.cfg, settings)
        r_valid = settings.source_validity_min_radius_over_D or (source.r_a / d)
        out["decision_table"] = {
            "note": "what each (core angle, threshold) WOULD give at "
                    f"source_validity_min_radius_over_D = {r_valid:.4g}; for "
                    "choosing them before any DSMC result, not against one",
            "rows": interface.radius_table(
                source, settings, angles_deg=(10.0, 20.0, 30.0, 45.0, 60.0),
                thresholds=sorted({0.02, 0.03, 0.05, settings.bird_threshold}),
                r_min_over_D=r_valid, r_max_m=r_max),
        }
    return out


def manifest_entry(root: Path, name: str, derived, assessment: dict,
                   *, group: str) -> dict:
    """One ``interface-manifest.yaml`` row."""
    selection, source, settings = derived.selection, derived.source, derived.settings
    d = source.diameter_m
    entry = {
        "name": name,
        "path": derived.case_dir.resolve().relative_to(Path(root).resolve()).as_posix(),
        "group": group,
        "mode": "interface",
        "model": interface.MODEL_NAME,
        "Kn": float(derived.cfg.exit.knudsen),
        "purpose": settings.purpose,
        "run_mode": settings.run_mode,
        "bird_threshold": float(settings.bird_threshold),
        "bird_definition": settings.bird_definition,
        "core_half_angle_deg": settings.core_half_angle_deg,
        "selection_status": selection.status,
        "R_interface_over_D": None if selection.radius_m is None
        else selection.radius_m / d,
        "r_a_over_D": source.r_a / d,
        "max_core_bird_at_r_min": selection.max_core_bird_at_r_min,
        "axis_bird_at_r_min": selection.axis_bird_at_r_min,
        "core_mass_fraction": selection.core_fraction,
        "status": assessment["status"],
        "numerical_status": assessment["numerical"]["status"],
        "coupling_validation": assessment["coupling_validation"]["status"],
        "runnable": assessment["runnable"],
        "unresolved": assessment["unresolved"],
        "exit": {"D_m": d, "S0": float(derived.cfg.exit.speed_ratio),
                 "T0_K": float(derived.cfg.exit.T0_K), "n0_per_m3": source.n_e,
                 "U0_m_per_s": source.U_e,
                 "mdot_target_kg_per_s": source.targets.mdot_kg_per_s,
                 "h_target_J_per_kg": source.targets.h_J_per_kg},
    }
    if derived.meshable:
        entry.update({
            "n_cells": int(derived.planned.mesh_plan.retained_cells),
            "core_cell_size_m": derived.planned.mesh_plan.core_cell_size_m,
            "cell_over_interface_mfp": derived.planned.mesh_plan.cell_over_mfp,
            "deltaT_s": float(derived.planned.run.delta_t_s),
            "n_equivalent_particles": float(derived.planned.run.n_equivalent_particles),
            "average_start_s": float(derived.planned.run.average_start_s),
            "end_time_s": float(derived.planned.run.end_time_s),
            "n_steps": int(derived.planned.run.n_steps),
            "estimated_parcels": float(derived.planned.parcels),
        })
    return entry


def describe(derived, assessment: dict) -> list[str]:
    """Report lines for the console."""
    s, src = derived.selection, derived.source
    d = src.diameter_m
    lines = [
        f"          interface {s.status}: {s.reason}",
        f"          source    r_a = {src.r_a / d:.4f} D, M_a = {src.M_a:.4f}, "
        f"T_t = {src.T_t:.2f} K, V_inf = {src.V_inf:.2f} m/s, core fraction "
        f"{s.core_fraction if s.core_fraction is None else round(s.core_fraction, 4)}",
    ]
    if derived.meshable:
        p, r = derived.planned.mesh_plan, derived.planned.run
        lines += [
            f"          mesh      R = {p.radius_m / d:.4f} D, {p.retained_cells:,} "
            f"cells (lattice {p.lattice.n_cells:,}), core {p.core_cell_size_m:.4e} m, "
            f"cell/lambda {p.cell_over_mfp:.3g}",
            f"          step      deltaT {r.delta_t_s:.4e} s ({r.delta_t_source}), "
            f"Courant {r.courant:.3f}, weight {r.n_equivalent_particles:.4e}, "
            f"~{derived.planned.parcels:.2e} parcels",
        ]
    lines.append(f"          verdict   {assessment['status'].upper()}"
                 f" (numerical: {assessment['numerical']['status']}, coupling: "
                 f"{assessment['coupling_validation']['status']})")
    for reason in assessment["unresolved"]:
        lines.append(f"          UNRESOLVED {reason}")
    return lines


def write_report(case_dir: Path, document: dict) -> Path:
    path = Path(case_dir) / REPORT_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("# GENERATED by plumetools.cai2012ext.interface_case -- do "
                    "not edit.\n" + yaml.safe_dump(_plain(document), sort_keys=False,
                                                   default_flow_style=False),
                    encoding="utf-8", newline="\n")
    return path


def _plain(value):
    """numpy scalars and tuples to plain YAML types."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def label_text(assessment: dict) -> str:
    return ("This interface case was run by its runCase.sh job.\n"
            f"Status: {assessment['status']} (purpose {assessment['purpose']}, "
            f"run_mode {assessment['run_mode']}). Nothing from it is a "
            "validated production result.\n\nUnresolved:\n"
            + "".join(f"  - {line}\n" for line in assessment["unresolved"]))


def collect(case_dir: Path) -> dict:
    """Concise evidence from a run: verdicts, metrics, resources, log tails.

    Reads only small files the case already wrote -- never a field -- so it
    can be run on the cluster and the output copied to the laptop.
    """
    case_dir = Path(case_dir)
    out = {"case": case_dir.resolve().name}
    for key, relative in (("interface_report", REPORT_FILE),
                          ("inflow", Path("results") / "interface-inflow.yaml"),
                          ("overlap", Path("results") / "overlap.yaml"),
                          ("comparison", Path("results") / "comparison.yaml"),
                          ("diagnostics", Path("results") / "diagnostics_summary.yaml"),
                          ("stationarity", Path("results") / "stationarity.yaml")):
        path = case_dir / relative
        if path.is_file():
            try:
                doc = yaml.safe_load(path.read_text(encoding="utf-8"))
            except yaml.YAMLError:
                doc = {"unreadable": str(path)}
            out[key] = _summarise(key, doc)
    marker = case_dir / "NON-PRODUCTION"
    out["non_production_marker"] = marker.read_text(encoding="utf-8") \
        if marker.is_file() else None
    for log in ("log.dsmcFoam", "log.snappyHexMesh", "log.checkMesh",
                "output.log", "error.log"):
        path = case_dir / log
        if path.is_file():
            out.setdefault("logs", {})[log] = log_summary(path)
    return out


def log_summary(path: Path) -> dict:
    """A log's line count, last 15 lines, last 5 memory lines, last ExecutionTime.

    One streamed pass holding only those few lines: a ``log.dsmcFoam`` from a
    long run is gigabytes, and reading it whole to keep twenty lines of it
    needed several times its size in memory. Lines are split exactly as
    ``str.splitlines`` splits the whole text.
    """
    count = 0
    tail, memory, execution = deque(maxlen=15), deque(maxlen=5), deque(maxlen=1)
    with open(path, encoding="utf-8", errors="replace", newline="") as handle:
        for physical in handle:
            for line in physical.splitlines():
                count += 1
                tail.append(line)
                if "MB" in line or "memory" in line.lower():
                    memory.append(line)
                if line.startswith("ExecutionTime"):
                    execution.append(line)
    return {"lines": count, "tail": list(tail), "memory_lines": list(memory),
            "execution_time": list(execution)}


def _summarise(key: str, doc):
    if not isinstance(doc, dict):
        return doc
    if key == "interface_report":
        keep = ("purpose", "non_production", "selection", "assessment",
                "reference_cell", "estimated_parcels", "ideal_sphere_flux_balance")
        summary = {k: doc[k] for k in keep if k in doc}
        if "assessment" in summary:
            summary["assessment"] = {k: summary["assessment"].get(k) for k in
                                     ("status", "runnable", "unresolved",
                                      "numerical", "coupling_validation")}
        if "selection" in summary:
            summary["selection"].pop("scan", None)
        return summary
    return doc
