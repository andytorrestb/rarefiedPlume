"""Overlap-validation evidence: what a production coupling verdict may rest on.

A production interface case is a candidate only when a truncated validation
run at the same coupling has been compared with the full-domain benchmark and
passed. :mod:`plumetools.cai2012ext.overlap` writes that comparison to
``results/overlap.yaml``, with the block :func:`evidence_provenance` builds; a
production case names it in ``interface.validation_evidence``, and
:func:`coupling_validation` decides whether it counts.

What has to match
-----------------
The evidence validates one *coupling*: one model, at one Bird threshold and
definition, injecting one core with one closure, from one source problem.
Every item of :data:`REQUIRED` must be recorded in the evidence and equal the
production case's (``purpose`` must be ``interface-validation``: a pass of any
other kind of run is not overlap evidence). The core angle is included because
the overlap region is confined to the injected core and the excluded tail
shapes what reaches it, so a pass at one angle says nothing about another.

What may differ
---------------
The Knudsen number, and with it ``R``, the mesh, the time step and the particle
weight. The validation runs at Cai's Kn = 0.01, where a full-domain benchmark
exists, precisely so that its verdict can be carried to the low-Kn cases. These
are recorded in the evidence for the reader, never compared.

Before 2026-10-27 a ``verdict: pass`` with a matching ``bird_threshold`` was
enough, so a pass under another Bird definition, core angle, closure, tail
treatment or gas -- or a file that recorded nothing else -- validated the case.
Evidence written before then lacks the provenance and is reported as
``insufficient-provenance``; re-running ``Alloverlap.py`` on the validation case
rewrites it from the same fields.
"""

from __future__ import annotations

import math

import yaml

from plumetools.cai2012ext import interface, provenance

#: Provenance items the evidence must record and the production case must match.
REQUIRED = ("model", "purpose", "bird_definition", "bird_threshold",
            "core_half_angle_deg", "tails", "angular_exponent_offset",
            "source_problem")

#: The only kind of run whose verdict is overlap evidence.
EVIDENCE_PURPOSE = "interface-validation"


def source_problem(cfg) -> dict:
    """The exit and gas the closure is built from: Cai's ``D``, ``S0``, ``T0``.

    Not the Knudsen number -- see the module docstring.
    """
    return {"D_m": float(cfg.nozzle.diameter_m),
            "S0": float(cfg.exit.speed_ratio),
            "T0_K": float(cfg.exit.T0_K),
            "gas": provenance.gas_parameters(cfg)}


def evidence_provenance(derived) -> dict:
    """The ``interface`` block of ``overlap.yaml``: the coupling that was tested.

    The first five keys are the ones evidence has always carried; the rest are
    what :func:`coupling_validation` needs to decide whether the evidence fits a
    production case, and the Knudsen number for the reader.
    """
    settings, selection = derived.settings, derived.selection
    d = derived.source.diameter_m
    return {
        "bird_threshold": settings.bird_threshold,
        "bird_definition": settings.bird_definition,
        "core_half_angle_deg": settings.core_half_angle_deg,
        "R_interface_over_D": None if selection.radius_m is None
        else selection.radius_m / d,
        "core_mass_fraction": selection.core_fraction,
        "model": interface.MODEL_NAME,
        "purpose": settings.purpose,
        "tails": settings.tails,
        "angular_exponent_offset": settings.angular_exponent_offset,
        "source_problem": source_problem(derived.cfg),
        "knudsen": float(derived.cfg.exit.knudsen),
    }


def _same(recorded, wanted) -> bool:
    if isinstance(wanted, dict):
        return (isinstance(recorded, dict) and set(recorded) == set(wanted)
                and all(_same(recorded[k], wanted[k]) for k in wanted))
    if isinstance(wanted, bool) or isinstance(recorded, bool):
        return recorded is wanted
    if isinstance(wanted, (int, float)) and isinstance(recorded, (int, float)):
        return math.isclose(float(recorded), float(wanted), rel_tol=1e-12)
    return recorded == wanted


def compatibility(recorded: dict, derived) -> tuple[list[str], list[dict]]:
    """``(missing, mismatches)`` of evidence provenance against a production case."""
    wanted = evidence_provenance(derived)
    wanted["purpose"] = EVIDENCE_PURPOSE
    missing = [key for key in REQUIRED if key not in recorded]
    mismatches = [{"item": key, "evidence": recorded[key], "case": wanted[key]}
                  for key in REQUIRED
                  if key in recorded and not _same(recorded[key], wanted[key])]
    return missing, mismatches


def coupling_validation(derived) -> dict:
    """Whether passed overlap evidence for THIS coupling exists.

    Statuses: ``not-applicable`` (not a production run), ``pending`` (none
    named), ``missing``, ``unreadable``, the evidence's own verdict when it is
    not ``pass``, ``threshold-mismatch``, ``incompatible`` (another coupling),
    ``insufficient-provenance`` (it does not say which coupling), ``passed``.
    """
    settings = derived.settings
    if settings.purpose != "production":
        return {"required": False, "status": "not-applicable",
                "note": f"a {settings.purpose} run does not need it"}
    path = settings.validation_evidence
    if path is None:
        return {"required": True, "status": "pending",
                "note": "no overlap evidence named (interface.validation_evidence)"}
    evidence = (derived.case_dir / path).resolve()
    if not evidence.is_file():
        return {"required": True, "status": "missing", "evidence": str(evidence),
                "note": "the named overlap result is not on this machine"}
    try:
        doc = yaml.safe_load(evidence.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        return {"required": True, "status": "unreadable", "evidence": str(evidence),
                "note": str(exc)}
    verdict = str(doc.get("verdict", "missing"))
    recorded = doc.get("interface") or {}
    threshold = recorded.get("bird_threshold")
    result = {"required": True, "evidence": str(evidence),
              "evidence_threshold": threshold, "verdict": verdict}
    if verdict != "pass":
        return {**result, "status": verdict}
    missing, mismatches = compatibility(recorded, derived)
    if any(m["item"] == "bird_threshold" for m in mismatches):
        status = "threshold-mismatch"
    elif mismatches:
        status = "incompatible"
    elif missing:
        status = "insufficient-provenance"
    else:
        status = "passed"
    result["status"] = status
    if mismatches:
        result["mismatches"] = mismatches
        result["note"] = ("the evidence validated a different coupling: "
                          + ", ".join(m["item"] for m in mismatches))
    elif missing:
        result["missing_provenance"] = missing
        result["note"] = ("the evidence does not record "
                          + ", ".join(missing) + "; re-run Alloverlap.py on the "
                          "validation case to write it")
    return result
