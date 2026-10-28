"""Overlap-validation evidence: what a production coupling verdict may rest on.

A production interface case is a candidate only when a truncated validation
run at the same coupling has been compared with the full-domain benchmark and
passed (:mod:`plumetools.cai2012ext.overlap` writes that ``overlap.yaml``).
This module owns that contract: which evidence a case names, how it is read,
and when it counts. Split out of :mod:`plumetools.cai2012ext.interface_case`.
"""

from __future__ import annotations

import math

import yaml

def coupling_validation(derived) -> dict:
    """Whether overlap evidence at THIS threshold has passed."""
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
    threshold = (doc.get("interface") or {}).get("bird_threshold")
    same = threshold is not None and math.isclose(
        float(threshold), float(settings.bird_threshold), rel_tol=1e-12)
    status = "passed" if verdict == "pass" and same else (
        "threshold-mismatch" if verdict == "pass" else verdict)
    return {"required": True, "status": status, "evidence": str(evidence),
            "evidence_threshold": threshold, "verdict": verdict}
