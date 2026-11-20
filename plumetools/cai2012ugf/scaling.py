"""The scaling inventory and fits of the cai2012-unigas family (spec §10).

The frozen kn-extension inventory
(:func:`plumetools.cai2012ext.scaling_inventory.inventory`) builds every row:
configuration, exit state, field source, fit readiness. It cannot judge these
cases -- its suitability requires the dsmcFoam pre-run feasibility, which USP
cells fail by design -- and it keys strata on physics alone, so the two solvers
would pool. This module relabels each row:

* ``solver`` is ``uniGasFoam hybrid USP-DSMC``;
* ``stratum`` gains ``+uniGasFoam``, so a fit never mixes solvers;
* ``suitable`` is the post-run verdict's ``production_ready`` (spec §7), and
  the row must still be compatible and fit-ready; the reasons are the
  verdict's unresolved items plus the known source limitations.

Fitting and reporting are the frozen workflow's functions, unchanged. Outputs
go under this family's ``results/scaling/``.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from plumetools.cai2012ext import scaling_inventory as si
from plumetools.cai2012ext import scaling_workflow as sw
from plumetools.cai2012ugf.verdict import read_verdict

SOLVER = "uniGasFoam hybrid USP-DSMC"
STRATUM_TAG = "+uniGasFoam"
CRITERIA = "cai2012ugf post-run verdict (spec §7)"


def relabel(row: dict) -> dict:
    """One frozen inventory row, judged by this family's verdict."""
    row = dict(row)
    row["solver"] = SOLVER
    if row.get("stratum"):
        row["stratum"] = f"{row['stratum']}{STRATUM_TAG}"
    old_reasons = set((row.get("suitability") or {}).get("reasons") or [])
    verdict = read_verdict(Path(row["path"]))
    reasons = list(verdict.get("unresolved") or [])
    for limitation in row.get("known_source_limitations") or []:
        if limitation not in reasons:
            reasons.append(limitation)
    suitable = (verdict.get("production_ready") is True and not reasons
                and bool(row.get("compatible")) and bool(row.get("fit_ready")))
    row["suitable"] = suitable
    row["suitability"] = {"suitable": suitable,
                          "status": verdict.get("status", "post-run-verdict-pending"),
                          "reasons": reasons, "criteria": CRITERIA}
    row["warnings"] = [w for w in row.get("warnings") or [] if w not in old_reasons] + reasons
    return row


def inventory(study_root: Path) -> dict:
    """The frozen inventory of this family, relabelled row by row."""
    inv = si.inventory(Path(study_root), None)
    inv["cases"] = [relabel(row) for row in inv["cases"]]
    strata = {}
    for row in inv["cases"]:
        if row.get("compatible"):
            strata.setdefault(row["stratum"], []).append(row["name"])
    inv["compatibility_strata"] = strata
    inv["solver"] = SOLVER
    inv["suitability_criteria"] = CRITERIA
    return inv


def main(argv=None, *, study_root=None) -> int:
    """``inventory``, ``fit --full`` and ``report``, as the frozen workflow's."""
    parser = argparse.ArgumentParser(description="Scaling of the cai2012-unigas family.")
    parser.add_argument("mode", nargs="?", default="inventory",
                        choices=["inventory", "fit", "report"])
    parser.add_argument("--config")
    parser.add_argument("--output")
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    root = Path(study_root or Path.cwd()).resolve()
    try:
        cfg = sw.load_config(root, args.config)
        output = Path(args.output or root / cfg["output"]).resolve()
        if args.mode == "report":
            import json
            report = json.loads((output / "assessment.json").read_text(encoding="utf-8"))
            sw.write_summary(output, report)
            if not args.no_plots:
                from plumetools.cai2012ext.scaling_report import render_report
                render_report(output, report)
            print(f"Scaling report: {output}")
            return 0
        begin = time.perf_counter()
        inv = inventory(root)
        sw.write_json(output / "inventory.json", inv)
        sw.write_csv(output / "inventory.csv", inv["cases"])
        if args.mode == "fit":
            if not args.full:
                parser.error("fit requires --full")
            print(f"Cai evaluation workers: {sw.configure_workers(cfg['resources']['cores'])}",
                  flush=True)
            cases, skipped = sw.load_available(inv, cfg)
            sw.require_loaded(cases, skipped)
            report = sw.run_assessment(cases, cfg, output, inv=inv, skipped=skipped)
            report["timing_seconds"]["inventory_and_loading"] = time.perf_counter() - begin
            report["solver"] = SOLVER
            sw.write_json(output / "assessment.json", report)
            if not args.no_plots:
                from plumetools.cai2012ext.scaling_report import render_report
                render_report(output, report, cases=cases)
        print(f"Scaling {args.mode}: {output}")
        return 0
    except (ValueError, OSError, KeyError, MemoryError, si.InventoryError) as exc:
        parser.exit(2, f"Allscale: {exc}\n")
