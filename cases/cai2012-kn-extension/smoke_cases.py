#!/usr/bin/env python3
"""Generate and verify input dictionaries without meshing, DSMC or submission.

Existing solver results or meshes are refused. For those trees use
validate_cases.py instead of this preparation smoke check.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import yaml
from plumetools.cai2012 import dictionaries, mesh
from plumetools.cai2012.study import load_study
from plumetools.cai2012ext import generation, hpc, policy, provenance
from plumetools.cai2012ext.verification import verify_extension
from plumetools.cai2012ext.scaling_inventory import inventory


def digest_tree(root, entries):
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for entry in entries for path in sorted((root / entry["path"]).rglob("*"))
            if path.is_file()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--study", default=str(HERE / "study.yaml"))
    parser.add_argument("--report", default=None)
    args = parser.parse_args(argv)
    study_path = Path(args.study).resolve()
    root = study_path.parent
    start = time.perf_counter()
    study = load_study(study_path)
    policy.check_study(study)
    for case in study.expand():
        destination = root / study.case_path(case)
        if generation.has_solver_output(destination) or (destination / "constant" / "polyMesh").exists():
            parser.error(f"{destination} already holds mesh/solver data; use validate_cases.py for an existing HPC run")
    generator = [sys.executable, str(HERE / "generate_cases.py"), "--study", str(study_path)]
    generated = subprocess.run(generator, capture_output=True, text=True, timeout=180)
    if generated.returncode:
        print(generated.stdout + generated.stderr, file=sys.stderr)
        return generated.returncode
    manifest_path = root / study.manifest_name
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    written, guards = [], {}
    for entry in manifest["cases"]:
        directory = root / entry["path"]
        case = generation.derive(directory)
        case.report.raise_if_failed()
        refusal = hpc.guard(directory)
        if refusal:
            raise ValueError(f"Runtime policy refuses prepared case {entry['name']}: {refusal}")
        guards[entry["name"]] = "accepted by declared profile; suitability is separate"
        written.extend(mesh.write_mesh_setup(directory, case.cfg, case.geom, case.plan))
        written.extend(dictionaries.write_case_dictionaries(directory, case.cfg,
                                                            case.exit_state, case.run))
    verification = verify_extension(root, manifest, with_dictionaries=True)
    if not verification.passed:
        raise ValueError("Generated suite did not verify:\n" + "\n".join(verification.failures))
    before = digest_tree(root, manifest["cases"])
    with tempfile.TemporaryDirectory(prefix=".staging-smoke-", dir=root) as scratch:
        states = generation.classify(root, study, Path(scratch), study.expand())
        if any(state.status != generation.STATUS_IDENTICAL for state in states):
            raise ValueError("A repeated generation would change case inputs or runner files")
    repeated = subprocess.run(generator, capture_output=True, text=True, timeout=180)
    if repeated.returncode or digest_tree(root, manifest["cases"]) != before:
        raise ValueError("Repeated generation changed the prepared case trees")
    registration = root / "scaling-registration.yaml"
    scaling = inventory(root, registration if registration.exists() else None)
    active_paths = {str((root / entry["path"]).resolve()) for entry in manifest["cases"]}
    registered_paths = {row["path"] for row in scaling["cases"] if row["compatible"]}
    if not active_paths <= registered_paths:
        raise ValueError("Scaling inventory did not register every prepared direct-exit case")
    document = {"revision": provenance.git_provenance(REPO), "study": str(study_path),
                "profile": policy.suite_name(study), "manifest": str(manifest_path),
                "Kn_e": [entry["Kn"] for entry in manifest["cases"]],
                "case_count": len(manifest["cases"]), "dictionary_files_written": len(written),
                "verification_passes": len(verification.passes), "idempotent": True,
                "runtime_guards": guards, "scaling_registration_complete": True,
                "missing_whole_field_results": scaling["missing_results"],
                "cases": [{key: entry.get(key) for key in (
                    "name", "path", "Kn", "n0_per_m3", "production_status", "known_source_limitations")}
                    for entry in manifest["cases"]],
                "mesh_generated": False, "DSMC_executed": False, "scheduler_submitted": False,
                "seconds": time.perf_counter()-start}
    report_path = Path(args.report).resolve() if args.report else root / "results" / "full-suite-preparation" / "smoke.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(document, indent=2)+"\n", encoding="utf-8")
    print(f"PASS: {document['case_count']} cases; {len(written)} dictionary files; "
          f"{len(verification.passes)} verification checks; idempotence and scaling registration.")
    print("No meshes, DSMC runs, or scheduler jobs were launched.")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
