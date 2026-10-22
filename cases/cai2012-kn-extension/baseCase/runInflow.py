#!/usr/bin/env python3
"""Write this case's ``0/`` fields from case.yaml and the meshed nozzle patch.

    ./runInflow.py                write 0/ and verify
    ./runInflow.py --no-checks    write 0/ without running the DSMC checks

Writes:

    0/boundaryU                   U0 along +x on the nozzle
    0/boundaryT                   T0 on the nozzle
    0/boundaryNumberDensity_Ar    n0 on the nozzle, read per face
    0/{dsmcRhoN,fD,iDof,internalE,linearKE,momentum,q,rhoM,rhoN}

and refuses to write a case whose DSMC quality checks fail.

Nothing here is read from a dictionary. ``n0`` comes from the Knudsen number,
``U0`` from the speed ratio, and both are recomputed from ``case.yaml`` every
time -- so the fields the solver injects cannot drift from the configuration
that describes them.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012 import checks, foamfields, inflow, mesh
    from plumetools.cai2012.config import load_case_config
    from plumetools.cai2012.geometry import from_config
    from plumetools.mesh.boundary import read_boundary
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


def main(argv) -> int:
    run_checks = "--no-checks" not in argv[1:]
    unknown = [a for a in argv[1:] if a != "--no-checks"]
    if unknown:
        print(f"runInflow: unknown option(s) {unknown}", file=sys.stderr)
        return 2

    cfg = load_case_config(HERE)
    geom = from_config(cfg)
    exit_state = inflow.from_config(cfg)

    if not (HERE / "constant" / "polyMesh" / "boundary").is_file():
        print("runInflow: no constant/polyMesh. Generate it first:\n"
              "    ./Allmesh", file=sys.stderr)
        return 1
    patches = read_boundary(HERE)

    print()
    print("\n".join(exit_state.describe()))
    print()

    plan = mesh.plan(cfg, geom, exit_state)
    run = inflow.derive_run_settings(
        cfg, exit_state, geom,
        min_cell_size_m=plan.min_cell_size_m,
        exit_cell_volume_m3=plan.exit_cell_volume_m3)
    nozzle = mesh.predicted_nozzle_faces(plan, geom)

    name = cfg.mesh.patch_names["nozzle"]
    meshed_faces = patches[name].n_faces if name in patches else 0
    # The meshed area is what the solver injects through, so the theoretical
    # flow rate is quoted against it rather than against pi R0^2.
    meshed_area = nozzle.area_m2 * (meshed_faces / max(1, nozzle.n_faces))
    rate = exit_state.injection_rate_per_s(meshed_area)
    print(f"  nozzle patch             {name}: {meshed_faces} faces, "
          f"{meshed_area:.6e} m^2")
    print(f"  theoretical injection    {rate:.6e} molecules/s")
    print(f"                           = {rate / run.n_equivalent_particles:.6e} "
          f"parcels/s = "
          f"{rate * run.delta_t_s / run.n_equivalent_particles:.1f} parcels/step")
    print()

    report = None
    if run_checks:
        report = checks.run_checks(cfg, geom, exit_state, plan, run,
                                   patches=patches, nozzle=nozzle)
        print("\n".join(report.report()))
        print()
        try:
            report.raise_if_failed()
        except checks.DsmcCheckError as exc:
            print(f"runInflow: {exc}", file=sys.stderr)
            print("\nNo 0/ fields were written. Fix case.yaml and re-run "
                  "./Allmesh.", file=sys.stderr)
            return 1

    for path in foamfields.write_all(HERE, cfg, exit_state, patches):
        print(f"    wrote {path}")

    checks.write_summary(HERE, {
        "case": dict(cfg.meta),
        "exit": exit_state.as_dict(),
        "mesh": plan.as_dict(),
        "run": run.as_dict(),
        "nozzle": {
            "patch": name,
            "predicted_faces": nozzle.n_faces,
            "meshed_faces": meshed_faces,
            "area_m2": meshed_area,
            "exact_area_m2": nozzle.exact_area_m2,
            "area_error": nozzle.area_error,
        },
        "theoretical_injection_per_s": rate,
        "checks": report.as_dict() if report is not None else None,
    })
    print()
    print("    wrote case-summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
