#!/usr/bin/env python3
"""Generate this case's 0/ fields from case.yaml, and verify the inflow flux.

    python runInflow.py [case_dir]

Writes, into 0/:

    boundaryU                    per-face source-flow velocity
    boundaryT                    per-face temperature; wall temperature on the bodies
    boundaryNumberDensity_N2     per-face number density, read by plumeFieldInflow
    nine zeroed measurement fields that DSMCCloud constructs itself from

Run BEFORE dsmcInitialise, not after. dsmcInitialise constructs a DSMCCloud,
whose constructor reads 0/boundaryT and 0/boundaryU with MUST_READ -- so they
must already exist. (The legacy 3d-inflow case runs them the other way round and
its Allrun documents the hazard; here the ordering removes it.)
"""

from __future__ import annotations

import sys
from pathlib import Path

from plumetools.config import load_case_config
from plumetools.markelov1999 import checks as checks_module
from plumetools.markelov1999 import flux as flux_module
from plumetools.markelov1999 import foamfields, inflow as inflow_module
from plumetools.markelov1999 import mesh as mesh_module
from plumetools.markelov1999 import resolution as resolution_module
from plumetools.markelov1999.geometry import from_config
from plumetools.mesh import read_boundary


def main(argv) -> int:
    case_dir = Path(argv[1]) if len(argv) > 1 else Path(__file__).resolve().parent
    cfg = load_case_config(case_dir)
    geom = from_config(cfg)
    patches = read_boundary(case_dir)

    inflow = inflow_module.compute_inflow(case_dir, cfg)
    print("\n".join(inflow.summary()))
    print()

    written = foamfields.write_all(case_dir, inflow, cfg, patches)
    print(f"wrote {len(written)} field files to {case_dir / '0'}")
    for path in written[:3]:
        print(f"    {path.name}")
    print(f"    ... and {len(written) - 3} zeroed measurement fields")
    print()

    # Flux verification. Three routes: closed-form half-range moments, an
    # independent quadrature, and the same moments recomputed from the values as
    # they were written to disk.
    verification = flux_module.verify(inflow, cfg)
    print("\n".join(verification.report()))
    print()

    cell_sizes = mesh_module.surface_cell_sizes(cfg)
    estimate = resolution_module.estimate(
        cfg, geom, cell_sizes,
        n_equivalent_particles=cfg.dsmc.n_equivalent_particles)
    print("\n".join(estimate.report()))
    print()

    report = checks_module.run_checks(
        cfg, geom, inflow=inflow, patches=patches,
        estimate=estimate, cell_sizes=cell_sizes)
    print("\n".join(report.report()))

    summary = {
        "case": dict(cfg.meta),
        "geometry": geom.to_metadata(),
        "inflow": {
            "n_faces": inflow.n_faces,
            "area_m2": float(inflow.areas.sum()),
            "v_limit_m_per_s": inflow.v_limit,
            "temperature_K": float(inflow.T[0]),
            "number_density_min_per_m3": float(inflow.rhoN.min()),
            "number_density_max_per_m3": float(inflow.rhoN.max()),
            "area_weighted_number_density_per_m3":
                inflow_module.area_weighted_number_density(inflow),
        },
        "flux": {
            "analytical": verification.analytical.as_dict(),
            "quadrature": verification.quadrature.as_dict(),
            "from_fields": verification.from_fields.as_dict(),
            "relative_errors": verification.relative_errors(),
            "quadrature_agreement": verification.quadrature_errors(),
            "drift_only_number_per_s": verification.drift_only_number_per_s,
        },
        "resolution_estimate": {
            "n_equivalent_particles": estimate.n_equivalent_particles,
            "sizing_region": estimate.sizing_region,
            "target_particles_per_cell": estimate.target_particles_per_cell,
            "particles_per_cell": {
                r.region.name: r.particles_per_cell for r in estimate.regions},
            "cell_over_mean_free_path": {
                r.region.name: r.cell_over_mfp for r in estimate.regions},
            "total_particles": estimate.total_particles,
        },
        "checks": report.as_dict(),
    }
    path = checks_module.write_summary(case_dir, summary)
    print()
    print(f"wrote {path}")

    report.raise_if_failed()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv))
    except checks_module.DsmcCheckError as exc:
        print(f"\n{exc}", file=sys.stderr)
        raise SystemExit(1) from None
