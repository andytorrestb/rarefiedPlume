"""Generate the 0/ inflow boundary fields for this case.

Reads case.yaml, evaluates the plume source-flow model over the inflow patch, and
writes 0/boundaryU, 0/boundaryT and 0/boundaryNumberDensity_<species>.

Run AFTER dsmcInitialise, which creates 0/; this overwrites files in it.
See ./Allrun for the full ordering.

    python runInflow.py [case_dir]
"""

import sys
from pathlib import Path

import numpy as np

from plumetools.config import load_case_config
from plumetools.foamio.fields import write_inflow_fields, write_measurement_fields
from plumetools.inflow import area_weighted_number_density, compute_inflow
from plumetools.mesh import read_boundary


def main(argv):
    case_dir = Path(argv[1]) if len(argv) > 1 else Path(__file__).parent
    cfg = load_case_config(case_dir)
    inflow = compute_inflow(case_dir, cfg)
    write_inflow_fields(case_dir / "0", inflow, cfg)

    if cfg.output.dialect == "standard":
        # dsmcFoam constructs its cloud from nine measurement fields and aborts
        # with `cannot find file "0/q"` if any is missing. dsmcInitialise does not
        # create them; OpenFOAM's tutorials ship them in 0.orig instead.
        # Generated here so their boundary entries name this mesh's patches.
        n = len(write_measurement_fields(case_dir / "0", read_boundary(case_dir)))
        print(f"wrote {n} zeroed measurement fields (dsmcFoam requires them)")

    print(f"wrote {len(inflow.labels)} inflow faces to {case_dir / '0'}")
    print(f"  |U|  = {np.linalg.norm(inflow.U, axis=1).mean():.3f} m/s")
    print(f"  T    = {inflow.T[0]:g} K")
    print(f"  rhoN = {inflow.rhoN.min():.4e} .. {inflow.rhoN.max():.4e} 1/m^3")

    if cfg.output.dialect == "standard":
        # Standard dsmcFoam has no per-face number density: FreeStream takes one
        # scalar per species from constant/dsmcProperties. Print the value that
        # preserves total particle flux, and say plainly what it costs -- this is
        # the one place where the solver cannot represent the model.
        mean = area_weighted_number_density(inflow)
        spread = inflow.rhoN.max() / max(inflow.rhoN.min(), 1e-300)
        print()
        print("  Standard dsmcFoam cannot read a per-face number density.")
        print("  0/boundaryU and 0/boundaryT ARE per-face, so the plume's angular")
        print("  direction survives; its angular density profile does not.")
        print()
        print("  Put this in constant/dsmcProperties -> FreeStreamCoeffs:")
        print(f"      numberDensities {{ {cfg.gas.species_name}  {mean:.6e}; }}")
        print(f"  (area-weighted mean; collapses a {spread:.0f}x spread to a constant)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
