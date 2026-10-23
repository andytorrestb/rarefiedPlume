#!/usr/bin/env python3
"""INTERFACE MODE: write ``0/`` from the closure state on the meshed sphere.

    ./runInflow.py

Called by ``./Allrun`` before a fresh run. Writes, per face of the ``nozzle``
patch in patch order -- the downstream continuum interface, not the exit:

    0/boundaryU                   U(r) r_hat, the closure velocity
    0/boundaryT                   T(r), positive on every face
    0/boundaryNumberDensity_Ar    n(r, theta) in the core, exactly 0 beyond it
    0/{dsmcRhoN,fD,iDof,internalE,linearKE,momentum,q,rhoM,rhoN}

and ``results/interface-inflow.yaml``: the measured patch and the gross and
net mass and energy it injects against the core share of the exit's flow.

Refuses -- writing nothing -- if the mesh was not built for this configuration
(``interface-identity.yaml``), the patch is not the sphere, a written field does
not read back in patch order, or a scientific case's meshed flux misses its
target by more than ``meshed_flux_tolerance``. See
``plumetools/cai2012ext/interface_case.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    from plumetools.cai2012ext import interface_case
except ImportError:
    sys.exit("plumetools is not importable. From the repository root:\n"
             '    pip install -e ".[test]"')


if __name__ == "__main__":
    unknown = sys.argv[1:]
    if unknown:
        print(f"runInflow: unknown option(s) {unknown}", file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(interface_case.main(["inflow", str(HERE)]))
