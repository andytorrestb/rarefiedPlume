"""Generate the 0/ inflow boundary fields for this case.

Reads case.yaml, evaluates the plume source-flow model over the inflow patch, and
writes 0/boundaryU, 0/boundaryT and 0/boundaryNumberDensity_<species>.

Run AFTER dsmcInitialise+, which creates 0/; this overwrites three files in it.
See ./Allrun for the full ordering.

    python runInflow.py [case_dir]
"""

import sys
from pathlib import Path

from plumetools.config import load_case_config
from plumetools.foamio.fields import write_inflow_fields
from plumetools.inflow import compute_inflow


def main(argv):
    case_dir = Path(argv[1]) if len(argv) > 1 else Path(__file__).parent
    cfg = load_case_config(case_dir)
    inflow = compute_inflow(case_dir, cfg)
    write_inflow_fields(case_dir / "0", inflow, cfg)
    print(f"wrote {len(inflow.labels)} inflow faces to {case_dir / '0'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
