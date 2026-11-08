# Run a radial sizing case

Use WSL/Linux, Python with this checkout on its import path, and the installed
OpenFOAM v2512 hyStrath port. The owner's setup is an example, not a runner
dependency:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
export PATH="$HOME/OpenFOAM/andy-v2512/platforms/linux64GccDPInt32Opt/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/OpenFOAM/andy-v2512/platforms/linux64GccDPInt32Opt/lib:${LD_LIBRARY_PATH:-}"
cd /path/to/rarefiedPlume
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
export PYTHON=/path/to/python-with-numpy-scipy-pyyaml
"$PYTHON" cases/cai2012-radial-resolution/generate_cases.py --dry-run
```

Inspect the editable [study settings](../../cases/cai2012-radial-resolution/study.yaml).
The report gives seven initial plans, fixed domains, breakdown outcomes, cells,
weights, predicted populations, steps and budget gaps. At the default settings,
the densest case exceeds the ten-million-step budget. No resource limit relaxes
the mesh criterion. Generation writes dictionaries only:

```bash
"$PYTHON" cases/cai2012-radial-resolution/generate_cases.py
cd cases/cai2012-radial-resolution/Cases/Kn100/r0_ppc20
bash ./Allmesh
bash ./Allrun
bash ./Allpost
```

Each case starts empty and uses its fixed mapped inlet. `Allrun` refuses existing
output unless `--resume` is explicit. `bash ./Allrun --resume` preserves the cloud
and checks the newest complete checkpoint. Incomplete checkpoints produce the
shared recovery instructions and are not deleted. Changed inputs or mesh refuse
restart. Neither generation nor refinement launches a solver.

From the repository root, inspect the selected measured corrections:

```bash
"$PYTHON" cases/cai2012-radial-resolution/refine_case.py \
  cases/cai2012-radial-resolution/Cases/Kn100/r0_ppc20 --dry-run
"$PYTHON" cases/cai2012-radial-resolution/refine_case.py \
  cases/cai2012-radial-resolution/Cases/Kn100/r0_ppc20
"$PYTHON" cases/cai2012-radial-resolution/analyse.py
```

If necessary this creates the next `rN_ppc20` beside the parent. Run its own
Allmesh/Allrun/Allpost and reassess. `--correction mesh`, `parcels`, or `timestep`
isolates a correction; `--delta-t SECONDS` can request a smaller fixed timestep.
An exported portable `--plan DIRECTORY` is also accepted. A passing assessment
creates nothing. Missing measurements request evidence instead of guessing.

To change physical settings, use a separate study/output tree. For explicit
higher-occupancy experiments, copy study.yaml, change `targets.min_mean_parcels`,
and pass `generate_cases.py --study COPY --output NEW_TREE`. Matching generation
is idempotent and preserves results; incompatible existing inputs are refused.

Read the [physical and numerical method](../cai2012-radial-resolution.md) and the
[portable unrelated-case workflow](../mesh-refinement.md) before interpreting
streamtube totals as guidance for another geometry.
