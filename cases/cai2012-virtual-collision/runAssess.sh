#!/bin/bash
#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=24gb
#PBS -N virtual-collision-assess
# Submit manually from this directory: qsub runAssess.sh
# This is postprocessing only; no DSMC jobs or case generation are invoked.
set -eo pipefail
if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'PBS_O_WORKDIR is unset; submit from cases/cai2012-virtual-collision\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
test -f assess.py
# Follow the existing environment convention, allowing a site-specific override.
# The referenced existing family environment contains site paths; this template
# introduces none. Override VC_HPC_ENV for another cluster/environment.
VC_HPC_ENV=${VC_HPC_ENV:-../cai2012-kn-extension/hpc_env.sh}
test -f "$VC_HPC_ENV"
source "$VC_HPC_ENV"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python assess.py full --output results-hpc
