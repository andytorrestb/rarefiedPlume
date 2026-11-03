#!/bin/bash
#PBS -S /bin/bash
#PBS -q normalq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=64
#PBS -N knx-scaling
# Submit manually from this study directory. Analysis only: no DSMC or resubmission.
set -eo pipefail
if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'Submit from cases/cai2012-kn-extension with msub runScaling.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
test -f scaling.yaml
source ./hpc_env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python Allscale.py inventory
# fit --full already renders the report; a separate report step would repeat it.
python Allscale.py fit --full
