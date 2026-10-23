#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N knx-analyse
#
# The study-level table and figure, once every case's runPost.sh has finished:
# ./AllpostCases --analyse-only, i.e. ./analyse.py. Submit from here:
#
#     cd cases/cai2012-kn-extension && qsub runAnalyse.sh

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the family directory:\n    cd cases/cai2012-kn-extension && qsub runAnalyse.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f study.yaml ] || [ ! -f AllpostCases ]; then
    printf 'FAILED: %s is not the family directory (no study.yaml or AllpostCases). Submit from it:\n    cd cases/cai2012-kn-extension && qsub runAnalyse.sh\n' "$PWD" >&2
    exit 1
fi

# The family's hpc_env.sh: ./ from here.
HPC_ENV=""
for dir in . .. ../..; do
    if [ -f "$dir/hpc_env.sh" ]; then
        HPC_ENV=$dir/hpc_env.sh
        break
    fi
done
if [ -z "$HPC_ENV" ]; then
    printf 'FAILED: no hpc_env.sh in %s or the two directories above it\n' "$PWD" >&2
    exit 1
fi
source "$HPC_ENV"

printf '>>> Starting: the study-level analysis\n'
bash ./AllpostCases --analyse-only
printf '<<< Completed: the study-level analysis\n'
