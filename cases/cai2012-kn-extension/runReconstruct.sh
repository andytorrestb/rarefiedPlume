#PBS -S /bin/bash
#PBS -q normalq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -N knx-reconstruct
#
# Finish the reconstruction of every case whose runCase.sh ran out of walltime
# during reconstructPar: ./AllreconstructCases. Submit from here, and before
# runAnalyse.sh or any runPost.sh of those cases:
#
#     cd cases/cai2012-kn-extension && qsub runReconstruct.sh

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the family directory:\n    cd cases/cai2012-kn-extension && qsub runReconstruct.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f study.yaml ] || [ ! -f AllreconstructCases ]; then
    printf 'FAILED: %s is not the family directory (no study.yaml or AllreconstructCases). Submit from it:\n    cd cases/cai2012-kn-extension && qsub runReconstruct.sh\n' "$PWD" >&2
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

printf '>>> Starting: finish the reconstruction\n'
bash ./AllreconstructCases
printf '<<< Completed: finish the reconstruction\n'
