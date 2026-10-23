#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N knx-meshall
#
# Generate, validate and mesh all four cases, in one job. Submit from here:
#
#     cd cases/cai2012-kn-extension && qsub runMeshAll.sh
#
#   1. python generate_cases.py               Cases/ and manifest.yaml
#   2. python validate_cases.py --inputs-only the tree against the study
#   3. bash ./AllmeshCases                    mesh all four
#   4. python validate_cases.py               again, with the dictionaries
#
# Any failure fails the job. Then one job per case: cd Cases/<case> && qsub
# runCase.sh (README.md). To re-mesh a single case: Cases/<case>/runMesh.sh.

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the family directory:\n    cd cases/cai2012-kn-extension && qsub runMeshAll.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f study.yaml ] || [ ! -f generate_cases.py ]; then
    printf 'FAILED: %s is not the family directory (no study.yaml or generate_cases.py). Submit from it:\n    cd cases/cai2012-kn-extension && qsub runMeshAll.sh\n' "$PWD" >&2
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

printf '>>> Starting: generate the cases\n'
python generate_cases.py
printf '<<< Completed: generate the cases\n'

printf '>>> Starting: validate the inputs\n'
python validate_cases.py --inputs-only
printf '<<< Completed: validate the inputs\n'

printf '>>> Starting: mesh all four cases\n'
bash ./AllmeshCases
printf '<<< Completed: mesh all four cases\n'

printf '>>> Starting: validate the tree, dictionaries included\n'
python validate_cases.py
printf '<<< Completed: validate the tree, dictionaries included\n'
