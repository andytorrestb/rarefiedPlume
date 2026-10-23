#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N knx-mesh
#
# Re-mesh this one case. Submit from the case:
#
#     cd Cases/<case> && qsub runMesh.sh
#     cd Cases/Kn1_np1x && qsub -N knx-mesh-Kn1 runMesh.sh    # a per-case job name
#
# All four at once, with generation and validation, is ../../runMeshAll.sh.
# A case at Cai's Kn = 0.01, 0.1 or 100 is refused (hpc guard), as
# ../../AllmeshCases refuses a manifest that lists one.

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory:\n    cd Cases/<case> && qsub runMesh.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f case.yaml ]; then
    printf 'FAILED: %s has no case.yaml. Submit from a case directory:\n    cd Cases/<case> && qsub runMesh.sh\n' "$PWD" >&2
    exit 1
fi
CASE=${PWD##*/}

# The family's hpc_env.sh: ../../ from Cases/<case>/, ../ from baseCase/.
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

printf '>>> Starting: check the Knudsen number of %s\n' "$CASE"
python -m plumetools.cai2012ext.hpc guard .
printf '<<< Completed: check the Knudsen number of %s\n' "$CASE"

printf '>>> Starting: run Allmesh for %s\n' "$CASE"
bash ./Allmesh
printf '<<< Completed: run Allmesh for %s\n' "$CASE"
