#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-mesh
#
# Re-mesh this one case. Submit from the case:
#
#     cd Cases/<case> && qsub runMesh.sh
#
# All seven at once, with generation and validation, is ../../runMeshAll.sh.
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f case.yaml ] && [ -f unigas.yaml ] \
    || { printf 'FAILED: %s is not a cai2012-unigas case\n' "$PWD" >&2; exit 1; }
source ../../hpc_env.sh

printf '>>> Starting: run Allmesh for %s\n' "${PWD##*/}"
bash ./Allmesh
printf '<<< Completed: run Allmesh for %s\n' "${PWD##*/}"
