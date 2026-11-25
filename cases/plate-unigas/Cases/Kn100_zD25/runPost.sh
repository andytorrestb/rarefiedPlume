#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N plate-post
#
# Post-process this plate-unigas case once its run has finished: ./Allpost
# (reconstruct, cell centres and volumes, plate loads, the Cai comparison,
# slices) into results/. The whole node is requested for its memory.
#
#     cd Cases/<case> && qsub runPost.sh

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR
if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f case.yaml ] && [ -f unigas.yaml ] \
    || { printf 'FAILED: %s is not a plate-unigas case\n' "$PWD" >&2; exit 1; }
source ../../hpc_env.sh
python -m plumetools.plateugf finished . \
    || printf '    WARNING: %s has not reached its end time; post-processing the newest write\n' "${PWD##*/}"
printf '>>> Starting: run Allpost for %s\n' "${PWD##*/}"
bash ./Allpost
printf '<<< Completed: run Allpost for %s\n' "${PWD##*/}"
