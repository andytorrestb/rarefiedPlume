#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N cai2012-nPartStudy

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set; submit this job from the repository root.\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -x ./Allsetup ]; then
    printf 'FAILED: %s is not the repository root (./Allsetup is missing or not executable).\n' "$PWD" >&2
    exit 1
fi

source ./cases/cai2012-kn-extension/hpc_env.sh

# #./AllrunCases --cores 192 > output.log 2> error.log
# ./AllpostCases > output_post.log 2> error_post.log

./Allsetup