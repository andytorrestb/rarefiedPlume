#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-run
#
# Run this case on uniGasFoam: one PBS job per case, one node. Submit from the case:
#
#     cd Cases/<case> && qsub runCase.sh
#     cd Cases/Kn1_np1x && qsub -N ugf-run-Kn1 runCase.sh    # a per-case job name
#
# The rank count is dsmc.n_subdomains in case.yaml (192), written into
# system/decomposeParDict by ./Allmesh.
#
# In order:
#   1. ./Allrun >> output.log 2> error.log, timed. Allrun checks the built
#      uniGasFoam commit, runs the startup and sampling legs, refuses an unsafe
#      resume (lost particle clouds), reconstructs and converts the fields.
#   2. record the wall clock in ../../results/cost-model.yaml, unless the case
#      had already finished when this job started.
#
# A job that runs out of walltime is resubmitted BY HAND with the same qsub.
# Every case is post-run-verdict-pending until runPost.sh has judged it.
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f case.yaml ] && [ -f unigas.yaml ] \
    || { printf 'FAILED: %s is not a cai2012-unigas case\n' "$PWD" >&2; exit 1; }
CASE=${PWD##*/}
source ../../hpc_env.sh

ALREADY_FINISHED=0
if python -m plumetools.cai2012ugf.hpc finished .; then
    ALREADY_FINISHED=1
    printf '    %s had already finished: no wall clock will be recorded\n' "$CASE"
fi

printf '>>> Starting: run Allrun for %s\n' "$CASE"
START=$(date +%s)
bash ./Allrun >> output.log 2> error.log
ELAPSED=$(( $(date +%s) - START ))
printf '<<< Completed: run Allrun (%d s)\n' "$ELAPSED"

if [ "$ALREADY_FINISHED" -eq 0 ]; then
    printf '>>> Starting: record the wall clock (%d s)\n' "$ELAPSED"
    python -m plumetools.cai2012ugf.hpc record . "$ELAPSED"
    printf '<<< Completed: record the wall clock\n'
fi
