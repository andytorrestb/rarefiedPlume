#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N knx-run
#
# Run this case: one PBS job per case, one node. Submit from the case:
#
#     cd Cases/<case> && qsub runCase.sh
#     cd Cases/Kn1_np1x && qsub -N knx-run-Kn1 runCase.sh    # a per-case job name
#
# The rank count is dsmc.n_subdomains in case.yaml (192). ./Allmesh wrote it
# into system/decomposeParDict and ./Allrun reads it from there; nothing here
# overrides it.
#
# In order:
#   1. refuse a case at Cai's Kn = 0.01, 0.1 or 100           (hpc guard)
#   2. re-assess the case; if it is unsuitable for production, write a
#      NON-PRODUCTION marker into it, and the same text at the top of
#      output.log. It still runs.                               (hpc label)
#   3. ./Allrun >> output.log 2> error.log, timed
#   4. record the wall clock in ../../results/cost-model.yaml, unless the case
#      had already reached endTime when this job started (hpc finished, record)
#
# A job that runs out of walltime is resubmitted BY HAND with the same qsub:
# ./Allrun continues from the newest time directory. Nothing resubmits itself,
# and every time directory and parcel cloud is kept.

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory:\n    cd Cases/<case> && qsub runCase.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f case.yaml ]; then
    printf 'FAILED: %s has no case.yaml. Submit from a case directory:\n    cd Cases/<case> && qsub runCase.sh\n' "$PWD" >&2
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

printf '>>> Starting: assess %s for production\n' "$CASE"
python -m plumetools.cai2012ext.hpc label . | tee output.log
printf '<<< Completed: assess %s for production\n' "$CASE"

ALREADY_FINISHED=0
if python -m plumetools.cai2012ext.hpc finished .; then
    ALREADY_FINISHED=1
    printf '    %s had already reached endTime: no wall clock will be recorded\n' "$CASE"
fi

printf '>>> Starting: run Allrun\n'
START=$(date +%s)
bash ./Allrun >> output.log 2> error.log
ELAPSED=$(( $(date +%s) - START ))
printf '<<< Completed: run Allrun (%d s)\n' "$ELAPSED"

if [ "$ALREADY_FINISHED" -eq 0 ]; then
    printf '>>> Starting: record the wall clock (%d s)\n' "$ELAPSED"
    python -m plumetools.cai2012ext.hpc record . "$ELAPSED"
    printf '<<< Completed: record the wall clock\n'
fi
