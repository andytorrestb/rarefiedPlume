#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N knx-post
#
# Post-process this case once its run has finished: ../../AllpostCases' per-case
# chain (Allpost, Alldiagnose.py, Allcompare.py, Allaudit.py, Allconverge.py),
# without the imagery (runViz.sh) and without the study-level analysis, which
# waits for every case (../../runAnalyse.sh). Submit from the case:
#
#     cd Cases/<case> && qsub runPost.sh
#     cd Cases/Kn1_np1x && qsub -N knx-post-Kn1 runPost.sh    # a per-case job name

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory:\n    cd Cases/<case> && qsub runPost.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f case.yaml ] || [ ! -f ../../AllpostCases ]; then
    printf 'FAILED: %s is not a generated case (no case.yaml, or no ../../AllpostCases). Submit from a case directory:\n    cd Cases/<case> && qsub runPost.sh\n' "$PWD" >&2
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

printf '>>> Starting: post-process %s\n' "$CASE"
bash ../../AllpostCases --case "$CASE" --no-viz --no-analyse > output_post.log 2> error_post.log
printf '<<< Completed: post-process %s\n' "$CASE"
