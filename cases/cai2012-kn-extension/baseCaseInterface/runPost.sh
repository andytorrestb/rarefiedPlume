#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N knx-ipost
#
# INTERFACE MODE. Post-process this interface case once its run has finished.
# Interface cases are not in ../../manifest.yaml, so this runs the per-case
# chain itself rather than ../../AllpostCases:
#
#   Allpost          centreline n, U, T and the n contour, vs Cai's exit solution
#   Alldiagnose.py   rarefaction and resolution fields (no exit-plane check:
#                    the exit is not in the mesh)
#   Allcompare.py    r/D = 10 (or the case's diagnostic radius) and the plane,
#                    vs Cai's exit solution, with the interface mask recorded
#   Allconverge.py   stationarity across the sampled frames
#   Alloverlap.py    interface-VALIDATION cases only: against the full-domain
#                    Kn = 0.01 benchmark
#   collect          results/interface-evidence.yaml, the small file to bring
#                    back to the laptop
#
# Allaudit.py (the direct-exit occupancy regions) is not run here.
#
#     cd Interface/<case> && qsub -N knx-ipost-<case> runPost.sh

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from the case directory:\n    cd Interface/<case> && qsub runPost.sh\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
if [ ! -f case.yaml ] || [ ! -f interface.yaml ]; then
    printf 'FAILED: %s is not an interface case (no case.yaml or interface.yaml)\n' "$PWD" >&2
    exit 1
fi
CASE=${PWD##*/}

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

{
    printf '>>> Allpost %s\n' "$CASE"
    bash ./Allpost
    printf '>>> Alldiagnose.py\n'
    python Alldiagnose.py
    printf '>>> Allcompare.py\n'
    python Allcompare.py
    printf '>>> Allconverge.py\n'
    python Allconverge.py
    if grep -q '^overlap:' interface.yaml; then
        printf '>>> Alloverlap.py\n'
        python Alloverlap.py
    fi
    printf '>>> collect\n'
    python -m plumetools.cai2012ext.interface_case collect . \
        --output results/interface-evidence.yaml
} > output_post.log 2> error_post.log
printf '<<< Completed: post-process %s (results/interface-evidence.yaml)\n' "$CASE"
