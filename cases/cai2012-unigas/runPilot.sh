#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-pilot
#
# The pilot (spec §8): Kn = 0.001 on its PRODUCTION mesh, 192 ranks, for two
# domain transits -- not to steady state. It measures whether 3-D production on
# uniGasFoam is affordable: the time step and what sets it, the DSMC/USP split,
# the collision spacing in DSMC cells, the particle count, the rank imbalance,
# the disk per write and a projected wall time. Submit from here:
#
#     cd cases/cai2012-unigas && qsub runPilot.sh
#     qsub -v PILOT_KN=Kn0p01,PILOT_TRANSITS=1 runPilot.sh    # another Kn / length
#
#   1. python generate_cases.py --pilot $PILOT_KN   Pilot/<case>, outside the manifest
#   2. ./Allmesh                                      the production mesh
#   3. ./Allrun --stop-at-transits $PILOT_TRANSITS    reconstructs and converts at the end
#   4. cell centres, then python -m plumetools.cai2012ugf.pilot
#
# Read Pilot/<case>/results/pilot_report.md. Going ahead is the owner's call.
# A pilot that ran out of walltime is resubmitted BY HAND: step 1 refuses a
# Pilot/<case> that holds output, so resume with the case's own Allrun instead.
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from cases/cai2012-unigas\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f study.yaml ] && [ -f generate_cases.py ] \
    || { printf 'FAILED: %s is not the family directory\n' "$PWD" >&2; exit 1; }
source ./hpc_env.sh

PILOT_KN=${PILOT_KN:-Kn0p001}
PILOT_TRANSITS=${PILOT_TRANSITS:-2}

printf '>>> Starting: stage the pilot %s\n' "$PILOT_KN"
CASE_DIR=$(python generate_cases.py --pilot "$PILOT_KN")
printf '<<< Completed: stage %s\n' "$CASE_DIR"
cd "$CASE_DIR"

printf '>>> Starting: mesh\n'
bash ./Allmesh > output_mesh.log 2> error_mesh.log
printf '<<< Completed: mesh\n'

printf '>>> Starting: run %s transits\n' "$PILOT_TRANSITS"
START=$(date +%s)
bash ./Allrun --stop-at-transits "$PILOT_TRANSITS" >> output.log 2> error.log
printf '<<< Completed: run (%d s)\n' "$(( $(date +%s) - START ))"

printf '>>> Starting: cell centres and the pilot report\n'
postProcess -func writeCellCentres -latestTime > log.writeCellCentres 2>&1
python -m plumetools.cai2012ugf.pilot .
printf '<<< Completed: %s/results/pilot_report.md\n' "$PWD"
