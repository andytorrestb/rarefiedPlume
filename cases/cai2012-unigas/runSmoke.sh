#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=16
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-smoke
#
# Smoke test: the whole chain on a few cores, before the 192-rank pilot asks
# for a full node for 120 h. NOT a production or pilot result. Submit from here,
# AFTER runBuild.sh has finished:
#
#     cd cases/cai2012-unigas && msub runSmoke.sh
#     msub -v SMOKE_KN=Kn1,SMOKE_TRANSITS=0.02 runSmoke.sh     # another Kn / length
#     msub -l nodes=1:ppn=32 -v SMOKE_TRANSITS=1,SMOKE_TAG=1T runSmoke.sh
#                                   # -> Smoke/Kn0p1_np1x-1T; long enough for the
#                                   #    hybrid decomposition (every 100 steps)
#
# Then the post chain, including the Cai comparison: runSmokePost.sh.
#
#   1. python generate_cases.py --pilot $SMOKE_KN, moved to Smoke/<case>
#      (outside the manifest; the real pilot's Pilot/ is left alone)
#   2. ./Allmesh                                      the production mesh, serial
#   3. ./Allrun -np $PBS_NP --stop-at-transits X      build-commit check, decompose,
#                                                     startup leg, reconstruct, convert
#   4. cell centres + the pilot report (cost projection from this short run)
#
# Default Kn0p1: the hybrid USP <-> DSMC switching is exercised there, which
# Kn1 (mostly DSMC) would not. Read Smoke/<case>/results/pilot_report.md and
# Smoke/<case>/log.uniGasFoam.
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with msub from cases/cai2012-unigas\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f study.yaml ] && [ -f generate_cases.py ] \
    || { printf 'FAILED: %s is not the family directory\n' "$PWD" >&2; exit 1; }
source ./hpc_env.sh

SMOKE_KN=${SMOKE_KN:-Kn0p1}
SMOKE_TRANSITS=${SMOKE_TRANSITS:-0.05}
SMOKE_NP=${SMOKE_NP:-${PBS_NP:-16}}
SMOKE_TAG=${SMOKE_TAG:-}            # optional suffix: Smoke/<case>-<tag>, to keep runs apart
CASE_NAME="${SMOKE_KN}_np1x"
SMOKE_DIR="$PWD/Smoke/$CASE_NAME${SMOKE_TAG:+-$SMOKE_TAG}"

if [ -d "$SMOKE_DIR" ]; then
    if compgen -G "$SMOKE_DIR/processor*" >/dev/null || [ -f "$SMOKE_DIR/log.uniGasFoam" ]; then
        printf 'FAILED: %s holds solver output; remove it by hand to re-run the smoke test\n' "$SMOKE_DIR" >&2
        exit 1
    fi
    rm -rf "$SMOKE_DIR"
fi
if [ -e "Pilot/$CASE_NAME" ]; then
    printf 'FAILED: Pilot/%s exists; staging would replace it. Move it aside first\n' "$CASE_NAME" >&2
    exit 1
fi

printf '>>> Starting: stage the smoke case %s\n' "$CASE_NAME"
STAGED=$(python generate_cases.py --pilot "$SMOKE_KN")
mkdir -p Smoke
mv "$STAGED" "$SMOKE_DIR"
rmdir Pilot 2>/dev/null || true
printf 'SMOKE: %s, %s ranks, %s transits. Never a production or pilot result.\n' \
    "$CASE_NAME" "$SMOKE_NP" "$SMOKE_TRANSITS" > "$SMOKE_DIR/SMOKE"
printf '<<< Completed: stage %s\n' "$SMOKE_DIR"
cd "$SMOKE_DIR"

printf '>>> Starting: mesh\n'
START=$(date +%s)
bash ./Allmesh > output_mesh.log 2> error_mesh.log
printf '<<< Completed: mesh (%d s)\n' "$(( $(date +%s) - START ))"

printf '>>> Starting: run %s ranks, %s transits\n' "$SMOKE_NP" "$SMOKE_TRANSITS"
START=$(date +%s)
bash ./Allrun -np "$SMOKE_NP" --stop-at-transits "$SMOKE_TRANSITS" >> output.log 2> error.log
printf '<<< Completed: run (%d s)\n' "$(( $(date +%s) - START ))"

printf '>>> Starting: cell centres and the pilot report\n'
postProcess -func writeCellCentres -latestTime > log.writeCellCentres 2>&1
python -m plumetools.cai2012ugf.pilot .
printf '<<< Completed: %s/results/pilot_report.md\n' "$PWD"
