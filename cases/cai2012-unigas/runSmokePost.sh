#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=8
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-smoke-post
#
# The per-case post chain (Allpostugf: convert, Allpost, Alldiagnose, Allaudit,
# Allconverge, Allverdict, Allcompare vs Cai) on a smoke case. The smoke cases
# are outside the manifest, so AllpostCases / runPost.sh cannot reach them.
# Proves the plumbing only: a smoke run is far too short for a meaningful
# comparison or stationarity verdict. Submit from here:
#
#     cd cases/cai2012-unigas && msub runSmokePost.sh
#     msub -v SMOKE_CASE=Kn0p1_np1x-1T runSmokePost.sh
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with msub from cases/cai2012-unigas\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f study.yaml ] || { printf 'FAILED: %s is not the family directory\n' "$PWD" >&2; exit 1; }
source ./hpc_env.sh

SMOKE_CASE=${SMOKE_CASE:-Kn0p1_np1x}
[ -f "Smoke/$SMOKE_CASE/SMOKE" ] || { printf 'FAILED: Smoke/%s is not a smoke case\n' "$SMOKE_CASE" >&2; exit 1; }
cd "Smoke/$SMOKE_CASE"

printf '>>> Starting: Allpostugf on %s\n' "$SMOKE_CASE"
START=$(date +%s)
bash ./Allpostugf > output_post.log 2> error_post.log
printf '<<< Completed: Allpostugf (%d s); see %s/results\n' "$(( $(date +%s) - START ))" "$PWD"
