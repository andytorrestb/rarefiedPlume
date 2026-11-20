#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=32
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N ugf-build
#
# Build the uniGasFoam fork at the commit study.yaml pins, and record it.
# Submit from here:
#
#     cd cases/cai2012-unigas && qsub runBuild.sh
#     qsub -v UNIGAS_DIR=/path/to/uniGasFoam runBuild.sh     # another checkout
#
#   1. check out meta.unigas_commit in $UNIGAS_DIR (default
#      /home2/torresfa/scratch/OpenFOAM/uniGasFoam); a dirty checkout is refused
#   2. ./Allwmake with this job's cores (it copies the sources into
#      $WM_PROJECT_USER_DIR and builds there)
#   3. write the commit to $FOAM_USER_APPBIN/uniGasFoam.commit, which every
#      case's Allrun compares with its own meta.unigas_commit before running
set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

if [ -z "${PBS_O_WORKDIR:-}" ]; then
    printf 'FAILED: PBS_O_WORKDIR is not set. Submit with qsub from cases/cai2012-unigas\n' >&2
    exit 1
fi
cd "$PBS_O_WORKDIR"
[ -f study.yaml ] || { printf 'FAILED: %s has no study.yaml\n' "$PWD" >&2; exit 1; }
source ./hpc_env.sh

UNIGAS_DIR=${UNIGAS_DIR:-/home2/torresfa/scratch/OpenFOAM/uniGasFoam}
COMMIT=$(python -c "import yaml; print(yaml.safe_load(open('study.yaml'))['meta']['unigas_commit'])")
case "$COMMIT" in
    ''|None|pending-*) printf 'FAILED: study.yaml meta.unigas_commit is not set (%s)\n' "$COMMIT" >&2; exit 1 ;;
esac
[ -d "$UNIGAS_DIR/.git" ] || { printf 'FAILED: no git checkout at %s\n' "$UNIGAS_DIR" >&2; exit 1; }

printf '>>> Starting: check out uniGasFoam %s in %s\n' "$COMMIT" "$UNIGAS_DIR"
cd "$UNIGAS_DIR"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    printf 'FAILED: %s has uncommitted changes; refusing to switch commits\n' "$UNIGAS_DIR" >&2
    exit 1
fi
git fetch --all --quiet || true
git checkout --quiet "$COMMIT"
BUILT=$(git rev-parse HEAD)
printf '<<< Completed: check out %s\n' "$BUILT"

printf '>>> Starting: build (./Allwmake %s)\n' "${PBS_NP:-8}"
./Allwmake "${PBS_NP:-8}"
command -v uniGasFoam >/dev/null || { printf 'FAILED: uniGasFoam is not on PATH after the build\n' >&2; exit 1; }
printf '%s\n' "$BUILT" > "$FOAM_USER_APPBIN/uniGasFoam.commit"
printf '<<< Completed: build; %s/uniGasFoam.commit = %s\n' "$FOAM_USER_APPBIN" "$BUILT"
