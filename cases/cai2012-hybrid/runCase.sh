#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N hybrid-run-8M

# Select a profile with:
#   qsub -v PROFILE=verification-kn0p01 runCase.sh
set -eo pipefail

SCRIPT_DIR=${PBS_O_WORKDIR:-$(cd "$(dirname "$0")" && pwd)}
SCRIPT_DIR=$(cd "$SCRIPT_DIR" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
PROFILE=${PROFILE:-test-kn0p01-8M}
CASE_DIR=${1:-"$SCRIPT_DIR/runs/$PROFILE"}

if [ ! -f "$SCRIPT_DIR/Allrun" ] || [ ! -f "$SCRIPT_DIR/hybrid.yaml" ]; then
    printf 'FAILED: runCase.sh must be submitted from cases/cai2012-hybrid (PBS_O_WORKDIR=%s)\n' "$SCRIPT_DIR" >&2
    exit 1
fi

source "$REPO/cases/cai2012-kn-extension/hpc_env.sh"

if [ ! -f "$FOAM_USER_LIBBIN/libhybridPlumeThermo.so" ]; then
    echo ">>> Building libhybridPlumeThermo.so"
    "$REPO/applications/thermophysicalModels/Allwmake"
fi

if [ ! -f "$FOAM_USER_LIBBIN/libplumeDsmcBoundaryModels.so" ]; then
    echo ">>> Building libplumeDsmcBoundaryModels.so"
    "$REPO/applications/dsmcBoundaryModels/Allwmake"
fi

export PATH="$FOAM_APPBIN:$PATH"
echo ">>> Running profile $PROFILE in $CASE_DIR"
bash "$SCRIPT_DIR/Allrun" "$CASE_DIR" all \
    --config "$SCRIPT_DIR/hybrid.yaml" --profile "$PROFILE"
