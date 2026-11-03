#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -m abe
#PBS -N run-r-4d

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO_ROOT/cases/cai2012-kn-extension/hpc_env.sh"

# printf '>>> Starting: run Allmesh\n'
# ./Allmesh
# printf '<<< Completed: run Allmesh\n'
bash ./Allrun -np 192 > output.log 2> error.log
# ./AllpostCases > output_post.log 2> error_post.log
