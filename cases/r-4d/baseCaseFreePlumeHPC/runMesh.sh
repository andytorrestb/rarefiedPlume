#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -N mesh-r-4d

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO_ROOT/cases/cai2012-kn-extension/hpc_env.sh"

printf '>>> Starting: run Allmesh\n'
./Allmesh
printf '<<< Completed: run Allmesh\n'
# #./AllrunCases --cores 192 > output.log 2> error.log
# ./AllpostCases > output_post.log 2> error_post.log
