# cases/cai2012-unigas/hpc_env.sh -- the cluster environment, shared.
#
# SOURCED by every job script of this family. The settings live in ONE place,
# the kn-extension's hpc_env.sh (Spack, Lmod, OpenMPI, OpenFOAM v2512, the
# repository venv, ParaView): this family runs on the same machine.

hpc_env_here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$hpc_env_here/../cai2012-kn-extension/hpc_env.sh" || return 1
unset hpc_env_here
