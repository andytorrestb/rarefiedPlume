#PBS -S /bin/bash
#PBS -q longq
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=192
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N post-r-4d

set -eo pipefail
trap 'status=$?; printf "FAILED: exit=%d line=%d command=%s\n" "$status" "$LINENO" "$BASH_COMMAND" >&2; exit "$status"' ERR

printf '>>> Starting: load Spack environment\n'
source /home2/appman/spack/share/spack/setup-env.sh
printf '<<< Completed: load Spack environment\n'

printf '>>> Starting: initialize Lmod\n'
source "$(spack location -i lmod)/lmod/lmod/init/bash"
printf '<<< Completed: initialize Lmod\n'

printf '>>> Starting: purge modules\n'
module purge
printf '<<< Completed: purge modules\n'
printf '>>> Starting: add Spack module path\n'
module use /apps/spack/share/spack/modules/linux-rocky8-zen4
printf '<<< Completed: add Spack module path\n'
printf '>>> Starting: load Python module\n'
module load python/3.11.7-gcc-8.5.0-wfpoppf
printf '<<< Completed: load Python module\n'
printf '>>> Starting: load OpenMPI module\n'
module load openmpi/5.0.2-gcc-13.2.0-4pbtkdf
printf '<<< Completed: load OpenMPI module\n'
printf '>>> Starting: load FFmpeg module\n'
module load ffmpeg/6.1.1-gcc-13.2.0-75kx2tp
printf '<<< Completed: load FFmpeg module\n'

printf '>>> Starting: load OpenFOAM environment\n'
if source /home2/torresfa/scratch/OpenFOAM/OpenFOAM-v2512/etc/bashrc; then
    printf '<<< Completed: load OpenFOAM environment\n'
else
    status=$?
    printf 'FAILED: loading OpenFOAM environment returned exit %d\n' "$status" >&2
    exit "$status"
fi
printf '>>> Starting: activate Python virtual environment\n'
source /home2/torresfa/scratch/OpenFOAM/rarefiedPlume/.venv/bin/activate
printf '<<< Completed: activate Python virtual environment\n'

printf '>>> Starting: configure ParaView environment\n'
export PVPYTHON_PATH="${PVPYTHON_PATH:-/home2/torresfa/src/paraview/ParaView-6.2.0-MPI-Linux-Python3.12-x86_64/bin/pvpython}"
export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow
printf '<<< Completed: configure ParaView environment\n'

# printf '>>> Starting: run Allmesh\n'
# ./Allmesh
# printf '<<< Completed: run Allmesh\n'
# bash ./Allrun -np 192 > output.log 2> error.log
./Allpost > output_post.log 2> error_post.log
