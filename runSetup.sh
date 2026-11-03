#PBS -S /bin/bash
#PBS -q shortq
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -M torresfa@my.erau.edu
#PBS -m abe
#PBS -N cai2012-nPartStudy

set -eo pipefail

source /home2/appman/spack/share/spack/setup-env.sh
source "$(spack location -i lmod)/lmod/lmod/init/bash"

module purge
module use /apps/spack/share/spack/modules/linux-rocky8-zen4
module load python/3.11.7-gcc-8.5.0-wfpoppf
# module load openmpi/5.0.2-gcc-13.2.0-4pbtkdf
module load ffmpeg/6.1.1-gcc-13.2.0-75kx2tp

source /home2/torresfa/scratch/OpenFOAM/OpenFOAM-v2512/etc/bashrc
source /home2/torresfa/scratch/OpenFOAM/rarefiedPlume/.venv/bin/activate

export PVPYTHON_PATH=/home2/torresfa/src/paraview/ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64/bin/pvpython
export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow

# #./AllrunCases --cores 192 > output.log 2> error.log
# ./AllpostCases > output_post.log 2> error_post.log

./Allsetup