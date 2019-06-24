#!/bin/bash
#SBATCH -p GPU-All
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4

#export LD_LIBRARY_PATH=/sharedext4/lib64/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/sharedext4/lib64/intel64_lin:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/sharedext4/lib64/mpich/lib/:$LD_LIBRARY_PATH
source /sharedext4/softwares/source/vasp_gpu.env
/sharedext4/softwares/vasp/vasp.5.4.4.gpu/build/gpu/vasp
