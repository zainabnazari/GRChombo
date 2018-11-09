#!/bin/bash
#
#SBATCH  -N 2
#SBATCH  --ntasks-per-node=16 
#SBATCH --cpus-per-task=1
#SBATCH  --time=04:00:00 
#SBATCH  --partition long

module purge
module load gcc/4.9.0
module load gnu/openmpi/1.10.6
module load gnu-openmpi/phdf5/1.8.17
module load intel/2016
module load intel-mkl/2017
# add here you command line to run a parallel job
export OMP_NUM_THREADS=1
mpirun ./Main_ScalarField3d.Linux.64.mpicxx.ifort.DEBUG.OPT.MPI.ex params.txt
