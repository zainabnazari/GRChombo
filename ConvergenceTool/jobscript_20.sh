#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --constraint="ivybridge-ep"
#SBATCH --time=24:00:00 
#SBATCH --partition long

module purge
module load gcc/4.9.0 openmpi/1.10.2/intel/2016 intel/2016 gnu-openmpi/phdf5/1.8.17 intel-mkl/2017

# add here you command line to run a parallel job
export OMP_NUM_THREADS=1
mpirun ./ConvergenceTool3d.Linux.64.mpicxx.ifort.DEBUG.OPT.MPI.OPENMPCC.ex params.txt
