#!/bin/bash
#BATCH -N 1
#BATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00 
#SBATCH --partition testing

module purge
module load gcc/4.9.0 openmpi/1.10.2/intel/2016 intel/2016 gnu-openmpi/phdf5/1.8.17 intel-mkl/2017

# add here you command line to run a parallel job
export OMP_NUM_THREADS=2
mpirun -np 2 ./Main_ScalarField3d.Linux.64.mpicxx.ifort.OPTHIGH.MPI.OPENMPCC.ex params.txt


