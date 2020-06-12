#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --partition=skl_usr_prod
#SBATCH --mem=182000
#SBATCH -A ict20_hep
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=3
# add here you command line to run a parallel job
export OMP_NUM_THREADS=4
mpirun ./../../Main_ScalarField3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ex params.txt
