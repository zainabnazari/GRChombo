#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --partition=skl_usr_prod
#SBATCH --mem=180000
#SBATCH -A ict20_hep
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
# add here you command line to run a parallel job
export OMP_NUM_THREADS=2
mpirun ./ Main_ScalarField3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.ex params.txt
