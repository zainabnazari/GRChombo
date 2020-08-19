#!/bin/bash
#SBATCH --nodes=11
#SBATCH --ntasks-per-node=6
#SBATCH --partition=skl_usr_prod
#SBATCH --mem=182000
#SBATCH -A ict20_hep
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
# add here you command line to run a parallel job
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun ./../../Main_ScalarField3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ex params.txt
