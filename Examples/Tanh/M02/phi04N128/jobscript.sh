#!/bin/bash -l

#SBATCH -A dp016
#SBATCH -p cosma7
#SBATCH --nodes 10
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=2
#SBATCH -o output_file.%J.out
#SBATCH -e error_file.%J.err
#SBATCH -t 24:00:00
#SBATCH -J N128Tanh
#SBATCH --exclusive
#SBATCH --mail-type=ALL                          # notifications for job done & fail
#SBATCH --mail-user=katy.clough@physics.ox.ac.uk

#  The above runs a job on 4 nodes, with each node only executing a single
#  rank, so you have 28 or 16 threads per node depending on the compute nodes.

module purge
#load the modules used to build your program.
module load intel_comp/2019 intel_mpi/2019 parallel_hdf5/1.10.3

export OMP_NUM_THREADS=2

# Run the program
mpirun -np $SLURM_NTASKS ./../../Main_ScalarField3d.Linux.64.mpiicpc.ifort.DEBUG.OPT.MPI.ex params.txt
