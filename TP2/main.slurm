#!/bin/sh
#SBATCH --job-name=gravity 
#SBATCH --output=out/gravity.out 
#SBATCH --err=err/gravity.err 
#SBATCH --time=0:03:00 
#SBATCH -p mistral 
#SBATCH -N 1
#SBATCH --nodes=1                                                                                                                                                                                          
#SBATCH --ntasks-per-node=20    

module load slurm/14.03.0 compiler/gcc/6.3.0 mpi/openmpi/gcc/1.10.0-tm-mlx 

mpiexec -np 12 ./bin/main
