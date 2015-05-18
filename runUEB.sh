#!/bin/bash
#SBATCH --job-name=uprun
#SBATCH --output=upr.log
#SBATCH --error=upr.err
#SBATCH --open-mode=append
#SBATCH --time=120:00:00            
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
export LD_LIBRARY_PATH=/path to /pnclibs/lib:$LD_LIBRARY_PATH
export CPPFLAGS=-I/path to /pnclibs/include
export LDFLAGS=-L/path to /pnclibs/lib 
. /rc/tools/utils/dkinit
reuse .mpich-3.1.1-slurm  
###reuse GCC-4.8
###cd $PBS_O_WORKDIR
###echo "Working directory: $PBS_O_WORKDIR"
mpirun -np 1 ./uebparpio control.dat
