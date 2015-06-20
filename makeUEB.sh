#!/bin/bash
#SBATCH --job-name=ump
#SBATCH --output=ump.log
#SBATCH --error=ump.err
#SBATCH --open-mode=append
#SBATCH --time=00:10:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 1 
export LD_LIBRARY_PATH=/path to /pnclibs/lib:$LD_LIBRARY_PATH
export CPPFLAGS=-I/path to /pnclibs/include
export LDFLAGS=-L/path to /pnclibs/lib 
. /rc/tools/utils/dkinit
###reuse OpenMPI
reuse .mpich-3.1.1-slurm 
###reuse GCC-4.8
###cd $PBS_O_WORKDIR
###echo "Working directory: $PBS_O_WORKDIR"
mpicxx -std=c++0x -g main.cpp canopy.cpp matrixnvector.cpp nctestfunc.cpp snowdgtv.cpp snowdv.cpp snowxv.cpp uebdecls.cpp uebinputs.cpp -o uebparpio -I/path to /pnclibs/include -L/path to /pnclibs/lib -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lpnetcdf -ldl 
