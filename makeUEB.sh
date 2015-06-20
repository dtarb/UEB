#!/bin/bash
###------------This script works on USU HPC; may need changes on how libraries are loaded on a different platform
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
#
mpicxx -std=c++0x -g main.cpp canopy.cpp matrixnvector.cpp ncfunctions.cpp snowdgtv.cpp snowdv.cpp snowxv.cpp uebdecls.cpp uebinputs.cpp -o uebparpio -I/path to /pnclibs/include -L/path to /pnclibs/lib -lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lpnetcdf -ldl 
#On ms windows use the file 'ncfunctions_mswin.cpp' instead of 'ncfunctions.cpp'  
