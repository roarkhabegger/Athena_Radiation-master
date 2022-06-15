#!/bin/sh
#SBATCH --partition=univ2,astro3 
#SBATCH --time=0-04:00:00 #runtime in days-hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20 #cpus per node
#SBATCH --mem=128000 #RAM per node
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out
module load openmpi
module load hdf5
module load gcc
module load netcdf-c++4

make clean
make clean
sh my_conf_perturb.sh
make -j 20
