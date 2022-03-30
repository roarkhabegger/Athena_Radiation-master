#!/bin/sh
#SBATCH --partition=astro3 
#SBATCH --time=4-00:00:00 #runtime in days-hh:mm:ss
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=20 #cpus per node
#SBATCH --mem=128000 #RAM per node
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out

module load openmpi
module load hdf5
bas='/home/rhabegger/staging/'
base='ParkerInst'
run='finalA'


cd $bas/$base/$run

pwd
ls

export OMP_NUM_THREADS=20

mpiexec -n 10 astroAthena -i athinput.$base -t 95:00:00 > athena.01.log
