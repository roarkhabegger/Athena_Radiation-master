#!/bin/sh
#SBATCH --partition=univ2 
#SBATCH --time=4-00:00:00 #runtime in days-hh:mm:ss
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=20 #cpus per node
#SBATCH --mem=128000 #RAM per node
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out

module load openmpi
module load hdf5
bas='/home/rhabegger/staging/'
base='ParkerInst'
run='aU'


cd $bas/$base/$run

pwd
ls

./univAthena -i athinput.$base -m 320 > load.log
mpiexec -n 320 univAthena -i athinput.$base -t 95:00:00 > athena.01.log
