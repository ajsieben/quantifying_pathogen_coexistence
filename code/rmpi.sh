#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=3-00:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=asieben@uwyo.edu
#SBATCH --job-name=pathogen_coexistence

# Change to the relevant working directory
cd /project/coexistence/pathogen_coexistence/neutral_general

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 openmpi/3.1.0 r-rmpi/0.6-9-r353-py27

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet rmpi.R  rmpi.log
