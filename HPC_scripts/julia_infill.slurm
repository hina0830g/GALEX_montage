#!/bin/bash
#SBATCH --cpus-per-task=94 # Specify a number of CPUs for parallel computing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --account=hamden
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=hina0830@arizona.edu

export JULIA_NUM_THREADS=10
module load julia

julia pointsource_infill.jl
