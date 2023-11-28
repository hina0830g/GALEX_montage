#!/bin/bash
#SBATCH --cpus-per-task=5 # Specify a number of CPUs for parallel computing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --account=hamden
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=hina0830@arizona.edu
#SBATCH --array=1-3

module load python/3.9/3.9.10
source /xdisk/hamden/hina0830/test/bin/activate

CurrentCoordinates="$( sed "${SLURM_ARRAY_TASK_ID}q;d" Coordinates )"
python3 query_test.py --pair-list $CurrentCoordinates

