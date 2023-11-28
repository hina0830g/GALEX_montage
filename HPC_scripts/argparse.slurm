#!/bin/bash
#SBATCH --cpus-per-task=48 # Specify a number of CPUs for parallel computing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --account=hamden
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=hina0830@arizona.edu
#SBATCH --array=1-3

module load python/3.9/3.9.10
source /xdisk/hamden/hina0830/venv/bin/activate

CurrentCoordinates="$( sed "${SLURM_ARRAY_TASK_ID}q;d" Coordinates )"

time python3 query_new.py --pair-list $CurrentCoordinates
time python3 pixel_nan_parallel.py --pair-list $CurrentCoordinates
time python3 mask_exptime_parallel.py --pair-list $CurrentCoordinates
time python3 mask_test.py --pair-list $CurrentCoordinates
time python3 montage_parallel.py --pair-list $CurrentCoordinates
time python3 combine.py --pair-list $CurrentCoordinates
