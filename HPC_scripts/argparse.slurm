#!/bin/bash
#SBATCH --cpus-per-task=16 # Specify a number of CPUs for parallel computing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --account=hamden
#SBATCH --partition=standard
#SBATCH --mail-type=END
#SBATCH --mail-user=hina0830@arizona.edu
#SBATCH --array=1-4

module load python/3.9
source /xdisk/hamden/hina0830/venv39/bin/activate

input_path="/xdisk/hamden/hina0830/venv39/input/coordinates1"
CurrentCoordinates="$( sed "${SLURM_ARRAY_TASK_ID}q;d" "$input_path" )"

time python3 query.py --pair-list $CurrentCoordinates
time python3 preprocessing.py --pair-list $CurrentCoordinates
time python3 object_detection.py --pair-list $CurrentCoordinates

module load julia
export JULIA_NUM_THREADS=10

srun julia pointsource_infill.jl $SLURM_ARRAY_TASK_ID $input_path
time python3 montage.py --pair-list $CurrentCoordinates
