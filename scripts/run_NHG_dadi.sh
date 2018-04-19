#!/bin/bash -l
#SBATCH --array=0-19
#SBATCH --mem=16G
#SBATCH -t 6-12:00:00
#SBATCH -J NHG_dadi

module load python
module load numpy
module load dadi
module load scipy
module load matplotlib

output=../DADI/results/NA_HAW_GUA_${SLURM_ARRAY_TASK_ID}

python run_NA_HAW_GUA.py > $output
