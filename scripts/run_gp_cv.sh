#!/bin/bash -l
#SBATCH --array=1-90
#SBATCH --mem=20G
#SBATCH -t 6-12:00:00
#SBATCH -J monGPcv

# runs gp_cv.R to do a leave one out cv for every sample in flt_genos

cd ../gp/
mkdir r_${SLURM_ARRAY_TASK_ID}
cd r_${SLURM_ARRAY_TASK_ID}

Rscript ~/monarch/github/F-H_2018/scripts/gp_cv.R ${SLURM_ARRAY_TASK_ID} ../BayesB_leave_one_out_cv.txt

