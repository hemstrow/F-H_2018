#!/bin/bash -l
#SBATCH --array=1-40
#SBATCH --mem=12G
#SBATCH -t 6-12:00:00
#SBATCH -J WNAmonRFcv

# runs run_rf.R to run a random forest model multiple times, leaving one sample out each time for cross_validation

cd ~/monarch/github/F-H_2018/results/rf_cp
mkdir r_WNA_${SLURM_ARRAY_TASK_ID}
cd r_WNA_${SLURM_ARRAY_TASK_ID}

Rscript ~/monarch/github/F-H_2018/scripts/rf_refinement/run_rf.R ${SLURM_ARRAY_TASK_ID} rf_cv ${SLURM_ARRAY_TASK_ID} WNA


