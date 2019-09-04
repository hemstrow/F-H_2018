#!/bin/bash -l
#SBATCH --array=1,3,6-8,15-16,19,21,23,25,31,33,38
#SBATCH --mem=12G
#SBATCH -t 6-12:00:00
#SBATCH -J WNAmonRFrefine

# runs rf_refinement.R to run a refinement of a random forest model multiple times, leaving one sample out each time for cross_validation

cd ~/monarch/github/F-H_2018/results/rf_refinement
mkdir r_WNA_${SLURM_ARRAY_TASK_ID}
cd r_WNA_${SLURM_ARRAY_TASK_ID}

Rscript ~/monarch/github/F-H_2018/scripts/rf_refinement/run_rf_refinement.R ${SLURM_ARRAY_TASK_ID} rf_refinement ${SLURM_ARRAY_TASK_ID} WNA


