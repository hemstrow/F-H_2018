#!/bin/bash -l
#SBATCH --array=3,7,12,20,27,33,38,47,50
#SBATCH --mem=12G
#SBATCH -t 6-12:00:00
#SBATCH -J ENAmonRFrefine

# runs rf_refinement.R to run a refinement of a random forest model multiple times, leaving one sample out each time for cross_validation

cd ~/monarch/github/F-H_2018/results/rf_refinement
mkdir r_ENA_${SLURM_ARRAY_TASK_ID}
cd r_ENA_${SLURM_ARRAY_TASK_ID}

Rscript ~/monarch/github/F-H_2018/scripts/rf_refinement/run_rf_refinement.R ${SLURM_ARRAY_TASK_ID} rf_refinement ${SLURM_ARRAY_TASK_ID} ENA


