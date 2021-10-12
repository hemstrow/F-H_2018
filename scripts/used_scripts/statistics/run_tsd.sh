#!/bin/bash
#SBATCH --array=0-12
#SBATCH -t 10-12:00:00
#SBATCH --mem=60G
#SBATCH -J par_tsd

pops=(HAW GUA ROT SAI SAM FIJ NCA NOR QLD NSW VIC NZL NAM)

Rscript calc_tsd.R allsites_paralog_fix_snpR.RDS ${pops[${SLURM_ARRAY_TASK_ID}]}

