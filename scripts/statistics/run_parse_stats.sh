#!/bin/bash -l
#SBATCH -t 7-12:00:00
#SBATCH --mem=128G
#SBATCH -p bigmemm
#SBATCH -J fst


module load R/4.1.0
cd ~/monarch/github/F-H_2018/
Rscript ./scripts/used_scripts/statistics/parse_stats_outputs.R

