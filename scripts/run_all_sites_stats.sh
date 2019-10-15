#!/bin/bash -l
#SBATCH -t 10-12:00:00
#SBATCH --mem=240G
#SBATCH -p bigmemm
#SBATCH -J mon_stats

Rscript calc_all_sites_stats.R
