#!/bin/bash
#SBATCH -t 10-12:00:00
#SBATCH --mem=200G
#SBATCH -p bigmemm
#SBATCH -J mon_par_prep

module load R/4.1.0


Rscript Data_prep.R
