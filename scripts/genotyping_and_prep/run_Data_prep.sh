#!/bin/bash
#SBATCH -t 10-12:00:00
#SBATCH --mem=200G
#SBATCH -p bigmemm
#SBATCH -J mon_par_prep

Rscript Data_prep.R
