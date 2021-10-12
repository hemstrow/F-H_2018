#!/bin/bash -l
#SBATCH -t 7-12:00:00
#SBATCH --mem=128G
#SBATCH -p bigmemm
#SBATCH -J fst

Rscript run_fst.R

