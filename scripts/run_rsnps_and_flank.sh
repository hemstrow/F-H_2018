#!/bin/bash -l
#SBATCH --mem=48G
#SBATCH -t 4-12:00:00
#SBATCH -p bigmemm
#SBATCH -J get_rand

echo "Getting random snps..."
Rscript ~/monarch/github/F-H_2018/scripts/get_random_snps_gap.R

echo "Finding flanking sequences..."
bash ~/monarch/github/F-H_2018/scripts/get_flanking_seqs.sh
