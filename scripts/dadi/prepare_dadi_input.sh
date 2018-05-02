#!/bin/bash -l
#SBATCH --mem=96G
#SBATCH -t 60-12:00:00
#SBATCH -p bigmemm
#SBATCH -J prepDADI

echo "Getting random snps..."
Rscript ~/monarch/github/F-H_2018/scripts/dadi/filt_and_rsnps_gap.R

echo "Finding flanking sequences..."
bash ~/monarch/github/F-H_2018/scripts/dadi/get_flanking_seqs.sh

echo "Formatting for dadi..."
Rscript ~/monarch/github/F-H_2018/scripts/dadi/convert_dadi.R
