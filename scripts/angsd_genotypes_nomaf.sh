#!/bin/bash -l
#SBATCH -t 4-12:00:00
#SBATCH --mem 24G
#SBATCH -J mgen_nmaf
#SBATCH -n 8

minInd=130
out=~/monarch/github/F-H_2018/Raw_data/nomaf_raw_genotypes
bamlist=~/monarch/github/F-H_2018/Raw_data/gbamlist.txt

~/angsd/angsd -P 8 -bam $bamlist -GL 1 -out $out -doMaf 2 -doMajorMinor 1 -SNP_pval 0.00000001 \
	      -doGeno 4 -doPost 2 -postCutoff 0.95 -minQ 20 -minMapQ 20 -minInd $minInd
