#!/bin/bash -l
#SBATCH --mem=62G
#SBATCH -t 4-12:00:00
#SBATCH -n 8
#SBATCH -J monIBS

bamlist=~/monarch/github/F-H_2018/Raw_data/gbamlist_only_good.txt
out=~/monarch/github/F-H_2018/Data/monIBS_OnlyGood
nInd=$(wc $bamlist | awk '{print $1}')
minInd=$[$nInd*5/10]

~/angsd/angsd -bam ${bamlist} -out ${out}_clean -doMajorMinor 1 -minMapQ 20 -minQ 20 -SNP_pval 1e-12 -GL 1 -doMaf 1 -minMaf 0.05 -minInd ${minInd} -makeMatrix 1 -doCov 1 -doIBS 1 -doCounts 1 -P 8



