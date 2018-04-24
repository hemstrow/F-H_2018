#!/bin/bash -l
#SBATCH --mem=62G
#SBATCH -t 4-12:00:00
#SBATCH -J monGLF
#SBATCH -n 8

bamlist=~/monarch/github/F-H_2018/NGSadmix/gbamlist_only_good_sorted.txt
out=~/monarch/github/F-H_2018/NGSadmix/monGLF
nInd=$(wc $bamlist | awk '{print $1}')
minInd=$[$nInd*5/10]

~/angsd/angsd -bam ${bamlist} -out ${out} -doMajorMinor 1 -minMapQ 20 -minQ 20 -SNP_pval 1e-12 -doMaf 1 -minMaf 0.05 -minInd ${minInd} -doGlf 2 -P 8 -GL 1

