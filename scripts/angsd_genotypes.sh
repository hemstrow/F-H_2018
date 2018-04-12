#!/bin/bash -l
#SBATCH -t 4-12:00:00
#SBATCH --mem 20G


minInd=130

~/angsd/angsd -bam gbamlist.txt -GL 1 -out genotypes -doMaf 2 -doMajorMinor 1 -SNP_pval 0.00000001 -doGeno 4 -doPost 2 -postCutoff 0.95 -minQ 20 -minMapQ 10 -minInd $minInd -minMaf 0.05
