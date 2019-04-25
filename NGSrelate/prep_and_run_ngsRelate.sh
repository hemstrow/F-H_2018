#!/bin/bash -l
#SBATCH -t 6-24:00:00
#SBATCH --mem 20G
#SBATCH -n 4


minInd=45


~/angsd/angsd -b NA_gbamlist.txt -out relate_probs -gl 2 -domajorminor 1 -snp_pval 0.00000001 -domaf 1 -minmaf 0.05 -doGlf 3  -minQ 20 -minMapQ 20 

#~/angsd/angsd -bam NA_gbamlist.txt -GL 1 -out relate_probs -doMaf 1 -doMajorMinor 1 -SNP_pval 0.00000001 -doGeno 4 -doPost 2 -postCutoff 0.95 -minQ 20 -minMapQ 20 -minInd $minInd -minMaf 0.05

zcat relate_probs.mafs.gz | cut -f5 |sed 1d >freq

~/ngsTools/ngsRelate/ngsRelate  -g relate_probs.glf.gz -p 4 -n 91 -f freq  -O relate_out.txt
