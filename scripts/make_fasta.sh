#!/bin/bash -l
#SBATCH -t 6-12:00:00
#SBATCH --mem=48G
#SBATCH -J ERIP_aln

bam=SRR1552520.sort.flt.bam
outfile=D_erippus.fa

~/angsd/angsd -i $bam -doFasta 2 -doCounts 1 -out $outfile -minMapQ 10 -minQ 20 
