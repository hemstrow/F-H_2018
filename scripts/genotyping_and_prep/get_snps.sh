#!/bin/bash
#SBATCH --array=0-12
#SBATCH -t 10-12:00:00
#SBATCH --mem=2G
#SBATCH -J mp_snps

pops=(HAW GUA ROT SAI SAM FIJ NCA NOR QLD NSW VIC NZL NAM) ### array containing population names ###
ref=~/genomes/Dp_genome_v3.fasta ### reference fasta file used in alignment ###

### Calculate paralog probabilities and get a list of paralogous loci  ###

# find the pop
pop=${pops[${SLURM_ARRAY_TASK_ID}]}

#angsd -bam bamlists/${pop}_bamlist_fixed.txt -out results_snp/$pop -ref $ref -GL 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 20 -minQ 20
#gunzip results_snp/${pop}*.gz
cut -d$'\t' -f1-2  results_snp/${pop}.mafs | sed 1d > results_snp/${pop}.snp.pos

