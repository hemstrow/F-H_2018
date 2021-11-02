#!/bin/bash
#SBATCH --array=0-12
#SBATCH -t 10-12:00:00
#SBATCH --mem=40G
#SBATCH -J mon_paralog

pops=(HAW GUA ROT SAI SAM FIJ NCA NOR QLD NSW VIC NZL NAM) ### array containing population names ###
ref=~/genomes/Dp_genome_v3.fasta ### reference fasta file used in alignment ###

### Calculate paralog probabilities and get a list of paralogous loci  ###

# find the pop
pop=${pops[${SLURM_ARRAY_TASK_ID}]}

samtools mpileup -b bamlists/${pop}_bamlist_fixed.txt -l results_snp/${pop}.snp.pos -f ${ref} > results_paralogs/${pop}.depth
~/bin/ngsParalog/ngsParalog calcLR -infile results_paralogs/${pop}.depth > results_paralogs/${pop}.paralogs
