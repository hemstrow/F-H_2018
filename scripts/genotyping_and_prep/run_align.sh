#!/bin/bash -l
#SBATCH -t 6-12:00:00
#SBATCH --mem=128G
#SBATCH -p bigmemm
#SBATCH -J ERIP_aln

c1=SRR1552520
ref=~/genomes/Dp_genome_v3.fasta

bwa mem $ref ${c1}.fastq | samtools view -Sb - | samtools sort - ${c1}.sort
samtools view -b ${c1}.sort.bam | samtools rmdup -s - ${c1}.sort.flt.bam
