#!/bin/bash -l
#SBATCH -t 4-12:00:00
#SBATCH --mem 4G
#SBATCH -J chr_lengths
#SBATCH -p high


cat  ../../genomes/Dp_genome_v3.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > chr_lengths.txt
