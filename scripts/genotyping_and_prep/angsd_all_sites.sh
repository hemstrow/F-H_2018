#!/bin/bash -l
#SBATCH -t 4-12:00:00
#SBATCH --mem 60G
#SBATCH -J a.s.p.genotypes
#SBATCH -p high

minInd=130

angsd -bam run_gbamlist.txt -GL 1 -out all_sites_paralogs -doMaf 2 -doMajorMinor 1  -doGeno 4 -doPost 2 -postCutoff 0.95 -minQ 20 -minMapQ 20 -minInd $minInd -rf ../github/F-H_2018/results/paralogs/selected_clean_regions.txt

gunzip all_sites_paralogs.geno.gz
