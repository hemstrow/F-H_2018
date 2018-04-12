#!/bin/bash -l
#SBATCH --mem=48G
#SBATCH -t 4-12:00:00
#SBATCH -n 8
#SBATCH -p bigmemm
#SBATCH -J mon_SFS

ref=~/genomes/Dp_genome_v3.fasta
out=~/monarch/github/F-H_2018/DADI/

for POP in NA HAW GUA
do
	echo "Working on ${POP}. Calling ANGSD..."
	~/angsd/angsd -P 8 -bam ~/monarch/github/F-H_2018/Raw_data/${POP}_bamlist.txt -ref $ref -anc $ref -out ${out}${POP}_folded \
		-C 50 -baq 1 -minMapQ 20 -minQ 20 -doCounts 1  -GL 2 -doSaf 1 -fold 1

	echo "calling realSFS..."
	~/angsd/misc/realSFS ${out}${POP}_folded.saf.idx -maxIter 1000 -P 8 > ${out}${POP}.smallFolded.sfs

done

