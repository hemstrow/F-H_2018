#!/bin/bash -l
#SBATCH --mem=48G
#SBATCH -t 6-12:00:00
#SBATCH -n 16
#SBATCH -p bigmemm
#SBATCH -J mSFSess

ref=~/genomes/Dp_genome_v3.fasta
out=~/monarch/github/F-H_2018/dadi/safs/ess

for POP in NAM HAW GUA ROT
do
	echo "Working on ${POP}. Calling ANGSD..."
	~/angsd/angsd -P 16 -bam ~/monarch/github/F-H_2018/raw_data/${POP}_ess_bamlist.txt -ref $ref -anc $ref -out ${out}${POP}_folded \
		-C 50 -baq 1 -minMapQ 20 -minQ 20 -doCounts 1  -GL 2 -doSaf 1 -fold 1

done

echo "calling realSFS for each combination..."
echo "NAM HAW..."
~/angsd/misc/realSFS ${out}NAM_folded.saf.idx ${out}HAW_folded.saf.idx -maxIter 1000 -P 16 > ${out}NH_2d.smallFolded.sfs
echo "NAM GUA..."
~/angsd/misc/realSFS ${out}NAM_folded.saf.idx ${out}GUA_folded.saf.idx -maxIter 1000 -P 16 > ${out}NG_2d.smallFolded.sfs
echo "NAM ROT..."
~/angsd/misc/realSFS ${out}NAM_folded.saf.idx ${out}ROT_folded.saf.idx -maxIter 1000 -P 16 > ${out}NR_2d.smallFolded.sfs
echo "HAW GUA..."
~/angsd/misc/realSFS ${out}HAW_folded.saf.idx ${out}GUA_folded.saf.idx -maxIter 1000 -P 16 > ${out}HG_2d.smallFolded.sfs
echo "HAW ROT..."
~/angsd/misc/realSFS ${out}HAW_folded.saf.idx ${out}ROT_folded.saf.idx -maxIter 1000 -P 16 > ${out}HR_2d.smallFolded.sfs
echo "GUA ROT..."
~/angsd/misc/realSFS ${out}GUA_folded.saf.idx ${out}ROT_folded.saf.idx -maxIter 1000 -P 16 > ${out}GR_2d.smallFolded.sfs
