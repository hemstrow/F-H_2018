#!/bin/bash -l
#SBATCH --mem=96G
#SBATCH -t 4-12:00:00
#SBATCH -n 8
#SBATCH -p bigmemm
#SBATCH -J 3d_unf_SFS

#takes three idx files from -doSaf and generates a 3dsfs via realsfs

idx1=../DADI/NA_folded.saf.idx
idx2=../DADI/HAW_folded.saf.idx
idx3=../DADI/GUA_folded.saf.idx
out=../DADI/NHG_folded

~/angsd/misc/realSFS $idx1 $idx2 $idx3 -maxIter 1000 -P 8 > ${out}.3dsmall.sfs
