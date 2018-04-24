#!/bin/bash -l
#SBATCH --mem=96G
#SBATCH -t 40-12:00:00
#SBATCH -n 8
#SBATCH -p bigmemm
#SBATCH -J 3d_unf_SFS

#takes three idx files from -doSaf and generates a 3dsfs via realsfs

idx1=../DADI/safs/NA.saf.idx
idx2=../DADI/safs/HAW.saf.idx
idx3=../DADI/safs/GUA.saf.idx
out=../DADI/safs/NHG_unfolded

~/angsd/misc/realSFS $idx1 $idx2 $idx3 -maxIter 100 -P 8 > ${out}.3dsmall.sfs
