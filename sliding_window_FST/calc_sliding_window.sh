#!/bin/bash
#SBATCH -t 6-24:00:00
#SBATCH --mem=20G
#SBATCH -J ENA_WNA_fst

# copy this script into a new directory under monarchs
# edit the paths for angsd and realSFS. RealSFS is in the angsd/misc folder
# edit the paths to the bamlists, they are in the Raw_data folder of the github repo
# don't run this in the github, run it in the general monarch folder
# angsd is in ~mlyjones/bin/angsd_dir

# to run this, type sbatch calc_sliding_window.sh

angsdpath=~/angsd/angsd
sfspath=~/angsd/misc/realSFS

#this is with 2pops
#first calculate per pop saf for each population
${angsdpath} -b ../Raw_data/ENA_bamlist.txt -anc /home/hemstrow/genomes/Dp_genome_v3.fasta -out ENA -dosaf 1 -gl 1 -fold 1
${angsdpath} -b ../Raw_data/WNA_bamlist.txt -anc /home/hemstrow/genomes/Dp_genome_v3.fasta -out WNA -dosaf 1 -gl 1 -fold 1

#calculate the 2dsfs prior
ll ${sfspath} ENA.saf.idx WNA.saf.idx > ENA.WNA.ml

#prepare the fst for easy window analysis etc
${sfspath} fst index ENA.saf.idx WNA.saf.idx -sfs ENA.WNA.ml -fstout ENA.WNA

#get the global estimate
${sfspath} fst stats ENA.WNA.fst.idx > ENA_WNA_global_fst.txt

#below is not tested that much, but seems to work
${sfspath} fst stats2 ENA.WNA.fst.idx -win 50000 -step 10000 > ENA_WNA_sliding_window.txt
