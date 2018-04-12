#!/bin/bash -l
#SBATCH --mem=48G
#SBATCH -t 6-12:00:00
#SBATCH -J get_flank

#this function takes an ac_formated snp dataset (or just something where the first column is the chr and the second is the pos) and finds the flanking sequence in both a reference and ancestor genome.
#This is needed for dadi, at least for an unfolded spectra!

snp_ac=~/monarch/github/F-H_2018/Data/rand_gap_snps.txt
anc=~/monarch/Erp_BRA_16005/D_erippus.fa
ref=~/genomes/Dp_genome_v3.fasta

wc=$(wc -l $snp_ac | awk '{print $1}')

x=1

q=""

echo "Creating samtools query..."

while [ $x -le $wc ]
do

	if (( $x % 5 == 0 ))           # no need for brackets
	then
		echo "Progress: $x out of ${wc}."
	fi

        string="sed -n ${x}p $snp_ac"
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1, $2}')
        set -- $var
        chr=$1
        pos=$2

	end=`expr $pos + 1`
	start=`expr $pos - 1`

	q="$q ${chr}:${start}-${end}"
        x=$(( $x + 1 ))
done

echo "Querying reference fasta."
samtools faidx $ref $q | awk -v ADD="---" -v PATTERN="^>" 'L { print L; if((L ~ PATTERN) && ($0 ~ PATTERN)) print ADD }; { L=$0 } END { print L }' | grep -v "^>" > ref_flank.txt

echo "Querying ancestral fasta."
samtools faidx $anc $q | awk -v ADD="---" -v PATTERN="^>" 'L { print L; if((L ~ PATTERN) && ($0 ~ PATTERN)) print ADD }; { L=$0 } END { print L }' | grep -v "^>" > anc_flank.txt
echo "Done."
