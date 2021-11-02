#!/bin/bash -l
#SBATCH --array=1-281
#SBATCH --mem=4G
#SBATCH -t 6-12:00:00
#SBATCH -J Align


# Aligns reads to the monarch genome and filters the result. Before running, please adjust the SBATCH arguments above, specifically the --array argument. This argument should
# be 1-the number of rows in your alignment list (aka the number of fastq files). If you have data from 200 individuals, for example, it should read --array=1-200
# The name of the job (-J), time (-t), and memory (--m) usually don't need to be changed.

list=../fastqs/fastq_list.txt # A file containing the names of the files containing data for each individual. The RA and RB names are seperated by a tab. The .fastq at the end is cut off.
read_dir=../fastqs/ # the directory containing the fastq files. Note that a final / is required after the directory (eg ~/myproject/fastqs/ not ~/myproject/fastqs).
ref=../genome/Dp_genome_v3.fasta # the path to the genome for the alignment


string="sed -n ${SLURM_ARRAY_TASK_ID}p ${list}" 
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')   
set -- $var
c1=$1
c2=$2

bwa mem $ref ${read_dir}${c1}.fastq ${read_dir}${c2}.fastq | samtools view -Sb - | samtools sort - -n -o ${c1}.sort.bam # align and sort by name
samtools fixmate -r -m ${c1}.sort.bam ${c1}.fixmate.bam # fixmate
samtools sort -o ${c1}.psort.bam ${c1}.fixmate.bam # sort by position
samtools markdup -r ${c1}.psort.bam ${c1}.markdup.bam # remove dups
samtools view -q 5 -b ${c1}.markdup.bam > ${c1}.q1.bam # remove poorly mapped
samtools sort -n -o ${c1}.namesort.bam ${c1}.q1.bam # sort by name again
samtools fixmate -m ${c1}.namesort.bam ${c1}.fixmate.bam # filter bad mates again
samtools view -f 0x2 -b ${c1}.fixmate.bam > ${c1}.flt.bam # remove improper pairs
samtools sort -o ${c1}.sort.flt.bam ${c1}.flt.bam # sort
samtools index ${c1}.sort.flt.bam # index

# clean
rm ${c1}.markdup.bam
rm ${c1}.fixmate.bam
rm ${c1}.q1.bam
rm ${c1}.psort.bam
rm ${c1}.namesort.bam
