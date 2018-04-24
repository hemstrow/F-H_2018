#!/bin/bash -l
#SBATCH --array=0-129
#SBATCH --mem=20G
#SBATCH -t 6-12:00:00
#SBATCH -J f_NGSa

#runs the full data set for k 1 to 13, 10 reps of each

input=~/monarch/github/F-H_2018/NGSadmix/monGLF.beagle.gz
output=~/monarch/github/F-H_2018/NGSadmix/full/out
mmaf=0.05

k=`expr $SLURM_ARRAY_TASK_ID + 10`
k=`expr $k / 10`
r=`expr $SLURM_ARRAY_TASK_ID % 10`
r=`expr $r + 1`

output=${output}_K${k}_r${r}

echo "K = $k"
echo "run $r"

~/angsd/misc/NGSadmix -likes $input -K $k -minMaf $mmaf -o $output 

