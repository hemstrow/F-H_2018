#!/bin/bash -l
#SBATCH --array=0-49
#SBATCH --mem=20G
#SBATCH -t 6-12:00:00
#SBATCH -J NA_fSTR

module load python
source ~/.bashrc

#runs the NA data set for k 1 to 5, 10 reps of each

input=../Data/fastSTRUCTURE/NAonly
output=../fastSTRUCTURE/NA/monarch_NA

k=`expr $SLURM_ARRAY_TASK_ID + 10`
k=`expr $k / 10`
r=`expr $SLURM_ARRAY_TASK_ID % 10`
r=`expr $r + 1`

output=${output}_${r}

echo "K = $k"
echo "run $r"

python ~/bin/fastStructure/structure.py -K $k --input=${input} --output=${output} --full --format=str 
