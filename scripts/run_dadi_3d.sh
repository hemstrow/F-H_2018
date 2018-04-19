#!/bin/bash -l
#SBATCH --array=0-49
#SBATCH --mem=4G
#SBATCH -t 6-12:00:00
#SBATCH -J 3d_dadi

#runs run_3d_mods.py with the parameters given in the parmfile, space seperated. Make sure that the correct number of array tasks are requested!

idir=
infile=
outfile=
parmfile=

#get parameters
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${file}`
IFS=$' ' read -a params <<< $line
model="${params[0]}"
maxiters="${params[1]}"
pops="${params[2]}"
fs="${params[3]}"
fp="${params[4]}"
ip="${params[5]}"
ub="${params[6]}"
lb="${params[7]}"
proj="${params[8]}"
optimizer="${params[9]}"

#run
python run_3d_mods.py $idir $infile $outfile $model $maxiters $pops $fs $fp $ip $ub $lb $proj $optimizer
