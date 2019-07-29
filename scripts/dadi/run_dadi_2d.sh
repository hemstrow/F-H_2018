#!/bin/bash -l
#SBATCH --array=1-4200
#SBATCH --mem=5G
#SBATCH -t 6-12:00:00
#SBATCH -J NHmondadi

#runs run_2d_mods.py with the parameters given in the parmfile, space seperated. Make sure that the correct number of array tasks are requested!

module load bio
module load matplotlib
module load dadi

idir=~/monarch/github/F-H_2018/data/dadi_inputs/
infile=dadi_10kgap_snps.txt
outfile=NH_r4_adjusted_portik.txt
parmfile=~/monarch/github/F-H_2018/dadi/parmfiles/NH_r4_adjusted_portik.txt
pyscript=~/monarch/github/F-H_2018/scripts/dadi/run_2d_mods.py



#get parameters
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${parmfile}`
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

echo "Beginning. Call:"
echo "$pyscript $idir $infile $outfile $model $maxiters $pops $fs $fp $ip $ub $lb $proj $optimizer"

python $pyscript $idir $infile run_${SLURM_ARRAY_TASK_ID}_${outfile} $model $maxiters $pops $fs $fp $ip $ub $lb $proj $optimizer
