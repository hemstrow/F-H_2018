#!/bin/bash -l
#SBATCH --array=1-800
#SBATCH --mem=2G
#SBATCH -t 12-12:00:00
#SBATCH -J MonDadi2

#runs run_3d_mods.py with the parameters given in the parmfile, space seperated. Make sure that the correct number of array tasks are requested!

module load numpy
module load scipy
module load dadi
module load matplotlib

idir=~/monarch/github/F-H_2018/data/dadi_inputs/
infile=dadi_10kgap_snps.txt
outfile=1st_pass_optim_HGR_new_bounds.txt
parmfile=~/monarch/github/F-H_2018/dadi/parmfiles/1st_passHGR_2d_parms_new_bounds.txt
pyscript=~/monarch/github/F-H_2018/scripts/dadi/run_3d_mods.py



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
