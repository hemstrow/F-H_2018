#!/bin/bash


# all sites
echo "#!/bin/bash -l
#SBATCH -t 10-12:00:00
#SBATCH --mem=60G
#SBATCH -J par_stats_all
module load R/4.1.0
Rscript calc_stats.R ../../../results/paralogs/allsites_paralog_fix_snpR.RDS" > all_sites.sh
sbatch all_sites.sh
rm all_sites.sh
