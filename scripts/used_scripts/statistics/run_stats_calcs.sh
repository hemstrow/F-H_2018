#!/bin/bash


# all sites
echo "#!/bin/bash -l
#SBATCH -t 10-12:00:00
#SBATCH --mem=60G
#SBATCH -J par_stats_all
#SBATCH -p high
Rscript calc_stats.R allsites_paralog_fix_snpR.RDS FALSE FALSE" > all_sites.sh
sbatch all_sites.sh
rm all_sites.sh

# no maf
echo "#!/bin/bash -l
#SBATCH -t 10-12:00:00
#SBATCH --mem=40G
#SBATCH -J par_stats_nomaf
#SBATCH -p high
Rscript calc_stats.R nomaf_paralog_fix_snpR.RDS FALSE FALSE" > nomaf.sh
sbatch nomaf.sh
rm nomaf.sh

# maf
echo "#!/bin/bash -l
#SBATCH -t 10-12:00:00
#SBATCH --mem=20G
#SBATCH -J par_stats_maf
#SBATCH -p high
Rscript calc_stats.R maf_paralog_fix_snpR.RDS TRUE FALSE" > maf.sh
sbatch maf.sh
rm maf.sh
