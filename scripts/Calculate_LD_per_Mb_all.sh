#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=LD_per_Mb_all.oe
#SBATCH --job-name="LD_per_Mb_all"



source activate CT_PopGen

# VCFI="/home/bwang97/vast-eande106/data/c_tropicalis/WI/variation/20231201/vcf/WI.20231201.hard-filter.isotype.vcf.gz"

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data

mkdir -p LD_all_isotype
cd LD_all_isotype



plink --vcf $VCFI \
--allow-no-sex --allow-extra-chr \
--maf 0.05 --geno 0.2 --r2 --ld-window 999999 \
--ld-window-r2 0 --out C_tro_all





time cat C_tro_all.ld | awk '{ wind=$1 " " int($2/100000) " " $4 " " int($5/100000); used[wind] += $7; count[wind]++; } END { for (d in count) { print d, used[d]/count[d] } }' > C_tro_all.summary



