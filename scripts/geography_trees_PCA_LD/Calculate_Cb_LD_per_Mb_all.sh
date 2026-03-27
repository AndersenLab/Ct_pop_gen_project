#!/bin/bash
#SBATCH -A eande106_bigmem
#SBATCH -p bigmem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mem-per-cpu=25G
#SBATCH --output=Calculate_Cb_LD_per_Mb_all.oe
#SBATCH --job-name="CbLdAll"


source activate CT_PopGen
VCFI="../../data/VCF/Cb/WI.20250626.hard_filter.715_isotype.vcf.gz"

cd ../../processed_data/
mkdir -p LD_three_species
cd LD_three_species

plink --vcf $VCFI \
--allow-no-sex --allow-extra-chr \
--threads 15 \
--maf 0.05 --geno 0.2 --r2 --ld-window 999999 \
--ld-window-r2 0 --out C_bri_all

time cat C_bri_all.ld | awk '{ wind=$1 " " int($2/100000) " " $4 " " int($5/100000); used[wind] += $7; count[wind]++; } END { for (d in count) { print d, used[d]/count[d] } }' > C_bri_all.summary

rm C_bri_all.ld
