#!/bin/bash
#SBATCH -A eande106_bigmem
#SBATCH -p bigmem
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 26
#SBATCH --mem=780G
#SBATCH --output=Calculate_Cb_LD.oe
#SBATCH --job-name="Calculate_Cb_LD"



source activate CT_PopGen

VCFI="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_VCF/variation/WI.20240726.hard-filter.vcf.gz"

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data

mkdir -p Cb_diversity/LD
cd Cb_diversity/LD



plink --vcf $VCFI \
--allow-no-sex --allow-extra-chr \
--maf 0.05 --geno 0.2 --r2 --ld-window 999999 \
--ld-window-r2 0 --out Cb_all --threads 26





time cat Cb_all.ld | awk '{ wind=$1 " " int($2/100000) " " $4 " " int($5/100000); used[wind] += $7; count[wind]++; } END { for (d in count) { print d, used[d]/count[d] } }' > Cb_all.summary



