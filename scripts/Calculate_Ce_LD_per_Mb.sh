#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=Ce_LD_per_Mb.oe
#SBATCH --job-name="CeLD"



source activate CT_PopGen


VCFI_Ce="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/Ce_VCF/WI.20231213.hard-filter.isotype.vcf.gz"

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data

mkdir -p Ce_LD_all_isotype
cd Ce_LD_all_isotype



plink --vcf $VCFI_Ce \
--allow-no-sex --allow-extra-chr \
--maf 0.05 --geno 0.2 --r2 --ld-window 999999 \
--ld-window-r2 0 --out Ce_all





time cat Ce_all.ld | awk '{ wind=$1 " " int($2/100000) " " $4 " " int($5/100000); used[wind] += $7; count[wind]++; } END { for (d in count) { print d, used[d]/count[d] } }' > Ce_all.summary



