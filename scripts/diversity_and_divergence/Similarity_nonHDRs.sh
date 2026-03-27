#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --output=Similarity_nonHDRs.oe
#SBATCH --job-name="SimiNonHDRS"


source activate CT_PopGen

cd ../../processed_data/
mkdir -p gt_check_no_HDRs
cd gt_check_no_HDRs

out_dir="../../processed_data/gt_check_no_HDRs"
IVCF="../../processed_data/pi_theta_exclude_HDRs/WI.20250627.hard-filter.isotype.noHDR.vcf.gz"

bcftools gtcheck --no-HWE-prob -e 0 $IVCF > gtcheck_noHDR.tsv
