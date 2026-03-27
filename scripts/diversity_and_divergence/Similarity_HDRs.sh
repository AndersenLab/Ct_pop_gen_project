#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --output=Similarity_HDRs.oe
#SBATCH --job-name="SimiHDRs"


source activate CT_PopGen

cd ../../processed_data/
mkdir -p gt_check_only_HDRs
cd gt_check_only_HDRs

out_dir="../../processed_data/gt_check_only_HDRs"
IVCF="../../processed_data/pi_theta_only_HDRs/WI.20250627.hard-filter.isotype.onlyHDR.vcf.gz"

bcftools gtcheck --no-HWE-prob -e 0 $IVCF > gtcheck_onlyHDR.tsv
