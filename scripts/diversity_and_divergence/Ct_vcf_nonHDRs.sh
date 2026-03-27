#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=Ct_vcf_nonHDRs.oe
#SBATCH --job-name="VcfNonHDRs"


source activate CT_PopGen

cd ../../processed_data/
mkdir -p pi_theta_exclude_HDRs
cd pi_theta_exclude_HDRs

out_dir="../../processed_data/pi_theta_exclude_HDRs"
raw_VCF="../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"

bcftools view \
  -T ^HDRs_file.tsv \
  -Oz \
  -o WI.20250627.hard-filter.isotype.noHDR.vcf.gz \
  $raw_VCF

bcftools index WI.20250627.hard-filter.isotype.noHDR.vcf.gz
