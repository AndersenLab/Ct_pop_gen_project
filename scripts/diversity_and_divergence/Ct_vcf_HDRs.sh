#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=Ct_vcf_HDRs.oe
#SBATCH --job-name="VcfHDRs"


source activate CT_PopGen

cd ../../processed_data/
mkdir -p pi_theta_only_HDRs
cd pi_theta_only_HDRs

out_dir="../../processed_data/pi_theta_only_HDRs"
raw_VCF="../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"

bcftools view \
  -T ../pi_theta_exclude_HDRs/HDRs_file.tsv \
  -Oz \
  -o WI.20250627.hard-filter.isotype.onlyHDR.vcf.gz \
  $raw_VCF

bcftools index --threads 10 WI.20250627.hard-filter.isotype.onlyHDR.vcf.gz
