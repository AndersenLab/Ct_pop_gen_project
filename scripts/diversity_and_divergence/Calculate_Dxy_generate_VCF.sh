#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --array=0-7
#SBATCH --output=Calculate_Dxy_generate_VCF_%A_%a.oe
#SBATCH --job-name="CtDxyVcf"


source activate CT_PopGen

cd ../../processed_data/
mkdir -p Dxy_LAC
cd Dxy_LAC
mkdir -p vcf

vcf_raw="../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"

sample_lists=(LAC_*.txt)
sample_list=${sample_lists[$SLURM_ARRAY_TASK_ID]}
prefix=$(basename ${sample_list} .txt)

vcf_out="./vcf/WI.20250627.hard-filter.isotype_${prefix}.vcf.gz"

echo "Processing ${sample_list}"
echo "Output ${vcf_out}"

bcftools view -S <(cut -f1 ${sample_list}) ${vcf_raw} \
  | bcftools view -i 'COUNT(GT="0/0") > 0 && COUNT(GT="1/1") > 0' \
  -Oz -o ${vcf_out}
bcftools index -t ${vcf_out}
