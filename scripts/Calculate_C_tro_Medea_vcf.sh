#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=C_tro_Medea_vcf_20240508.oe
#SBATCH --job-name="CtMedVcf58"

# source activate vcf_zarr
source activate CT_PopGen

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
mkdir -p Medea
cd Medea

VCFI="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz"


vcftools --gzvcf $VCFI --chr III --from-bp 1270000 --to-bp 1510000 --recode --out Medea_Chr3_127_151
vcftools --gzvcf $VCFI --chr V --from-bp 12940000 --to-bp 13320000 --recode --out Medea_Chr5_1294_1332


