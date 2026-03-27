#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=geo_vcf.oe
#SBATCH --job-name="geo_vcf"


source activate CT_PopGen

geo_info="../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv"
VCFI="../../data/VCF/WI.20250331.hard-filter.isotype.vcf.gz"

cd ../../processed_data/
mkdir -p pi_theta_d_geo/vcf_and_zarr
cd pi_theta_d_geo/vcf_and_zarr

categories=("Taiwan" "Hawaii" "Caribbean" "Africa" "South America" "Central America" "Australia" "Indonesia" "Micronesia")
for category in "${categories[@]}"; do
    bcftools view -S <(awk -v category="$category" -F ',' '{if($5==category) print $2}' "$geo_info") -o "${category}.vcf" "$VCFI"
done
