#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=Ct_vcf_nonHDRs_geo.oe
#SBATCH --job-name="nonHDRGeo"


source activate CT_PopGen
geo_info="../../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv"
VCFI="../../../processed_data/pi_theta_exclude_HDRs/WI.20250627.hard-filter.isotype.noHDR.vcf.gz"

cd ../../processed_data/
mkdir -p geo_vcf_noHDR
cd geo_vcf_noHDR

categories=("Taiwan" "Hawaii" "Caribbean" "Africa" "South America" "Central America" "Australia" "Indonesia" "Micronesia")

for category in "${categories[@]}"; do
    bcftools view \
      -S <(awk -v category="$category" -F ',' 'NR>1 {if($5==category) print $2}' "$geo_info") \
      -Oz \
      -o "${category}.noHDR.vcf.gz" \
      "$VCFI"

    bcftools index -f "${category}.noHDR.vcf.gz"
done
