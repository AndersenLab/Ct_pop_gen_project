#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=LD_per_Mb_geo.oe
#SBATCH --job-name="LD_per_Mb_geo"


source activate CT_PopGen
# I have installed bcftools and plink in conda

geo_info="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/indep_isotype_info_geo.csv"
# VCFI="/home/bwang97/vast-eande106/data/c_tropicalis/WI/variation/20231201/vcf/WI.20231201.hard-filter.isotype.vcf.gz"

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
mkdir -p LD_geo
cd LD_geo


# Loop through each category in the geo_info file
categories=("Taiwan" "Hawaii" "Caribbean" "Africa" "South America" "Central America" "Australia")

for category in "${categories[@]}"; do
    # Use bcftools view with awk to filter species based on the category
    bcftools view -S <(awk -v category="$category" -F ',' '{if($4==category) print $1}' "$geo_info") -o "${category}.vcf" "$VCFI"
done




for vcf_file in *.vcf; do

plink --vcf "$vcf_file" \
        --allow-no-sex --allow-extra-chr \
        --maf 0.05 --geno 0.2 --r2 --ld-window 999999 \
        --ld-window-r2 0 --out "${vcf_file%.vcf}_LD"
done



for LD_file in *.ld; do
    cat "$LD_file" | awk '{ wind=$1 " " int($2/100000) " " $4 " " int($5/100000); used[wind] += $7; count[wind]++; } END { for (d in used) { print d, used[d]/count[d] } }' > "C_tro_${LD_file%.LD}_LD.summary"
done
