#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=Ct_heterozygous_SNVs.oe
#SBATCH --job-name="CtHets"

source activate CT_PopGen

cd ./processed_data/
mkdir -p HDR_stats
cd HDR_stats

HDR_TABLE="../../tables/TableS7_HDR_CT_allStrain_5kbclust_20251201.tsv"
SOFT_VCF="../../data/VCF/WI.20250627.soft-filter.isotype.vcf.gz"

PREFIX="Ct_soft_filtered_hets"

# 1. Table S6 ranges -> BED
awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3}' ${HDR_TABLE} \
  | sort -k1,1 -k2,2n \
  > ${PREFIX}.TableS7_HDR_ranges.bed

# 2. Count het genotype calls per biallelic SNV site
bcftools view -m2 -M2 -v snps ${SOFT_VCF} \
  | bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
  | awk 'BEGIN{
      OFS="\t";
      print "CHROM","POS","n_het","n_called","n_total","has_het"
    }
    {
      n_het=0;
      n_called=0;
      n_total=NF-2;

      for (i=3; i<=NF; i++) {
        gt=$i;

        if (gt != "./." && gt != ".|." && gt != ".") {
          n_called++;
        }

        if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
          n_het++;
        }
      }

      print $1,$2,n_het,n_called,n_total,(n_het > 0 ? 1 : 0);
    }' \
  > ${PREFIX}.het_counts_per_site.tsv

# 3. SNV sites -> BED
awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2-1, $2, $1":"$2}' \
  ${PREFIX}.het_counts_per_site.tsv \
  > ${PREFIX}.SNV_sites.bed

# 4. Get site IDs that are in HDR ranges
### -wa -u: each SNV site is reported once even if it overlaps multiple HDR ranges.
bedtools intersect \
  -a ${PREFIX}.SNV_sites.bed \
  -b ${PREFIX}.TableS7_HDR_ranges.bed \
  -wa -u \
  | awk '{print $4}' \
  | sort -u \
  > ${PREFIX}.HDR_site_ids.txt

# 5. Add HDR/non-HDR label
awk 'BEGIN{OFS="\t"}
  NR==FNR {
    hdr[$1]=1;
    next;
  }
  FNR==1 {
    print $0,"region";
    next;
  }
  {
    site_id=$1":"$2;
    region=(site_id in hdr ? "HDR" : "non-HDR");
    print $0,region;
  }' \
  ${PREFIX}.HDR_site_ids.txt \
  ${PREFIX}.het_counts_per_site.tsv \
  > ${PREFIX}.het_counts_per_site.HDR_labeled.tsv

