#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --output=extract1to1_4fold.oe
#SBATCH --job-name="Ct4foldExtract"

source activate CT_PopGen

cd ./processed_data/HDR_stats

ORTHO="N2_NIC58_1to1_orthologs.tsv"
DEGEN="degenotate_NIC58/degeneracy-all-sites.bed"
HARD_VCF="../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"

PREFIX="Ct_1to1_ortholog_4fold"

echo "ORTHO: ${ORTHO}"
echo "DEGEN: ${DEGEN}"
echo "HARD_VCF: ${HARD_VCF}"

ls -lh ${ORTHO}
ls -lh ${DEGEN}
ls -lh ${HARD_VCF}
ls -lh ${HARD_VCF}.tbi

# 1. Extract Ct 1:1 ortholog transcript IDs
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $3}' ${ORTHO} \
  | sed 's/^transcript://' \
  | sed 's/^transcript_//' \
  | sed 's/^protein://' \
  | sed 's/^protein_//' \
  | sed 's/[[:space:]].*$//' \
  | sort -u \
  > ${PREFIX}.Ct_1to1_transcript_ids.clean.txt

echo "Number of Ct 1:1 transcripts:"
wc -l ${PREFIX}.Ct_1to1_transcript_ids.clean.txt
head ${PREFIX}.Ct_1to1_transcript_ids.clean.txt

# 2. Extract all 4-fold degenerate sites from degenotate output
awk 'BEGIN{FS=OFS="\t"} $5 == 4 {
  split($4,a,":");
  tx=a[2];
  print $1,$2,$3,tx,$4,$5,$6,$7
}' ${DEGEN} \
  | sort -k1,1 -k2,2n \
  > ${PREFIX}.all_4fold_sites.with_tx.bed

echo "All 4-fold sites:"
wc -l ${PREFIX}.all_4fold_sites.with_tx.bed
head ${PREFIX}.all_4fold_sites.with_tx.bed


# 3. Keep only 4-fold sites in Ct 1:1 ortholog transcripts
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
  keep[$1]=1;
  next;
}
{
  tx=$4;
  if (tx in keep) {
    print $0;
  }
}' \
  ${PREFIX}.Ct_1to1_transcript_ids.clean.txt \
  ${PREFIX}.all_4fold_sites.with_tx.bed \
  > ${PREFIX}.1to1_4fold_sites.with_tx.bed

echo "1:1 ortholog 4-fold sites:"
wc -l ${PREFIX}.1to1_4fold_sites.with_tx.bed
head ${PREFIX}.1to1_4fold_sites.with_tx.bed

# Unique genomic BED for bcftools
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' \
  ${PREFIX}.1to1_4fold_sites.with_tx.bed \
  | sort -k1,1 -k2,2n \
  | uniq \
  > ${PREFIX}.1to1_4fold_sites.3col.bed

echo "Unique genomic 1:1 ortholog 4-fold sites:"
wc -l ${PREFIX}.1to1_4fold_sites.3col.bed
head ${PREFIX}.1to1_4fold_sites.3col.bed

# 4. Extract variable SNVs at these 1:1 ortholog 4-fold sites
bcftools view \
  -R ${PREFIX}.1to1_4fold_sites.3col.bed \
  -m2 -M2 -v snps \
  -Oz \
  -o ${PREFIX}.hard_filter.variable_sites.vcf.gz \
  ${HARD_VCF}

tabix -p vcf ${PREFIX}.hard_filter.variable_sites.vcf.gz

echo "Variable SNVs at 1:1 ortholog 4-fold sites:"
bcftools view -H ${PREFIX}.hard_filter.variable_sites.vcf.gz | wc -l

echo "Output files for downstream pi:"
ls -lh ${PREFIX}.1to1_4fold_sites.with_tx.bed
ls -lh ${PREFIX}.1to1_4fold_sites.3col.bed
ls -lh ${PREFIX}.hard_filter.variable_sites.vcf.gz
ls -lh ${PREFIX}.hard_filter.variable_sites.vcf.gz.tbi

echo "Done."

