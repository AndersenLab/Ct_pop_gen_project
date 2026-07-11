#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Ct_make_1to1_4fold_sites_file_for_pixy.oe
#SBATCH --job-name="Ct4foldSites"

source activate CT_PopGen

cd ../../processed_data/HDR_stats

ORTHO="N2_NIC58_1to1_orthologs.tsv"
DEGEN="degenotate_NIC58/degeneracy-all-sites.bed"

PREFIX="Ct_1to1_ortholog_4fold"

echo "ORTHO: ${ORTHO}"
echo "DEGEN: ${DEGEN}"

ls -lh ${ORTHO}
ls -lh ${DEGEN}

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

echo "1:1 ortholog 4-fold BED sites:"
wc -l ${PREFIX}.1to1_4fold_sites.with_tx.bed
head ${PREFIX}.1to1_4fold_sites.with_tx.bed

# Pixy --sites_file needs two columns:
# CHROM POS
#
# The degenotate BED is 0-based half-open.
# Pixy sites_file uses VCF-style 1-based POS.
# Therefore POS = BED start + 1.
awk 'BEGIN{FS=OFS="\t"} {
  print $1, $2 + 1
}' ${PREFIX}.1to1_4fold_sites.with_tx.bed \
  | sort -k1,1 -k2,2n \
  | uniq \
  > ${PREFIX}.1to1_4fold_sites.pixy_sites.txt

echo "Pixy sites_file:"
wc -l ${PREFIX}.1to1_4fold_sites.pixy_sites.txt
head ${PREFIX}.1to1_4fold_sites.pixy_sites.txt

echo "Done."
