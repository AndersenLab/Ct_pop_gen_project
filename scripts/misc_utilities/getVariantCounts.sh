#!/bin/bash

source activate bcftools

#usage: query alt sites for each sample and count variants across each genomic bin
#while IFS= read -r strain; do sbatch --export=strain=$strain ../../scripts/getVariantCounts.sh; done < ../sample_list/CT_all_samples.txt
#the species hard-filtered vcf is soft-linked to processed_data/variant_counts/

bcftools view -s $strain WI.20250627.hard-filter.isotype.vcf.gz |\
bcftools filter -i 'GT="alt"' -Oz -o $strain.vcf.gz
bedtools coverage -a ../genomic_bins/ONT_NIC58_1kb_bins.bed -b $strain.vcf.gz -counts > $strain.variant_counts.tsv
rm $strain.vcf.gz
