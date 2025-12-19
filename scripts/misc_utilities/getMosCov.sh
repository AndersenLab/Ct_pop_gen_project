#!/bin/bash

#usage: loop through a sample list, retrieve coverage stats for each sample using alignment BAMs and 1kb windows file
#while IFS= read -r strain; do sbatch --export=strain=$strain ../../scripts/getMosCov.sh; done < ../sample_list/CT_all_samples.txt 

source activate bcftools

mosdepth $strain alignments/${strain}.bam -b ../genomic_bins/ONT_NIC58_1kb_bins.bed -t 4 -T 1,2,5 -n
