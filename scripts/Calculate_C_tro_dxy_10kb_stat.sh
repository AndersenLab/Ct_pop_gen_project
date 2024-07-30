#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=avg_dxy_chunk.oe
#SBATCH --job-name=t_a_dxy_Mt



cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/dxy

# mv chunk_all_pairs_pop_10kb_chr_MtDNA_dxy.txt test_chunk_all_pairs_pop_10kb_chr_MtDNA_dxy.txt
# mv C_tro_dxy_MtDNA.summary test_C_tro_dxy_MtDNA.summary

for chunk_dxy_file in chunk_*; do
time cat "$chunk_dxy_file" | awk '{ wind=$3 " " $4 " " $5; used[wind] += $6; count[wind]++; } END { for (d in used) { print d, used[d]/count[d] } }' > "C_tro_species_wide_${chunk_dxy_file%.txt}.summary";done
cat C_tro_species_wide_*.summary > species_wide_dxy_summary.txt



for chunk_dxy_file in chunk_*; do
time cat "$chunk_dxy_file" | awk '{ wind=$1 " " $2; used[wind] += $6; count[wind]++; } END { for (d in used) { print d, used[d]/count[d] } }' > "C_tro_all_518_paired_${chunk_dxy_file%.txt}.summary";done
cat C_tro_all_518_paired_*.summary > all_518_paired_dxy_summary.txt
