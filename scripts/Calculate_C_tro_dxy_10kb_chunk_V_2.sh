#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=C_tro_dxy_chunk_V_2.oe
#SBATCH --job-name="dxycV2"


source activate pixy

### calculate genetic distance (Dxy) using pixy
### pixy is a command-line tool for painlessly and correctly estimating average nucleotide diversity within (π) and between (dxy) populations
### https://github.com/ksamuk/pixy
### https://pixy.readthedocs.io/en/latest/examples.html#basic-usage


### please use python version 3.8.1 for installation when creating conda environment 
### conda create -n "pixy" python=3.8.1
### conda install -c bioconda htslib
### conda install -c bioconda samtools=1.9 --force-reinstall -y
 



cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
mkdir -p dxy


VCFI="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz"
VCF_Index="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz.tbi"
geo_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/geo_populations_file.txt"
all_pairs_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/all_pairs_pop_file.txt"
Out_dir="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/dxy"





# Run Pixy for the specified chromosome
pixy --stats dxy \
--vcf $VCFI \
--populations $all_pairs_pop_file \
--output_folder $Out_dir \
--window_size 10000 \
--bypass_invariant_check yes \
--output_prefix chunk_all_pairs_pop_10kb_chr_V_2 \
--n_cores 48 \
--chromosomes V \
--interval_start 7660001 \
--interval_end 15340000\



