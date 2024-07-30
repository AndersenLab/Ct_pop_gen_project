#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --output=C_tro_fst_geo_20240511.oe
#SBATCH --job-name="fstGeo511"


source activate pixy

### calculate genetic distance (Dxy) using pixy
### pixy is a command-line tool for painlessly and correctly estimating average nucleotide diversity within (π) and between (dxy) populations
### https://github.com/ksamuk/pixy
### https://pixy.readthedocs.io/en/latest/examples.html#basic-usage


### please use python version 3.8.1 for installation when creating conda environment 
### conda create -n "pixy" python=3.8.1
### conda install -c bioconda htslib
### conda install -c bioconda samtools=1.9 --force-reinstall -y
 



cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/fst_dxy/
mkdir -p fst/geo
cd fst/geo


# cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/
# head -n 10 "/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/test_Hawaii_pop_file.txt" >test_Hawaii_10_pop_file.txt

# Hawaii_10_pairs_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/test_Hawaii_10_pop_file.txt"




file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/fst_dxy/Fst_dxy_bed_file.txt"

VCFI="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz"
geo_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/geo_populations_file.txt"
# all_pairs_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/all_pairs_pop_file.txt"
# Hawaii_pairs_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/test_Hawaii_pop_file.txt"
Out_dir="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/fst_dxy/fst/geo"



# Run Pixy for the specified chromosome
pixy --stats fst \
--vcf $VCFI \
--populations $geo_pop_file \
--output_folder $Out_dir \
--bypass_invariant_check yes \
--output_prefix C_tro_fst_geo \
--n_cores 15 \
--bed_file $file

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/fst_dxy/fst/geo
sed -i 's/Central America/Central_America/g' C_tro_fst_geo_fst.txt
sed -i 's/South America/South_America/g' C_tro_fst_geo_fst.txt


