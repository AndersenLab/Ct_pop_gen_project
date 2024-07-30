#!/bin/bash

# set chromosomes
chromosomes=("I" "II" "III" "IV" "X")

for chromosome in "${chromosomes[@]}"; do
    # set interval_start, interval_end, output_prefix_suffix for first case
    if [ "$chromosome" == "I" ]; then
        interval_start=6000001
        interval_end=12300000
    elif [ "$chromosome" == "II" ]; then
        interval_start=6000001
        interval_end=12740000
    elif [ "$chromosome" == "III" ]; then
        interval_start=6000001
        interval_end=11890000
    elif [ "$chromosome" == "IV" ]; then
        interval_start=6000001
        interval_end=13960000
    elif [ "$chromosome" == "X" ]; then
        interval_start=6000001
        interval_end=15980000
    else
        interval_start=0
        interval_end=0
    fi
    
    # submit job for first case
    sbatch <<EOT
#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=C_tro_dxy_chunk_${chromosome}_2.oe
#SBATCH --job-name="dxyc_${chromosome}_2"


source activate pixy

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
mkdir -p dxy

VCFI="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz"
VCF_Index="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz.tbi"
geo_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/geo_populations_file.txt"
all_pairs_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/all_pairs_pop_file.txt"
Out_dir="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/dxy"

pixy --stats dxy \
--vcf \$VCFI \
--populations \$all_pairs_pop_file \
--output_folder \$Out_dir \
--window_size 10000 \
--bypass_invariant_check yes \
--output_prefix chunk_all_pairs_pop_10kb_chr_"$chromosome"_2 \
--n_cores 48 \
--chromosomes $chromosome \
--interval_start $interval_start \
--interval_end $interval_end
EOT

    # set interval_start, interval_end, output_prefix_suffix for second case
    if [ "$chromosome" == "I" ]; then
        interval_start=1
        interval_end=6000000
    elif [ "$chromosome" == "II" ]; then
        interval_start=1
        interval_end=6000000
    elif [ "$chromosome" == "III" ]; then
        interval_start=1
        interval_end=6000000
    elif [ "$chromosome" == "IV" ]; then
        interval_start=1
        interval_end=6000000
    elif [ "$chromosome" == "X" ]; then
        interval_start=1
        interval_end=6000000
    else
        interval_start=0
        interval_end=0
    fi
    
    # submit job for second case
    sbatch <<EOT
#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=C_tro_dxy_chunk_${chromosome}_1.oe
#SBATCH --job-name="dxyc_${chromosome}_1"


source activate pixy

cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
mkdir -p dxy

VCFI="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz"
VCF_Index="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz.tbi"
geo_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/geo_populations_file.txt"
all_pairs_pop_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Geo_info/all_pairs_pop_file.txt"
Out_dir="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/dxy"

pixy --stats dxy \
--vcf \$VCFI \
--populations \$all_pairs_pop_file \
--output_folder \$Out_dir \
--window_size 10000 \
--bypass_invariant_check yes \
--output_prefix chunk_all_pairs_pop_10kb_chr_"$chromosome"_1 \
--n_cores 48 \
--chromosomes $chromosome \
--interval_start $interval_start \
--interval_end $interval_end
EOT

done
