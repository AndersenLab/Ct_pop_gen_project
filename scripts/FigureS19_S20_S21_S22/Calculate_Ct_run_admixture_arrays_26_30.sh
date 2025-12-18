#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=Calculate_Ct_run_admixture_arrays_26_30_%A_%a.oe
#SBATCH --job-name="CtAdm"
#SBATCH --array=1-50





cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Ct_admixture
cd Ct_admixture

out_folder="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_admixture/"



####### process seeds #######
seed_file="${out_folder}seeds.txt"
if [ -f "$seed_file" ]; then
    echo "Loading seeds from $seed_file..."
    readarray -t seeds < "$seed_file"
else
    echo "Error: $seed_file not found."
    exit 1
fi

####### Job Array logic #######
pops=($(seq 26 30))  # K range
total_pops=${#pops[@]}
total_seeds=${#seeds[@]}

# calculate current job's K value and seed
task_id=$SLURM_ARRAY_TASK_ID
task_idx=$((task_id - 1))  # convert as 0-based index

pop_idx=$((task_idx / total_seeds))
seed_idx=$((task_idx % total_seeds))

current_pop=${pops[$pop_idx]}
current_seed=${seeds[$seed_idx]}

####### ADMIXTURE analysis #######
echo "Processing K=${current_pop} with seed=${current_seed}"

# check if the file has already exist
if [ -f "${out_folder}LD_0.9_${current_pop}_${current_seed}.P" ] && \
   [ -f "${out_folder}LD_0.9_${current_pop}_${current_seed}.Q" ]; then
    echo "Existing results found. Skipping K=${current_pop} seed=${current_seed}"
    exit 0
fi

# run ADMIXTURE
source activate ADMIXTURE

admixture --cv=10 -s "$current_seed" \
          "${out_folder}LD_0.9.ped" \
          $current_pop \
          -j10 | tee "${out_folder}log_${current_pop}_${current_seed}.out"

# rename the output file
if [ $? -eq 0 ]; then
    mv "${out_folder}LD_0.9.${current_pop}.P" "${out_folder}LD_0.9_${current_pop}_${current_seed}.P"
    mv "${out_folder}LD_0.9.${current_pop}.Q" "${out_folder}LD_0.9_${current_pop}_${current_seed}.Q"
else
    echo "ADMIXTURE failed for K=${current_pop} seed=${current_seed}"
    exit 1
fi



#### check if all complete
#### for i in {2..25}; do   ll | grep "_${i}_" | wc -l; done

