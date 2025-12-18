#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Ct_admixture_prepare_files.oe
#SBATCH --job-name="CtAdmPre"






cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Ct_admixture
cd Ct_admixture

out_folder="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_admixture/"



############ copy .ped file and rename it
cp -n "$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/INPUTFILES/eigenstrat_input.ped" "${out_folder}LD_0.9.ped"




### generate seeds
cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Ct_admixture
cd Ct_admixture

shuf -i 1000-2000 -n 10 > seeds.txt


