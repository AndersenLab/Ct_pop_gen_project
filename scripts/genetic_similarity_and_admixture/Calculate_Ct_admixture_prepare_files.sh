#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_Ct_admixture_prepare_files.oe
#SBATCH --job-name="CtAdmPre"

cd ../../processed_data
mkdir -p Ct_admixture
cd Ct_admixture

out_folder="../../processed_data/Ct_admixture/"
ln -s "../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/INPUTFILES/eigenstrat_input.ped" "${out_folder}LD_0.9.ped"

### generate seeds
cd ../../processed_data
mkdir -p Ct_admixture
cd Ct_admixture
shuf -i 1000-2000 -n 10 > seeds.txt
