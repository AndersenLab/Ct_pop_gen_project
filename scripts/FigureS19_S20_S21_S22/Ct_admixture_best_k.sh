#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Ct_admixture_best_k.oe
#SBATCH --job-name="CtAdmBeK"




cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Ct_admixture_best_k
cd Ct_admixture_best_k

ln -s $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_admixture/*.Q ./
ln -s $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/scripts/admixture_best_k/concat_Qs.sh ./
bash concat_Qs.sh





ln -s $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_admixture/*.out ./
ln -s $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/scripts/admixture_best_k/generate_CV_matrix.sh ./
bash generate_CV_matrix.sh




