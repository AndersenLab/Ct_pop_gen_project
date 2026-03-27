#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Ct_admixture_best_k.oe
#SBATCH --job-name="CtAdmBeK"


cd ../../processed_data
mkdir -p Ct_admixture_best_k
cd Ct_admixture_best_k

ln -s ../../../processed_data/Ct_admixture/*.Q ./
ln -s ../../../scripts/genetic_similarity_and_admixture/concat_Qs.sh ./
bash concat_Qs.sh

ln -s ../../../processed_data/Ct_admixture/*.out ./
ln -s ../../../scripts/genetic_similarity_and_admixture/generate_CV_matrix.sh ./
bash generate_CV_matrix.sh

