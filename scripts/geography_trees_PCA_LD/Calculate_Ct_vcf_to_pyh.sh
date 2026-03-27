#!/bin/bash
#SBATCH --job-name="CtVCF2Pyh"
#SBATCH --output=Calculate_Ct_vcf_to_pyh.out
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 6


cd ../../processed_data/
mkdir -p LD_pruned_trees
cd LD_pruned_trees

# define input files 
input_file="../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/INPUTFILES/eiganstrat_input.vcf.gz"

mkdir -p LD_0.9
cd LD_0.9

python ../../scripts/vcf2phylip_master/vcf2phylip.py -i $input_file

mv eiganstrat_input.min4.phy phy_file_0.9.phy
