#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=LD_0.6_singletons_eiganstrat_input.oe
#SBATCH --job-name="LD_pruned_0.6_matrix_to_tree"



module load anaconda3/2022.05
conda activate tree

iqtree -s /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/LD_pruned_trees/LD_0.6_singletons_eiganstrat_input.min4.phy -bb 1000 -T 24
