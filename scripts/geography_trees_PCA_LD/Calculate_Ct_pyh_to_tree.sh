#!/bin/bash
#SBATCH --job-name="CtPyh2Tree"
#SBATCH --output=Calculate_Ct_pyh_to_tree.oe
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 48


cd ../../processed_data/LD_pruned_trees

# define input files 
input_file="../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy"

# run iqtree 
source activate tree
iqtree -s $input_file -mset GTR -bb 1000 -T 48

