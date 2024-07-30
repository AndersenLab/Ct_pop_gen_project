#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=LD_all_singletons_eiganstrat_input.oe
#SBATCH --job-name="vcf_to_matrix_LD6_9"



python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/vcf2phylip_master/vcf2phylip.py -i /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/LD_pruned/LD_0.6_singletons_eiganstrat_input.vcf.gz
python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/vcf2phylip_master/vcf2phylip.py -i /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/LD_pruned/LD_0.7_singletons_eiganstrat_input.vcf.gz
python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/vcf2phylip_master/vcf2phylip.py -i /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/LD_pruned/LD_0.8_singletons_eiganstrat_input.vcf.gz
python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/vcf2phylip_master/vcf2phylip.py -i /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/LD_pruned/LD_0.9_singletons_eiganstrat_input.vcf.gz

