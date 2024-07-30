#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=Cb_zarr_to_pi_theta_d.oe
#SBATCH --job-name="CbZPtd"




# Need to change the details in the script /caeno_diversity_stats-main/scripts/scripts/pi_theta_d.py
## I changed line 90 (the path to the .pkl file) accordingly:
## divergent_bool = pd.read_pickle("/home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/inputs/05.10.22_ce_div_bools.pkl")
## remember to be carful with indent 




cd /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_diversity



# copy the python file into bwang97 dir 
# cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/scripts/pi_theta_d.py /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts
# cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/src/vcf_stats/vcf_stat.py /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts




# copy the Cb chrom_length file into bwang97 dir 
cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/inputs/05.07.22_cb_chrom_dict /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_diversity

# # copy the .pkl file into bwang97 dir
# cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/inputs/05.10.22_ce_div_bools.pkl /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/pi_theta_d/geo_pi_theta_d

# # copy the bed file into bwang97 dir
# cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/inputs/05.01.22_divergent_regions_all.bed /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/pi_theta_d/geo_pi_theta_d


#Define the inputs and outputs
out_dir="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_diversity"
chrom_length="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_diversity/05.07.22_cb_chrom_dict"
zarr="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_diversity/zarr/WI.20240726.hard-filter.vcf.gz.zarr"
script_path="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts"

# activate conda env
source activate vcf_stats


### raw  script ###
# python scripts/scripts/pi_theta_d.py --zarr $zarr --chrom_lengths $chrom_length --mask_missing T --out $out_dir

### run the  script ###
python $script_path/pi_theta_d.py --zarr $zarr --chrom_lengths $chrom_length --mask_missing T --out $out_dir


