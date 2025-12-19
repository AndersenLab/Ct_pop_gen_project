#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Ct_zarr_to_pi_theta_d.oe
#SBATCH --job-name="CtZ2PTD"






cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/
mkdir -p pi_theta_d
cd pi_theta_d



#Define the inputs and outputs
out_dir="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/pi_theta_d"
chrom_length="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/genome_feature/06.04.22_ct_chrom_dict"
# div_bool=""
rep_bool="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/make_Ct_repeats_bed_file/Ct_repeat_region.pkl"
##### using the previous zarr file
zarr="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/pi_theta_d/zarr/WI.20250627.hard-filter.isotype.vcf.gz.zarr"
script_path="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/scripts/pi_theta_d"

# activate conda env
source activate vcf_stats


### raw script ###
##### python scripts/scripts/pi_theta_d.py --zarr $zarr --chrom_lengths $chrom_length --mask_missing T --out $out_dir


### run the script ###
python $script_path/pi_theta_d.py \
--zarr $zarr \
--chrom_lengths $chrom_length \
--mask_missing T \
--mask_repeats T \
--repeats_file $rep_bool \
--out $out_dir \
--filter_monomorphic T 



# --maf_threshold 0.05 \ default is none

# --mask_divergent T \
# --divergent_file $div_bool \


#### the default missing_threshold is 0.8, 
#### This means that if more than 80% of the strains are missing genotype information for any site, 
#### it is excluded from the calculation (considered neither variant or invariant).

