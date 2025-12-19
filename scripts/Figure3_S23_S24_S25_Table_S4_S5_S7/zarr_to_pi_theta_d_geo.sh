#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=zarr_to_pi_theta_d_geo.oe
#SBATCH --job-name="Ct_z_ptd"




# make new dir
cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
mkdir -p pi_theta_d_geo
cd pi_theta_d_geo


#Define the inputs and outputs
out_dir="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/pi_theta_d_geo"
chrom_length="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/genome_feature/06.04.22_ct_chrom_dict"
rep_bool="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/make_Ct_repeats_bed_file/Ct_repeat_region.pkl"
zarr="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/pi_theta_d_geo/vcf_and_zarr"
script_path="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/scripts/pi_theta_d_by_geo"

# activate conda env
source activate vcf_stats


### raw script ###
# python scripts/scripts/pi_theta_d.py --zarr $zarr --chrom_lengths $chrom_length --mask_missing T --out $out_dir


## make it a for-loop

# find all dir end with .zarr
cd $zarr
mv Central\ America.zarr Central_America.zarr
mv South\ America.zarr South_America.zarr
mv Malay\ Archipelago.zarr Malay_Archipelago.zarr


zarr_dirs=$(find . -type d -name "*.zarr")

# loop every dir
for zarr_dir_name in $zarr_dirs; do
    # extract dir name (.zarr)
    dir_name=$(basename "$zarr_dir_name"| sed 's/\.zarr$//')
    # define export dir
    out_dir_geo="$out_dir/$dir_name"
    # Create output directory if it doesn't exist
    mkdir -p "$out_dir_geo"
    # run the script
    python "$script_path/pi_theta_d.py" --zarr "$zarr_dir_name" --chrom_lengths "$chrom_length" --mask_missing T --mask_repeats T --repeats_file $rep_bool --filter_monomorphic T --out "$out_dir_geo"
done


