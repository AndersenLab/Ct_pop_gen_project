#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=zarr_to_pi_theta_d_geo.oe
#SBATCH --job-name="Ct_z_ptd"


cd ../../processed_data
mkdir -p pi_theta_d_geo
cd pi_theta_d_geo

out_dir="../../processed_data/pi_theta_d_geo"
chrom_length="../../data/genome_feature/06.04.22_ct_chrom_dict"
rep_bool="../../processed_data/make_Ct_repeats_bed_file/Ct_repeat_region.pkl"
zarr="../../processed_data/pi_theta_d_geo/vcf_and_zarr"
script_path="./pi_theta_d_python"

source activate vcf_stats

cd $zarr
mv Central\ America.zarr Central_America.zarr
mv South\ America.zarr South_America.zarr

zarr_dirs=$(find . -type d -name "*.zarr")

for zarr_dir_name in $zarr_dirs; do
    dir_name=$(basename "$zarr_dir_name"| sed 's/\.zarr$//')
    out_dir_geo="$out_dir/$dir_name"
    mkdir -p "$out_dir_geo"
    python "$script_path/pi_theta_d.py" --zarr "$zarr_dir_name" --chrom_lengths "$chrom_length" --mask_missing T --mask_repeats T --repeats_file $rep_bool --filter_monomorphic T --out "$out_dir_geo"
done
