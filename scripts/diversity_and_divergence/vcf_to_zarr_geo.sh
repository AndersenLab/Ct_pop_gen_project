#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=vcf_to_zarr_geo.oe
#SBATCH --job-name="V2ZaGeo"


source activate CT_PopGen

out_dir="../../processed_data/"
script_path="./pi_theta_d_python"

mkdir -p $out_dir/pi_theta_d_geo
mkdir -p $out_dir/pi_theta_d_geo/vcf_and_zarr

for file in "$out_dir"/pi_theta_d_geo/vcf_and_zarr/*.vcf; do
    bgzip -f "$file"
done

for vcf_file in "$out_dir"/pi_theta_d_geo/vcf_and_zarr/*.vcf.gz; do
    bcftools index -f "$vcf_file"
done

source activate vcf_zarr
vcf_files=($out_dir/pi_theta_d_geo/vcf_and_zarr/*.vcf.gz)
for vcf_file in "${vcf_files[@]}"; do
    vcf_name=$(basename "$vcf_file")
    zarr="$out_dir/pi_theta_d_geo/vcf_and_zarr/${vcf_name%.vcf.gz}.zarr"
    echo "How many CPUs you asked for ${SLURM_NPROCS}"
    echo "Processing file: $vcf_file"
    echo "Output Zarr file: $zarr"
    python $script_path/parallel_testing.py \
        --nproc "${SLURM_NPROCS}" \
        --vcf "$vcf_file" \
        --out "$zarr"
    echo "Finished processing $vcf_file"
done
