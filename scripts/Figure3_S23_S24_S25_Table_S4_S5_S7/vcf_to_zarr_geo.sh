#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=vcf_to_zarr_geo.oe
#SBATCH --job-name="V2ZaGeo"

# source activate vcf_zarr
source activate CT_PopGen



#Define the inputs and outputs
## use full paths

out_dir="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/"
vcf_raw_dir="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/LD_geo"
script_path="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/scripts/pi_theta_d"


# upload the updated python file parallel_testing.py into bwang97 dir 

### Copy and Index the VCF ###
#make the output directories
mkdir -p $out_dir/pi_theta_d_geo
mkdir -p $out_dir/pi_theta_d_geo/vcf_and_zarr

# get the file name of the VCF
for vcf_file in "$vcf_raw_dir"/*.vcf; do
    vcf_name=$(basename "$vcf_file")
    cp "$vcf_file" "$out_dir/pi_theta_d_geo/vcf_and_zarr/$vcf_name"
done



# compress the files
for file in "$out_dir"/pi_theta_d_geo/vcf_and_zarr/*.vcf; do
    bgzip -f "$file"
done



# index the VCF
for vcf_file in "$out_dir"/pi_theta_d_geo/vcf_and_zarr/*.vcf.gz; do
    bcftools index -f "$vcf_file"
done




source activate vcf_zarr

## use this page to install scikit-allel https://github.com/cggh/scikit-allel/issues/203
## conda install -c conda-forge scikit-allel
# (DON'T! ###not necessary###)module load python/3.8.6

# define the zarr output path relative to the conversion directory
# all vcf files
vcf_files=($out_dir/pi_theta_d_geo/vcf_and_zarr/*.vcf.gz)

# for loop
for vcf_file in "${vcf_files[@]}"; do
    # extract file names
    vcf_name=$(basename "$vcf_file")

    # zarr file path
    zarr="$out_dir/pi_theta_d_geo/vcf_and_zarr/${vcf_name%.vcf.gz}.zarr"

    # export following message
    echo "How many CPUs you asked for ${SLURM_NPROCS}"
    echo "Processing file: $vcf_file"
    echo "Output Zarr file: $zarr"

    # run python script
    python $script_path/parallel_testing.py \
        --nproc "${SLURM_NPROCS}" \
        --vcf "$vcf_file" \
        --out "$zarr"

    echo "Finished processing $vcf_file"
done
