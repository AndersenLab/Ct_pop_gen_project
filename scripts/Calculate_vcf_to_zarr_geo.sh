#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --output=vcf_to_zarr_geo.oe
#SBATCH --job-name="vzg"

source activate CT_PopGen

#Define the inputs and outputs
## use full paths
out_dir=/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
vcf_raw_dir="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/LD_geo"

# copy the python file into bwang97 dir 
cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/src/vcf_to_zarr/parallel_testing.py /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts

### Copy and Index the VCF ###
#make the output directories
mkdir -p $out_dir/pi_theta_d
mkdir -p $out_dir/pi_theta_d/geo_vcf_and_zarr

# get the file name of the VCF
for vcf_file in "$vcf_raw_dir"/*.vcf; do
    vcf_name=$(basename "$vcf_file")
    cp "$vcf_file" "$out_dir/pi_theta_d/geo_vcf_and_zarr/$vcf_name"
done



# compress the files
for file in "$out_dir"/pi_theta_d/geo_vcf_and_zarr/*.vcf; do
    bgzip -f "$file"
done



# index the VCF
for vcf_file in "$out_dir"/pi_theta_d/geo_vcf_and_zarr/*.vcf.gz; do
    bcftools index -f "$vcf_file"
done




source activate vcf_zarr


# define the zarr output path relative to the conversion directory
# all vcf files
vcf_files=($out_dir/pi_theta_d/geo_vcf_and_zarr/*.vcf.gz)

# for loop
for vcf_file in "${vcf_files[@]}"; do
    # extract file names
    vcf_name=$(basename "$vcf_file")

    # zarr file path
    zarr="$out_dir/pi_theta_d/geo_vcf_and_zarr/${vcf_name%.vcf.gz}.zarr"

    # export following message
    echo "How many CPUs you asked for ${SLURM_NPROCS}"
    echo "Processing file: $vcf_file"
    echo "Output Zarr file: $zarr"

    # run python script
    python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/parallel_testing.py \
        --nproc "${SLURM_NPROCS}" \
        --vcf "$vcf_file" \
        --out "$zarr"

    echo "Finished processing $vcf_file"
done
