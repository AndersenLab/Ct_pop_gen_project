#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 60:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=vcf_to_zarr_HVRs.oe
#SBATCH --job-name="v2z_HVRs"




# source activate vcf_zarr

# source activate CT_PopGen

# module load bcftools/1.15.1


#Define the inputs and outputs
## use full paths

vcf_file="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/VCF/WI.20231201.hard-filter.isotype.vcf.gz"

#make the output directories
mkdir -p /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/pi_theta_d/HVRs_vcf_and_zarr
out_dir=/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/pi_theta_d/HVRs_vcf_and_zarr


# copy the python file into bwang97 dir 
# cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/src/vcf_to_zarr/parallel_testing.py /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts





source activate vcf_zarr

## use this page to install scikit-allel https://github.com/cggh/scikit-allel/issues/203
## conda install -c conda-forge scikit-allel
# (DON'T! ###not necessary###)module load python/3.8.6


# zarr file path
    # zarr="$out_dir/pi_theta_d/geo_vcf_and_zarr/${vcf_name%.vcf.gz}.zarr"

    # export following message
    echo "How many CPUs you asked for ${SLURM_NPROCS}"
    echo "Processing file: $vcf_file"
    echo "Output Zarr file: $out_dir"

    # run python script
    python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/parallel_testing.py \
        --nproc "${SLURM_NPROCS}" \
        --vcf "$vcf_file" \
        --out "$out_dir"

    echo "Finished processing $vcf_file"

