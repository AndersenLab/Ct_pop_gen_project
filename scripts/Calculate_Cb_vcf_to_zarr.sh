#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --output=Cb_vcf_to_zarr.oe
#SBATCH --job-name="cbvzg"


source activate CT_PopGen



# # copy the python file into bwang97 dir 
# cp /home/bwang97/vast-eande106/projects/Ryan/cb_pop_gen/scripts/src/vcf_to_zarr/parallel_testing.py /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts


#Define the inputs and outputs
## use full paths
out_dir=/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data
vcf_raw="/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/Cb_VCF/variation/WI.20240726.hard-filter.vcf.gz"


### Copy and Index the VCF ###
#make the output directories
mkdir -p $out_dir/Cb_diversity/vcf
mkdir -p $out_dir/Cb_diversity/zarr


# get the file name of the VCF
vcf_name=$(basename $vcf_raw)
cp $vcf_raw $out_dir/Cb_diversity/vcf/$vcf_name


# index the VCF
bcftools index $out_dir/Cb_diversity/vcf/$vcf_name




source activate vcf_zarr

#define the zarr output path relative to the conversion directory
zarr=$out_dir/Cb_diversity/zarr/$vcf_name.zarr



echo How many CPUs you asked for ${SLURM_NPROCS}
python /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts/parallel_testing.py \
--nproc ${SLURM_NPROCS} \
--vcf $out_dir/Cb_diversity/vcf/$vcf_name \
--out $zarr


