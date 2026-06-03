#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --output=Ct_1to1_4fold_vcf_to_zarr.oe
#SBATCH --job-name="Ct4fV2Zarr"

source activate CT_PopGen

cd ./processed_data/HDR_stats
mkdir -p pi_theta_d_1to1_4fold
cd pi_theta_d_1to1_4fold

out_dir="../../../processed_data/HDR_stats/pi_theta_d_1to1_4fold"
raw_VCF="../../../processed_data/HDR_stats/Ct_1to1_ortholog_4fold.hard_filter.variable_sites.vcf.gz"

mkdir -p $out_dir/vcf
mkdir -p $out_dir/zarr

vcf_name=$(basename $raw_VCF)

ln -s $raw_VCF $out_dir/vcf/$vcf_name
bcftools index $out_dir/vcf/$vcf_name
vcf_input="$out_dir/vcf/$vcf_name"

source activate vcf_zarr
zarr=$out_dir/zarr/$vcf_name.zarr
echo How many CPUs you asked for ${SLURM_NPROCS}
python ../../../scripts/pi_theta_d_python/parallel_testing.py \
--nproc ${SLURM_NPROCS} \
--vcf $vcf_input \
--out $zarr

