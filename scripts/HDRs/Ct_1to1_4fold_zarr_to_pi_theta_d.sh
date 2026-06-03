#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Ct_1to1_4fold_zarr_to_pi_theta_d.oe
#SBATCH --job-name="Ct4fZ2PTD"

cd ./processed_data/HDR_stats
mkdir -p pi_theta_d_1to1_4fold
cd pi_theta_d_1to1_4fold

out_dir="../../../processed_data/HDR_stats/pi_theta_d_1to1_4fold"
chrom_length="../../../data/genome_feature/06.04.22_ct_chrom_dict"
rep_bool="../../../processed_data/make_Ct_repeats_bed_file/Ct_repeat_region.pkl"

zarr="../../../processed_data/HDR_stats/pi_theta_d_1to1_4fold/zarr/Ct_1to1_ortholog_4fold.hard_filter.variable_sites.vcf.gz.zarr"

script_path="../../../scripts/pi_theta_d_python/parallel_testing.py"

source activate vcf_stats

python $script_path/pi_theta_d.py \
--zarr $zarr \
--chrom_lengths $chrom_length \
--mask_missing T \
--mask_repeats T \
--repeats_file $rep_bool \
--out $out_dir \
--filter_monomorphic T

