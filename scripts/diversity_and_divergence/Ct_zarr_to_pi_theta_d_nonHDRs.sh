#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=Ct_zarr_to_pi_theta_d_nonHDRs.oe
#SBATCH --job-name="CtZ2PtdNonHDRs"


cd ../../processed_data/
mkdir -p pi_theta_exclude_HDRs
cd pi_theta_exclude_HDRs

out_dir="../../../processed_data/pi_theta_exclude_HDRs"
chrom_length="../../../data/genome_feature/06.04.22_ct_chrom_dict"
rep_bool="../../../processed_data/make_Ct_repeats_bed_file/Ct_repeat_region.pkl"
zarr="../../../processed_data/pi_theta_exclude_HDRs/zarr/WI.20250627.hard-filter.isotype.noHDR.vcf.gz.zarr"
script_path="../../../scripts/diversity_and_divergence/pi_theta_d_python"

source activate vcf_stats

python $script_path/pi_theta_d.py \
--zarr $zarr \
--chrom_lengths $chrom_length \
--mask_missing T \
--mask_repeats T \
--repeats_file $rep_bool \
--out $out_dir \
--filter_monomorphic T 
