########################################################################
############################## PCA all ###################################
########################################################################

# This script need to be run in a cluster with a SLURM scheduler
# Here is a link to the pipeline in use for further detail:
# https://github.com/AndersenLab/post-gatk-nf

source activate /data/eande106/software/conda_envs/nf24_env

cd ../../../processed_data/geo_info/
mkdir -p Taiwan
cd Taiwan

Taiwan_isotypes_GIS="../../../processed_data/geo_info/Taiwan/Ct_TW_isotype_GIS.txt"
tail -n +2 $Taiwan_isotypes_GIS | cut -f1 > Taiwan_GIS_strain_list.txt

strain_list="../../../processed_data/geo_info/Taiwan/Taiwan_GIS_strain_list.txt"
out_dir="../../../processed_data/geo_info/Taiwan"
raw_VCF_path="../../../data/VCF"
raw_VCF="../../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"

nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $raw_VCF \
--pops $strain_list \
--species c_tropicalis \
--vcf_folder $raw_VCF_path \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output $out_dir \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal
