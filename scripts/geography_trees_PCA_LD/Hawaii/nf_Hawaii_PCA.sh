########################################################################
############################## PCA all ###################################
########################################################################

# This script need to be run in a cluster with a SLURM scheduler
# Here is a link to the pipeline in use for further detail:
# https://github.com/AndersenLab/post-gatk-nf

cd ../../../processed_data/geo_info
mkdir -p Hawaii
cd Hawaii

Hawaii_isotypes_GIS="../../../processed_data/geo_info/Hawaii/Ct_HW_isotype_GIS.txt"
tail -n +2 $Hawaii_isotypes_GIS | cut -f1 > Hawaii_GIS_strain_list.txt

strain_list="../../../processed_data/geo_info/Hawaii/Hawaii_GIS_strain_list.txt"
out_dir="../../../processed_data/geo_info/Hawaii"
raw_VCF_path="../../../data/VCF"
raw_VCF="../../../data/VCF/WI.20250331.hard-filter.isotype.vcf.gz"

source activate /data/eande106/software/conda_envs/nf24_env

nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $raw_VCF \
--pops $strain_list \
--species c_tropicalis \
--vcf_folder $raw_VCF_path \
--eigen_ld 0.9 \
--postgatk false \
--delly false \
--output $out_dir
