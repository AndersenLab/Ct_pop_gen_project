########################################################################
############################## PCA all ###################################
########################################################################

# This script need to be run in a cluster with a SLURM scheduler
# Here is a link to the pipeline in use for further detail:
# https://github.com/AndersenLab/post-gatk-nf

source activate CT_PopGen

cd ../../processed_data/
mkdir -p Ct_pruned_VCF_and_PCA
cd Ct_pruned_VCF_and_PCA

bcftools query -l "../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz" > sample_list.txt

source activate /data/eande106/software/conda_envs/nf24_env

cd ../../processed_data/Ct_pruned_VCF_and_PCA

nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf ../../data/VCF/WI.20250627.hard-filter.isotype.vcf.gz \
--pops ../../processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt \
--species c_tropicalis \
--vcf_folder ../../data/VCF \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output ../../processed_data/Ct_pruned_VCF_and_PCA \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal
