# This script need to be run in a cluster with a SLURM scheduler
# Here is a link to the pipeline in use for further detail:
# https://github.com/AndersenLab/post-gatk-nf

cd ../../processed_data/PCA_by_chrom
mkdir -p ${chr}
cd ${chr}

source activate /data/eande106/software/conda_envs/nf24_env

nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf ../../processed_data/PCA_by_chrom/by_chrom_vcfs/${chr}.vcf.gz \
--pops ../../processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt \
--species c_tropicalis \
--vcf_folder ../../processed_data/PCA_by_chrom/by_chrom_vcfs \
--eigen_ld 0.9 \
--postgatk false \
--singletons false \
--delly false \
--output ../../processed_data/PCA_by_chrom/${chr} \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal
