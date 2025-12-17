########################################################################
############################## PCA all ###################################
########################################################################

############# copy these scripts and run them on the login nodes ########



### generate link
ln -s $HOME/vast-eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF
ln -s $HOME/vast-eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz.tbi $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF



raw_VCF_path="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF"
raw_VCF="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"



source activate CT_PopGen

cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/
mkdir -p Ct_pruned_VCF_and_PCA
cd Ct_pruned_VCF_and_PCA

bcftools query -l "$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF/WI.20250627.hard-filter.isotype.vcf.gz" > sample_list.txt





#####
tmux new -s PCA

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env

cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_pruned_VCF_and_PCA


nextflow run -latest andersenlab/post-gatk-nf/main.nf \
-r delly --pca true \
--pca_vcf $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF/WI.20250627.hard-filter.isotype.vcf.gz \
--pops $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt \
--species c_tropicalis \
--vcf_folder $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF \
--eigen_ld 0.9,0.8,0.7,0.6 \
--postgatk false \
--singletons false \
--delly false \
--output $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Ct_pruned_VCF_and_PCA \
--eigen_par_outlier_removal eigpar_outliers_removed \
--eigen_par_no_removal eigpar_no_removal



#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA
