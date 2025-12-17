




############# copy these scripts and run them on the login nodes ########



tmux new -s PCA

#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env




cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Taiwan
cd Taiwan



Taiwan_isotypes_GIS="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Taiwan/Ct_TW_isotype_GIS.txt"
tail -n +2 $Taiwan_isotypes_GIS | cut -f1 > Taiwan_GIS_strain_list.txt

strain_list="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Taiwan/Taiwan_GIS_strain_list.txt"
out_dir="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Taiwan"
raw_VCF_path="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF"
raw_VCF="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF/WI.20250627.hard-filter.isotype.vcf.gz"




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




#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t PCA



