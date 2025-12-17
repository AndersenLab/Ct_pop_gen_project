
############# copy these scripts and run them on the login nodes ########


##### use tmux
tmux new -s heatmap_HF

cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/heatmap_hard_filtered


#### module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env


nextflow run -latest andersenlab/concordance-nf \
--vcf="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF/WI.20250627.hard-filter.isotype.vcf.gz" \
--bam_coverage="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/heatmap_hard_filtered/concordance_coverage_sample_sheet.txt" \
--species=c_tropicalis

#### add -resume if necessary


##### then: ctrl-B, d

# To reattach (reopen) the detached session,
tmux attach -t heatmap_HF









