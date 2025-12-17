#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --output=Calculate_nf_Ct_Hawaii_PCA.oe
#SBATCH --job-name="CtnfHaPCA"





cd $HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data
mkdir -p Hawaii
cd Hawaii



Hawaii_isotypes_GIS="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Hawaii/Ct_HW_isotype_GIS.txt"
tail -n +2 $Hawaii_isotypes_GIS | cut -f1 > Hawaii_GIS_strain_list.txt

strain_list="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Hawaii/Hawaii_GIS_strain_list.txt"
out_dir="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/processed_data/Hawaii"
raw_VCF_path="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF"
raw_VCF="$HOME/vast-eande106/projects/Bowen/PopGen_Tro_Project/2025_PopGen_Tro/data/VCF/WI.20250331.hard-filter.isotype.vcf.gz"



module load python/anaconda
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




