# Processed Data

This directory contains the processed datasets used for downstream population genomic analyses. Each subdirectory corresponds to a major analysis component. The files are derived from upstream pipelines (see `scripts/` and the main repository README for details).

## Directory Structure

processed_data/
- Geo_info/
- pi_theta_d/
- pi_theta_d_geo/
- Ct_admixture_k28/
- Ct_admixture_best_k/
- Ct_admixture/
- LD_pruned_trees/LD_0.9/
- Ct_pruned_VCF_and_PCA/
- PCA_by_chrom/
- heatmap_hard_filtered/


## Geo_info/

Processed files summarizing geographic info of the samples.

### Contains

- `Ct_indep_isotype_info_geo.csv`  
  Geographic information for each independent isotype.

- `Ct_isotype_geo_freq.csv`  
  Frequency of each isotype across geographic regions, with cosmopolitan strains treated as a group.

- `Ct_lineage_all.csv`  
  A table mapping isotypes to relatedness groups.

- `geo_and_lineage.csv`  
  A table showing geographic region, coordinates, and relatedness group of each isotype.

- `non_admixed_isotype_replicate_*.txt`  
  A table showing subpopulation assignment by ADMIXTURE, geographic region, coordinates, and relatedness group of each isotype.


## pi_theta_d/

Nucleotide diversity metrics per genomic window.


## pi_theta_d_geo

Nucleotide diversity metrics per genomic window for each geographic region



## Ct_admixture_k28

admixture data at best k = 28

### Contains

- `best_k_long_admix_pops_replicates_*.csv`  
  Final selected ADMIXTURE results with geography info for each isotype

- `K28_best_k_long_admix_pops_replicates_*.csv`  
  Final selected ADMIXTURE results with relatedness groups info for each isotype

- `K28_Processed_Ancestry_replicate_*.tsv`  
  ADMIXTURE Q-matrix files showing the ancestry proportion of each isotype across inferred genetic clusters


## Ct_admixture_best_k

admixture data for finding best K

### Contains

- `concat_Qfiles_K*.tsv`  
  Q-matrices for each value of K tested.

- `cv_matrix.tsv`  
  Cross-validation error matrix.



## Ct_admixture

Seeds, CV error, and admixture data for K from 26-30

### Contains

- `seeds.txt`  
  list of seeds used in admixture analysis
  
- `admix_replicates_CV.tsv`  
  A list of Cross-validation error.

- `LD_0.9_*_*.Q`
  Q-files for each value of K (26-30) tested.

## LD_pruned_trees/LD_0.9

LD-pruned tree files.

### Contains

- `phy_file_LD_0.9.phy.contree`  
  Consensus phylogenetic tree of all isotypes.

- `phy_file_LD_0.9.phy.contree`  
  Consensus phylogenetic tree of all isotypes (midpoint-rooted).


## Ct_pruned_VCF_and_PCA

Processed files summarizing PCA.

### Contains

- `sample_list.txt`  
  A list of all 622 isotypes used in this study.

- `TracyWidom_statistics_no_removal.tsv`
  TracyWidom statistics results for PCA analysis of all 622 isotypes.

- `eigenstrat_no_removal.evac`
  Eigenstrat values results for PCA analysis of all 622 isotypes.


## Ct_pruned_VCF_and_PCA

Processed files summarizing PCA by each chromosome.
### Contains

- `*/TracyWidom_statistics_no_removal.tsv`
  TracyWidom statistics results for PCA analysis by each chromosome.

- `*/eigenstrat_no_removal.evac`
  Eigenstrat values results for PCA analysis by each chromosome.


## Hawaii

Processed files summarizing PCA and GIS data for Hawaii.


## Taiwan

Processed files summarizing PCA and GIS data for Taiwan.


## heatmap_hard_filtered
Processed files summarizing genetic similarity estimates

### Contains

- `gtcheck.tsv`  
  Pairwise isotype genetic similarity estimates for all isotypes.

- `concordance_coverage_sample_sheet.txt`
  A list of all isotypes for calculating genetic similarity.


