# Scripts

This directory contains the scripts used to perform analyses and generate figures and tables for the C. tropicalis manuscript.

## Directory Structure

scripts/
- geography_trees_PCA_LD/
- genetic_similarity_and_admixture/
- diversity_and_divergence/
- HDRs/
- gene_enrichment/
- misc_utilities/

## Abbreviations

SF = Supplementary Figure  
SD = Supplementary Data  

---

## geography_trees_PCA_LD

This directory contains scripts for geography, PCA, trees, and linkage disequilibrium (LD) analyses.

### Core pipeline scripts

- `Calculate_Ct_vcf_to_pyh.sh`  
Converts VCF files to PHYLIP format for tree construction.

- `Calculate_Ct_pyh_to_tree.sh`  
Generates trees from PHYLIP alignments.

- `Trees.R`  
Plots trees.  
Generates Figure S4 and Figure S5.

- `nf_Ct_pruned_VCF_and_PCA.sh`  
Generates LD-pruned VCF files and performs PCA.

- `Calculate_Ct_pca_by_chrom.sh`  
Generates per-chromosome PCA input files.

- `nf_Ct_pca_by_chrom.sh`  
Runs PCA on per-chromosome VCF files and generates eigenvalues and PC scores.

- `PCA_by_chrom_LD_0.9.R`  
Plots chromosome-specific PCA results.  
Generates Figure S6.

- `PCA_figure_and_table.R`  
Generates PCA summary plots and associated TableS3.

- `Figure1.R`  
Assembles panels for Figure 1.

- `FigureS2.R`  
Assembles panels for Figure S2.

### Geographic visualization

- `Map_all_three_species.R`  
Plots geographic distributions across species.  
Generates Figure S1.

- `Geo_locations_isotypes.R`  
Plots geographic locations of isotypes.

- `Geo_locations_strains.R`  
Plots geographic locations of all strains.

- `Geo_distance_number_of_starins.R`  
Summarises geographic distances and number of strains across isotypes.

### LD analyses

- `LD_decay_three_species.R`  
Summarise LD decay across species.  
Generates Figure S3.

- `Calculate_Ct_LD_per_Mb_all.sh`  
Computes LD statistics for C. tropicalis.

- `Calculate_Cb_LD_per_Mb_all.sh`  
Computes LD statistics for C. briggsae.

- `Calculate_Ce_LD_per_Mb_all.sh`  
Computes LD statistics for C. elegans.

### Regional analyses

These subdirectories contain scripts for region-specific analyses, including geographic mapping, PCA, and genetic similarity visualization where applicable.

#### Africa/

- `Af_RGs_map.R`  
Plots relatedness group (RG) distribution in Africa.  
Generates Figure S26.

- `Af_similarity.R`  
Analyzes genetic similarity among African isotypes.  
Generates Figure S27.

#### Central_America/

- `CentralAmerica_RGs_map.R`  
Plots RG distribution in Central America.  
Generates Figure S19.

- `Central_America_genetic_similarity.R`  
Analyzes genetic similarity within Central American isotypes.  
Generates Figure S20.

#### Caribbean/

- `Caribbean_RGs_map.R`  
Plots RG distribution in Caribbean isotypes.  
Generates Figure S18.

#### South_America/

- `South_America_RGs_map.R`  
Plots RG distribution in South America.  
Generates Figure S21.

- `South_America_heatmap.R`  
Generates genetic similarity heatmap for South American isotypes.  
Generates Figure S22.

#### Hawaii/

- `Hw_PCA_RGs.R`  
Plot PCA and visualizes RG in Hawaiian isotypes.  

- `Hw_PCA_by_env.R`  
Analyzes PCA patterns in relation to environmental variables.  
Generates Figure S17.

- `Hw_PCA_group_map.R`  
Assemble PCA plot and generate map of RGs.
Generates Figure S16.

- `nf_Hawaii_PCA.sh`  
Pipeline script for PCA analysis in Hawaii dataset.

#### Taiwan/

- `Tw_PCA_RGs.R`  
Performs PCA and visualizes RG structure in Taiwanese isotypes.  


- `Tw_PCA_by_env.R`  
Analyzes PCA patterns in relation to environmental variables.  
Generates Figure S15.

- `Tw_PCA_group_map.R`  
Assemble PCA plot and generate map of RGs.
Generates Figure S14.

- `Taiwan_GIS.R`  
Sumarrise GIS data of Taiwanese sampling locations.

- `nf_Taiwan_PCA.sh`  
Pipeline script for PCA analysis in Taiwan dataset.

#### Micronesia/

- `Micronesia_RGs_map.R`  
Plots RG distribution in Micronesia.  
Generates Figure S23.

- `Micronesia_similarity.R`  
Analyzes genetic similarity among Micronesian isotypes.  
Generates Figure S24.

#### Indo/

- `Indo_RG_map.R`  
Plots RG distribution in Indonesia (Malay Archipelago).  
Generates Figure S25.

---

## genetic_similarity_and_admixture

This directory contains scripts for genetic similarity estimation, ADMIXTURE analyses, and isolation-by-distance analyses.

- `Calculate_Ct_admixture_prepare_files.sh`  
Prepares input files for ADMIXTURE analysis and generates random seeds for replicated runs.

- `Calculate_Ct_run_admixture_arrays.sh`  
Runs ADMIXTURE across a range of K values and seeds from 1 to 25 using the LD-pruned dataset.

- `Calculate_Ct_run_admixture_arrays_26_30.sh`  
Runs ADMIXTURE for selected K values from 26 to 30 using the LD-pruned dataset.

- `Calculate_Ct_run_admixture_post_processing.sh`  
Processes ADMIXTURE log files and extracts cross-validation (CV) errors into a summary table (`admix_replicates_CV.tsv`).

- `Calculate_Ct_admixture_best_k.sh`  
Collects ADMIXTURE output files, links the CV summarization script, and generates the CV matrix used to identify the optimal K.

- `generate_CV_matrix.sh`  
Parses ADMIXTURE log files and compiles cross-validation errors into `cv_matrix.tsv`.

- `concat_Qs.sh`  
Concatenates ADMIXTURE Q matrices across runs and writes group-specific Q files for downstream plotting.

- `Admixture_best_k.R`  
Plots ADMIXTURE CV errors across K values and identifies the optimal K.
Generate Figure S8.

- `Admixture_K26K30_rep5.R`  
Visualizes ADMIXTURE results for selected K values between K = 26 and K = 30.
Geenrate Figure S9.

- `Admixture_by_lRGs.R`  
Plots ADMIXTURE assignment patterns grouped by relatedness groups (RGs).  
Generates Figure S11.

- `Genetic_similarity_heatmap_hard_filtered.R`  
Calculates and visualizes pairwise genetic similarity among isotypes using hard-filtered variants, and exports the associated geographic and RGs metadata for downstream analyses.
Generate Figure 2 and Figure S7.

- `Similarity_histogram_hard_filtered.R`  
Plots the distribution of pairwise genetic similarity values across strains.
Generate Figure S37.

- `Isolation_by_distance.R`  
Analyzes the relationship between geographic distance and genetic distance among isotypes.
Generate Figure S12 and Figure S13.

---

## diversity_and_divergence

This directory contains scripts for calculating nucleotide diversity (π), Watterson’s θ, Tajima’s D, and Dxy, including analysis for HDRs or nonHDRs.

- `Ct_vcf_to_zarr.sh`  
Converts the genome-wide VCF file to Zarr format for downstream calculation of π, θ, and Tajima’s D.

- `Ct_zarr_to_pi_theta_d.sh`  
Calculates π, Watterson’s θ, and Tajima’s D from genome-wide Zarr files.

- `Ct_vcf_HDRs.sh`  
Generates VCF files restricted to hyper-divergent regions (HDRs).

- `Ct_vcf_nonHDRs.sh`  
Generates VCF files restricted to nonHDR regions.

- `Ct_vcf_to_zarr_nonHDRs.sh`  
Converts nonHDR VCF files to Zarr format for downstream diversity analyses.

- `Ct_zarr_to_pi_theta_d_nonHDRs.sh`  
Calculates π, Watterson’s θ, and Tajima’s D from nonHDR Zarr files.

- `Ct_vcf_nonHDRs_geo.sh`  
Generates nonHDR VCF files for geo-specific analyses.

- `vcf_to_zarr_geo.sh`  
Converts geo-specific VCF files to Zarr format.

- `zarr_to_pi_theta_d_geo.sh`  
Calculates π, Watterson’s θ, and Tajima’s D from Zarr files for geographic groups.

- `vcf_to_zarr_nonHDRs_geo.sh`  
Converts geo-specific nonHDR VCF files to Zarr format.

- `zarr_to_pi_theta_d_nonHDRs_geo.sh`  
Calculates π, Watterson’s θ, and Tajima’s D from nonHDR Zarr files for geographic groups.

- `Calculate_geo_vcf.sh`  
Generates VCF files for analyses of diversity and divergence across geographic groups.

- `Geo_pi_theta_d.R`  
Visualizes π, Watterson’s θ, and Tajima’s D across geographic groups.
Generate Table S4. 

- `Geo_pi_theta_d_Autosomes_X.R`  
Summarises diversity statistics across autosomes and the X chromosome for geographic groups.
Generate Table S5. 

- `pi_theta_nonHDRs.R`  
Analyzes π, Watterson’s θ, and Tajima’s D in nonHDR regions.
Generate Figure S35. 

- `pi_theta_d_fold_change_nonHDR_geo.R`  
Calculates fold changes in diversity statistics across geographic groups in nonHDR regions.
Generate Table S8. 

- `Calculate_Dxy_generate_VCF.sh`  
Generates VCF files required for pairwise Dxy calculations.

- `Generate_Dxy_pop_files.R`  
Generates population assignment files used for pairwise Dxy calculations among relatedness groups or geographic groups.

- `Calculate_Dxy_LAC.sh`  
Calculates pairwise Dxy among relatedness groups.

- `Dxy_HDRs_nonHDRs.R`  
Compares Dxy between relatedness groups within HDR or nonHDR regions.
Generate Figure S32. 

- `Generate_simple_HDRs_file.R`  
Generates simplified HDR annotation files for downstream diversity and divergence analyses.

- `characterize_HDRs.R`  
Summarize HDR statistics across wild isolate genomes.

- `Similarity_HDRs.sh`  
Calculates pairwise genetic similarity using variants located in HDRs.

- `Similarity_nonHDRs.sh`  
Calculates pairwise genetic similarity using variants located outside HDRs.

- `Heatmap_similarity_HDRs.R`  
Visualizes pairwise genetic similarity in HDR regions as a heatmap.
Generate Figure S33. 

- `Heatmap_similarity_nonHDRs.R`  
Visualizes pairwise genetic similarity in nonHDR regions as a heatmap.
Generate Figure S34. 

- `pi_theta_d_python/`  
Core Python scripts for calculating π, Watterson’s θ, and Tajima’s D.

---

## HDRs

This directory contains scripts for identifying and characterizing HDRs.

- `call_HDRs.R`  
Identifies HDRs across genomes.

- `TAs_overlap_HDRs.R`  
Analyzes overlap between TAs and HDRs.
Generate Figure S36.

- `visualize_spp_alignments.R`
Visualizes alignment data between divergent strain pairs across all three selfing _Caenorhabditis_ species

---

## gene_enrichment
This directory contains scripts for annotating the NIC58 proteome using InterProScan and performing functional enrichment analysis of IPR terms and Gene Ontology gene sets.

- `InterProScan.sh`
This script annotates NIC58 proteins with InterProScan terms, which outputs a table used in `IPR_GO_tropicalis.R`

- `IPR_GO_tropicalis.R`
Given the NIC58 reference genome and InterProsScan annotations of the NIC58 proteome, this script classifies genes that are found in the chromosomal arm domains, genes found in hyper-divergent regions in chromosomal arm domains, and performs a one-sided hyper-geometric test to identify statistically-enriched InterProScan functional protein domains in genes found in hyper-divergent regions in chromosomal arm domains in relation to genes found outside of hyper-divergent regions in chromosomal arm domains.

Generates Figure 4 and Table S10.

---

## misc_utilities

Utility scripts used throughout the pipeline.

- `getVariantCounts.sh`  
Counts variants from VCF files.

- `merge_count_files.sh`  
Merges count outputs.

- `getMosCov.sh`, `organize_mosdepth_files.sh`  
Processes coverage statistics.

- `gen_sample_list.sh`  
Generates sample lists.

- `merge_coords.sh`, `get_bins.sh`  
Processes genomic coordinates.

- `nucmer.sh`  
Performs genome alignment.

---

## utilities.R

Contains shared R functions used across multiple scripts.

