# Ct_pop_gen_project

This repository contains scripts and files used to perform the analyses and make the figures and tables associated with the *Caenorhabditis tropicalis* population genetics manuscript 

- The "data" folder contains the raw datasets from public database (*e.g.* CaeNDR)
- The "plots" folder contains all the figures generated by the scripts
- The "processed_data" folder contains all the intermediate results files generated during the analysis
- The "scripts" folder contains all the codes used for generating figures and tables
- The "tables" folder contains all the tables generated by the scripts

<br>

---

## Main Figures & Tables

#### Figure 1. Global distribution of 518 *C. tropicalis* isotype reference strains
- Figure 1 a-f: Plot_isotype_tropicalis_global_map.R 
- Figure 1 g-i: Plot_PCA_vcf.R

#### Figure 2. The ML tree of 518 *C. tropicalis* isotype reference strains generated from variants pruned with LD r<sup>2</sup> values less than or equal to 0.9
- Figure 2: Plot_tree_in_lat_with_label.R

#### Figure 3. Diversity statistic calculated across global *C. tropicalis* isotypes
- Figure 3: Plot_tro_pi_theta_d.R

#### Figure 4. LD decay for all *C. tropicalis* isotype reference strains across all autosomes
- Figure 4: Plot_LD_geo_all.R

#### Figure 5. Hyper-variable regions (HVRs) in *C. tropicalis*
- Figure 5 b, c: Plot_HVRs_TajimasD.R

#### Figure 6. Relatedness and frequency of 518 *C. tropicalis* isotype reference strains at the Medea regions
- Figure 6: Plot_Medea_vcf2tree.R

#### Figure 7. Human Impact Index of nematode sampling sites
- Figure 7: Plot_human_footprint_index.R

#### Table 1. Effective population size (N<sub>e</sub>) for *Caenorhabditis* spp. samples from different studies
- Table 1: Table_isotypes_Ne_Outcrossing rate_20240420.R

#### Table 2. Outcrossing rate of isotype reference strains from different sampling sites
- Table 2: Table_isotypes_Ne_Outcrossing rate_20240420.R


<br>
<br>

## Supplemental Figures & Tables

#### Figure S1. Global distribution of 690 *C. tropicalis* strains
- Figure S1: Plot_strains_tropicalis_global_map.R

#### Figure S2. Geographical distance between strains within an isotype
- Figure S2 a: Plot_Vincenty_geo_Ce.R
- Figure S2 b: Plot_Vincenty_geo_Ct.R

#### Figure S3-S5. The ML trees of 518 *C. tropicalis* isotype reference strains generated from LD-pruned variants with r<sup>2</sup> values less than or equal to 0.7-0.9
- Figure S3-S5: Plot_tree_in_lat_with_label.R

#### Figure S6. Unrooted maximum likelihood trees of 518 *C. tropicalis* isotype reference strains generated from LD-pruned variants
- Figure S6: Plot_equal angle_unrooted_tree.R

#### Figure S7. Scatter plot shows significant positive correlation between geographic distance and phylogenetic distance
- Figure S7: Plot_phy_geo.R

#### Figure S8. Genetic distance (*D<sub>xy</sub>*) correlation with geographic distance and regional variation
- Figure S8: Plot_dxy_geo_change_in_dxy_per_km_20240518.R

#### Figure S9. Pairwise fixation index (*F<sub>st</sub>*) among *C. tropicalis* from different geographic regions
- Figure S9: Plot_fst_geo.R

#### Figure S10. LD decay for all autosomes
- Figure S10: Plot_LD_geo_all_chromosome.R

#### Figure S12-S16. Relatedness of 518 *C. tropicalis* isotype reference strains at five different maternal effect regions
- Figure S12-S16: Plot_TAs_vcf2tree.R

#### Table S4 Diversity statistic (*π*, *θ*, Tajima's D) for all isotype reference strains and for each geographic region
- Table S4: Table_geo_p_theta.R

#### Table S5 Effective population size (N<sub>e</sub>) for all isotype reference strains and for each geographic region except for the two under-sampled regions, Africa and Australia.
- Table S5: Table_isotypes_Ne_Outcrossing rate_20240420.R


