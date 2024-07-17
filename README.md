# Ct_pop_gen_project

This repository contains scripts and files used to perform the analyses and make the figures and tables associated with the *Caenorhabditis tropicalis* population genetics manuscript 

- The "data" folder contains the raw datasets from public database (*e.g.* CaeNDR)
- The "plots" folder contains all the figures generated by the scripts
- The "processed_data" folder contains all the intermediate results files generated during the analysis.
- The "scripts" folder contains all the codes used for generating figures and tables.
- The "tables" folder contains all the tables generated by the scripts

---

### Main Figures & Tables

#### Figure 1. Global distribution of 518 *C. tropicalis* isotype reference strains.
- Figure 1 a-f: Plot_isotype_tropicalis_global_map.R 
- Figure 1 g-i: Plot_PCA_vcf.R

#### Figure 2. The ML tree of 518 *C. tropicalis* isotype reference strains generated from variants pruned with LD r<sup>2</sup> values less than or equal to 0.9.
- Figure 2: Plot_tree_in_lat_with_label.R

#### Figure 3. Diversity statistic calculated across global *C. tropicalis* isotypes.
- Figure 3: Plot_tro_pi_theta_d.R

#### Figure 4. LD decay for all *C. tropicalis* isotype reference strains across all autosomes. 
- Figure 4: Plot_LD_geo_all.R

#### Figure 5. Hyper-variable regions (HVRs) in *C. tropicalis*.
- Figure 5 b, c: Plot_HVRs_TajimasD.R

#### Figure 6. Relatedness and frequency of 518 *C. tropicalis* isotype reference strains at the Medea regions.
- Figure 6: Plot_Medea_vcf2tree.R

#### Figure 7. Human Impact Index of nematode sampling sites. 
- Figure 7: Plot_human_footprint_index.R

#### Table 1. Effective population size (Ne) for *Caenorhabditis* spp. samples from different studies.​​
- Table 1: Table_isotypes_Ne_Outcrossing rate_20240420.R

#### Table 2. Outcrossing rate of isotype reference strains from different sampling sites
- Table 2: Table_isotypes_Ne_Outcrossing rate_20240420.R


<br>
<br>

### Supplemental Figures & Tables

#### Figure S1. Global distribution of 690 *C. tropicalis* strains.
- Figure S1: Plot_strains_tropicalis_global_map.R

#### Figure S2. Geographical distance between strains within an isotype.
- Figure S2 a: Plot_Vincenty_geo_Ce.R
- Figure S2 b: Plot_Vincenty_geo_Ct.R

#### Figure S3-S5. The ML tree of 518 *C. tropicalis* isotype reference strains generated from LD-pruned variants with r<sup>2</sup> values less than or equal to 0.7-0.9.
- Figure S3-S5: Plot_tree_in_lat_with_label.R

#### 









