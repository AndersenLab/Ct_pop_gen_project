rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

het_data<-readr::read_tsv("../../processed_data/HDR_stats/Ct_soft_filtered_hets.het_counts_per_site.HDR_labeled.tsv")

head(het_data)
dim(het_data)

###### Check number of SNV sites ######
table(het_data$region)
sum(het_data$region == "HDR")
sum(het_data$region == "non-HDR")

###### Add non-heterozygous called genotype counts ######
het_data<-het_data %>%
  dplyr::mutate(nonhet_called = n_called - n_het)

###### Summary for HDR and non-HDR ######
summary_region<-het_data %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    total_SNV_sites = n(),
    sites_with_at_least_one_het = sum(has_het == 1),
    percent_sites_with_het = sites_with_at_least_one_het / total_SNV_sites * 100,
    total_het_genotype_calls = sum(n_het),
    total_called_genotypes = sum(n_called),
    total_nonhet_called_genotypes = sum(nonhet_called),
    het_call_rate_percent = total_het_genotype_calls / total_called_genotypes * 100,
    mean_het_calls_per_site = mean(n_het),
    median_het_calls_per_site = median(n_het)
  )

###### Genome-wide summary ######
summary_genomewide<-het_data %>%
  dplyr::summarise(
    region = "genome_wide",
    total_SNV_sites = n(),
    sites_with_at_least_one_het = sum(has_het == 1),
    percent_sites_with_het = sites_with_at_least_one_het / total_SNV_sites * 100,
    total_het_genotype_calls = sum(n_het),
    total_called_genotypes = sum(n_called),
    total_nonhet_called_genotypes = sum(nonhet_called),
    het_call_rate_percent = total_het_genotype_calls / total_called_genotypes * 100,
    mean_het_calls_per_site = mean(n_het),
    median_het_calls_per_site = median(n_het)
  )

###### Combine summary table ######
summary_all<-dplyr::bind_rows(summary_genomewide, summary_region)
summary_all

