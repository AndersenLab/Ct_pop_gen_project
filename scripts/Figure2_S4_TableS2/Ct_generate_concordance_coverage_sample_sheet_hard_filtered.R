rm(list = ls())

library(dplyr)

# read .txt file
strain_list <- read.table("../../processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt",
                          header = FALSE, stringsAsFactors = FALSE)



result <- strain_list %>%
  dplyr::mutate(
    bam = paste0(V1, ".bam"),
    bai = paste0(V1, ".bam.bai"),
    coverage = 0,
    percent_mapped = 0
  ) %>% 
  dplyr::rename("strain"="V1")

write.table(result, "../../processed_data/heatmap_hard_filtered/concordance_coverage_sample_sheet.txt", sep = "\t", row.names = FALSE, quote = FALSE)

