rm(list = ls())

library(tidyverse)
library(readr)
library(dplyr)
library(stringr)

infile <- "../../processed_data/gene_enrichment/N2_NIC58_orthogroups.tsv"
orthogroups <- readr::read_tsv(infile, show_col_types = FALSE)
colnames(orthogroups)

ce_col <- "c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein"
ct_col <- "c_tropicalis.NIC58_20251002.csq.longest.protein"

orthogroups_1to1 <- orthogroups %>%
  dplyr::mutate(
    ce_n = if_else(
      is.na(.data[[ce_col]]) | .data[[ce_col]] == "",
      0L,
      str_count(.data[[ce_col]], ",") + 1L
    ),
    ct_n = if_else(
      is.na(.data[[ct_col]]) | .data[[ct_col]] == "",
      0L,
      stringr::str_count(.data[[ct_col]], ",") + 1L
    )
  ) %>%
  dplyr::filter(ce_n == 1, ct_n == 1) %>%
  dplyr::transmute(
    Orthogroup,
    Ce_transcript = str_trim(.data[[ce_col]]),
    Ct_transcript = str_trim(.data[[ct_col]])
  )

out_file <- "../../processed_data/HDR_stats/N2_NIC58_1to1_orthologs.tsv"

readr::write_tsv(orthogroups_1to1, out_file)
nrow(orthogroups_1to1)





