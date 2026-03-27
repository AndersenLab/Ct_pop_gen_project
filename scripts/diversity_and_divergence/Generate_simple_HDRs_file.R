rm(list = ls())

library(dplyr)
library(readr)


HDRs_raw <- read_tsv("../../tables/HDR_CT_allStrain_5kbclust_20251201.tsv", show_col_types = FALSE)

HDRs_file<-HDRs_raw %>% 
  select(CHROM,minStart,maxEnd)

write_tsv(HDRs_file,
          "../../processed_data/pi_theta_exclude_HDRs/HDRs_file.tsv",
          col_names = FALSE)





