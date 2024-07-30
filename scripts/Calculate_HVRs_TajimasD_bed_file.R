# 

rm(list = ls())

library(dplyr)
library(ggplot2)
library(purrr)


HVR_table_raw<-read.table("../processed_data/HVRs/HDR_5kbclust_collapsed_wFreq.tsv",
                          header = TRUE)
HVR_bed<-HVR_table_raw %>% 
  select(CHROM,minStart,maxEnd) %>% 
  rename(chr=CHROM, strat=minStart, stop=maxEnd)

write.table(HVR_bed,
            "../processed_data/pi_theta_d/HVRs_pi_theta_d/CT_HVRs_240617.bed",
            sep = '\t',
            quote = FALSE,
            row.names = FALSE)




