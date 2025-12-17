rm(list = ls())

library(dplyr)

geo_raw<-read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv")

lineage_raw<-read.csv("../../processed_data/geo_info/Ct_lineage_all.csv")


Taiwan_geo_and_linegae<-geo_raw %>% 
  filter(geo %in% c("Taiwan")) %>% 
  select(-strain) %>% 
  left_join((lineage_raw %>% rename(isotype=sample)), by = c("isotype"))


write.table(Taiwan_geo_and_linegae,
            "../../processed_data/Taiwan/Ct_TW_isotype_GIS.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')  


