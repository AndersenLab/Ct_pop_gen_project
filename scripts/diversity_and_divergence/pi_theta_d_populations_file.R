rm(list=ls())

library(readr)
library(dplyr)

geo_info <- readr::read_csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv") %>% 
  dplyr::select(isotype,geo) %>%
  dplyr::filter(!is.na(isotype), !is.na(geo)) %>%
  dplyr::distinct()

table(geo_info$geo)

# population file by geo
write.table(geo_info,
            "../../processed_data/pi_theta_d/pop_file_geo.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

# population file for global estimates
global <- geo_info %>%
  mutate(geo = "GLOBAL")

table(global$geo)

write.table(global,
            "../../processed_data/pi_theta_d/pop_file_GLOBAL.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')
