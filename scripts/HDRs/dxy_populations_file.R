rm(list=ls())

library(readr)
library(dplyr)

# generate populations_file
lineage<- readr::read_csv("../../processed_data/geo_info/geo_and_lineage.csv") %>% 
  select(isotype,lineage)

table(lineage$lineage)

dir.create("../../processed_data/Dxy_LAC/")
write.table(lineage,"../../processed_data/Dxy_LAC/pop_file_lineage.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

