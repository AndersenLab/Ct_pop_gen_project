rm(list=ls())

library(readr)
library(dplyr)

lineage<- readr::read_csv("../../processed_data/geo_info/geo_and_lineage.csv") %>% 
  select(isotype,lineage)

LAC_Hw3<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Hw3"))

LAC_Af<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Af"))

LAC_HC<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","HC"))

LAC_Hw1<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Hw1"))

LAC_Mic1<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Mic1"))

LAC_Mic2<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Mic2"))

LAC_Tw1<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Tw1"))

LAC_Tw3<-lineage %>% 
  dplyr::filter(lineage %in% c("LAC","Tw3"))


dir.create("../../processed_data/Dxy_LAC/")
write.table(LAC_Hw3,"../../processed_data/Dxy_LAC/LAC_Hw3.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_Af,"../../processed_data/Dxy_LAC/LAC_Af.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_HC,"../../processed_data/Dxy_LAC/LAC_HC.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_Hw1,"../../processed_data/Dxy_LAC/LAC_Hw1.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_Mic1,"../../processed_data/Dxy_LAC/LAC_Mic1.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_Mic2,"../../processed_data/Dxy_LAC/LAC_Mic2.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_Tw1,"../../processed_data/Dxy_LAC/LAC_Tw1.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(LAC_Tw3,"../../processed_data/Dxy_LAC/LAC_Tw3.txt",
            quote = FALSE,  
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')


