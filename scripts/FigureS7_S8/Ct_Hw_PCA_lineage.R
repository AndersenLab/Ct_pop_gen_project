rm(list=ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(dendextend)
library(ggpubr)
library(cowplot)


# source
source("../utilities.R")



### input GIS data
GIS_data_raw<-read.table("../../processed_data/Hawaii/Ct_HW_isotype_GIS.txt",
                            sep = '\t',header = TRUE)

GIS_data <- GIS_data_raw %>% 
  dplyr::mutate(collection_island = 
                  ifelse(latitude >= 21.864700 & latitude <= 22.253097 &
                           longitude >= -159.823456 & longitude <= -159.242411,"Kauai",
                         ifelse(latitude >= 21.208833 & latitude <= 21.766432 &
                                  longitude >= -158.340215 & longitude <= -157.585817,"Oahu",
                                ifelse(latitude >= 21.025107 & latitude <= 21.262047 &
                                         longitude >= -157.348973 & longitude <= -156.680570,"Molokai",
                                       ifelse(latitude >= 20.554665 & latitude <= 21.055465 &
                                                longitude >= -156.707098 & longitude <= -155.936584,"Maui",
                                              ifelse(latitude >= 18.880884 & latitude <= 20.360849 &
                                                       longitude >= -156.139088 & longitude <= -154.676093,"Big Island",NA)))))) %>% 
  dplyr::mutate(`island age Ma` = case_when(collection_island == "Big Island" ~ 0.6,
                                            collection_island == "Maui" ~ 1,
                                            collection_island == "Molokai" ~ 1.8,
                                            collection_island == "Oahu" ~ 3.2,
                                            collection_island == "Kauai" ~ 5.1)) %>% 
  rename(isotype = strain)




# calculate % var explained. 
EVarNoRm_ld0.9 <- data.table::fread("../../processed_data/Hawaii/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") %>%
  dplyr::mutate(sum = sum(eigenvalue),
                VarExp = eigenvalue/sum,
                sigEV = ifelse(`p-value` < 0.05, T, F)) %>%
  dplyr::group_by(sigEV) %>%
  dplyr::mutate(sigVarExp = sum(VarExp)) %>%
  dplyr::ungroup()






cdf9no <- data.table::fread("../../processed_data/Hawaii/EIGESTRAT/LD_0.9/NO_REMOVAL/eigenstrat_no_removal.evac", skip = 1)


cdf9_2no <- cdf9no %>%
  dplyr::select(isotype = V1, PC1=V2, PC2=V3, PC3=V4, PC4=V5, PC5=V6, PC6 = V7)


#### add lineage info ###
lineage<-read.csv("../../processed_data/geo_info/Ct_lineage_all.csv") %>% 
  rename(isotype=sample)

cdf9_7no <- cdf9_2no %>%
  dplyr::full_join(GIS_data,by=c("isotype")) %>% 
  dplyr::left_join(lineage, by = c ("isotype"))

cdf9_7no$lineage %>% unique()



pca_TAC_ld0.9_no_rm <-  cdf9_7no


saveRDS(cdf9_7no, file = "../../processed_data/Hawaii/cdf9_7no.rds")




Evout <- EVarNoRm_ld0.9 %>%
  dplyr::select(N, VarExp) %>%
  dplyr::mutate(vperc = VarExp * 100) %>%
  dplyr::slice_head(n = 4)
Evout

# plot PC1 vs PC2
pc12 <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC1, y = PC2, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 14, color = "black")) + # presentation text
  labs(x = paste0("PC1(",round(Evout$vperc[1],2),"%)"), y = paste0("PC2(",round(Evout$vperc[2],2),"%)"),
       color = "cluster", shape = "island", 
       title = "hard-filtered LD=50-10-0.9\n n=56/56, ot=6.0")
pc12

pc12_blank <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC1, y = PC2, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21, stroke = 0.15) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 9, color = "black"), legend.position = "none") + 
  labs(x = paste0("PC1(",round(Evout$vperc[1],2),"%)"), y = paste0("PC2(",round(Evout$vperc[2],2),"%)"),color = "cluster", title = "")
pc12_blank

pc34 <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC3, y = PC4, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 14, color = "black"), legend.position = "right") + # presentation text
  labs(x = paste0("PC3(",round(Evout$vperc[3],2),"%)"), y = paste0("PC4(",round(Evout$vperc[4],2),"%)"),
       color = "cluster", 
       title = "hard-filtered LD=50-10-0.9\n n=56/56, ot=6.0")
pc34

pc34_blank <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC3, y = PC4, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21, stroke = 0.15) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 9, color = "black"), legend.position = "none") + 
  labs(x = paste0("PC3(",round(Evout$vperc[3],2),"%)"), y = paste0("PC4(",round(Evout$vperc[4],2),"%)"),color = "cluster", title = "")
pc34_blank


pc12_blank <- pc12_blank + guides(fill = guide_legend(title = "Group", ncol = 5))
pc34_blank <- pc34_blank + guides(fill = guide_legend(title = "Group", ncol = 5))

PCA_main <- ggarrange(
  pc12_blank, pc34_blank,
  labels = c("A", "B"),          
  label.x = 0.01, label.y = 0.95, 
  legend = "bottom",            
  common.legend = TRUE          
)

PCA_main
# ggsave("Ct_PCA_Hawaii_Lineage.pdf",
#        PCA_main,
#        width = 7.5, height = 4, units = "in")
saveRDS(PCA_main, file = "../../processed_data/Hawaii/Ct_PCA_Hawaii_Lineage.rds")





