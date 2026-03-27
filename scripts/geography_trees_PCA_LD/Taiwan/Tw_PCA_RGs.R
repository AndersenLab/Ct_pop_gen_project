rm(list=ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(dendextend)
library(ggpubr)
library(cowplot)

# source
source("../../utilities.R")

lineage_colors <- c(
  Tw1 = "#E41A1C",
  Tw2 = "#FF7F00",
  Tw3 = "#FFFF33",
  Tw4 = "#4DAF4A",
  Tw5 = "#40E0D0",
  Tw6 = "#377EB8",
  Tw7 ="#984EA3"
)

GIS_data_raw<-read.table("../../../processed_data/geo_info/Taiwan/Ct_TW_isotype_GIS.txt",
                            sep = '\t',header = TRUE)

GIS_data <- GIS_data_raw

EVarNoRm_ld0.9 <- data.table::fread("../../../processed_data/geo_info/Taiwan/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") %>%
  dplyr::mutate(sum = sum(eigenvalue),
                VarExp = eigenvalue/sum,
                sigEV = ifelse(`p-value` < 0.05, T, F)) %>%
  dplyr::group_by(sigEV) %>%
  dplyr::mutate(sigVarExp = sum(VarExp)) %>%
  dplyr::ungroup()

cdf9no <- data.table::fread("../../../processed_data/geo_info/Taiwan/EIGESTRAT/LD_0.9/NO_REMOVAL/eigenstrat_no_removal.evac", skip = 1)

cdf9_2no <- cdf9no %>%
  dplyr::select(isotype = V1, PC1=V2, PC2=V3, PC3=V4, PC4=V5, PC5=V6, PC6 = V7) 

lineage<-read.csv("../../../processed_data/geo_info/Ct_lineage_all.csv") %>% 
  rename(isotype=sample)

cdf9_7no <- cdf9_2no %>%
  dplyr::full_join(GIS_data,by=c("isotype")) %>%
  dplyr::left_join(lineage, by = c ("isotype"))

cdf9_7no$lineage %>% unique()
pca_TAC_ld0.9_no_rm <-  cdf9_7no
saveRDS(cdf9_7no, file = "../../../processed_data/geo_info/Taiwan/cdf9_7no.rds")

Evout <- EVarNoRm_ld0.9 %>%
  dplyr::select(N, VarExp) %>%
  dplyr::mutate(vperc = VarExp * 100) %>%
  dplyr::slice_head(n = 4)

pc12 <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC1, y = PC2, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 14, color = "black")) +
  labs(x = paste0("PC1(",round(Evout$vperc[1],2),"%)"), y = paste0("PC2(",round(Evout$vperc[2],2),"%)"),
       title = "hard-filtered LD=50-10-0.9\n n=54/54")

pc12_blank <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC1, y = PC2, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21, stroke = 0.15) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 9, color = "black"), legend.position = "none") + 
  labs(x = paste0("PC1(",round(Evout$vperc[1],2),"%)"), y = paste0("PC2(",round(Evout$vperc[2],2),"%)"),
       title = "")

pc34 <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC3, y = PC4, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text = element_text(size = 14, color = "black"), legend.position = "right") + # presentation text
  labs(x = paste0("PC3(",round(Evout$vperc[3],2),"%)"), y = paste0("PC4(",round(Evout$vperc[4],2),"%)"),
       title = "hard-filtered LD=50-10-0.9\n n=56/56, ot=6.0")

pc34_blank <- ggplot(pca_TAC_ld0.9_no_rm, aes(x=PC3, y = PC4, fill = factor(lineage))) +
  geom_point(size = 3, alpha = 1, shape = 21, stroke = 0.15) +
  scale_fill_manual("lineage", values = lineage_colors) +
  theme_bw() +
  theme(text = element_text(size = 9), axis.text = element_text(size = 9, color = "black"), legend.position = "none") + 
  labs(x = paste0("PC3(",round(Evout$vperc[3],2),"%)"), y = paste0("PC4(",round(Evout$vperc[4],2),"%)"),
       title = "")

pc12_blank <- pc12_blank + guides(fill = guide_legend(title = "Groups", ncol = 7))
pc34_blank <- pc34_blank + guides(fill = guide_legend(title = "Groups", ncol = 7))

PCA_main <- ggarrange(
  pc12_blank, pc34_blank,
  labels = c("a", "b"),          
  label.x = 0.01, label.y = 0.95, 
  legend = "bottom",            
  common.legend = TRUE          
)

saveRDS(PCA_main, file = "../../../processed_data/geo_info/Taiwan/Ct_PCA_Taiwan_Lineage.rds")

