rm(list=ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(dendextend)
library(ggpubr)
library(cowplot)
library(patchwork)
library(ggridges)

source("../utilities.R")

geo_info_raw<-read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw 

tracy_for_plot <- data.table::fread("../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/TracyWidom_statistics_no_removal.tsv") %>%
  dplyr::mutate(sum = sum(eigenvalue),
                VarExp = eigenvalue/sum,
                sigEV = ifelse(`p-value` < 0.05, T, F)) %>%
  dplyr::group_by(sigEV) %>%
  dplyr::mutate(sigVarExp = sum(VarExp)) %>%
  dplyr::ungroup()

cdf9no <- data.table::fread("../../processed_data/Ct_pruned_VCF_and_PCA/EIGESTRAT/LD_0.9/NO_REMOVAL/eigenstrat_no_removal.evac", skip = 1)

cdf9_2no <- cdf9no %>%
  dplyr::select(isotype = V1, PC1=V2, PC2=V3, PC3=V4, PC4=V5, PC5=V6, PC6 = V7)
cdf9_3no <- cdf9_2no %>% dplyr::select(-isotype)

cdf9_7no <- cdf9_2no %>%
  dplyr::left_join(geo_info,by=c("isotype")) 

pca_TAC_ld0.9_no_rm <- cdf9_7no 

###### Plot PCA function ######
plot_PCA<-function(PCA_input, tracy_for_plot_input, x_axis, y_axis){
  
  label_x<-tracy_for_plot_input %>% 
    dplyr::filter(N == as.numeric(sub("PC", "", x_axis))) %>% 
    dplyr::select(VarExp) %>%
    pull(VarExp) %>%
    `*`(100) %>%
    round(digits = 2) %>%
    as.character()
  
  label_y<-tracy_for_plot_input %>% 
    dplyr::filter(N == as.numeric(sub("PC", "", y_axis))) %>% 
    dplyr::select(VarExp) %>%
    pull(VarExp) %>%
    `*`(100) %>%
    round(digits = 2) %>%
    as.character()
  
  p1_1<-ggplot2::ggplot(PCA_input)+
    geom_point(shape=16, alpha=0.8, size=1.5, aes(x=.data[[x_axis]],y=.data[[y_axis]], color=geo))+
    scale_color_manual(values = geo.colours, name = "geo") +
    theme_bw() +
    labs(x=paste(x_axis," (",label_x,"%)",sep = ""),
         y=paste(y_axis," (",label_y,"%)",sep = ""))+
    theme(axis.title = element_text(size=8, face = "bold",color = "black"),
          axis.text = element_text(size=7, color = "black"),
          legend.position='none',
          panel.grid = element_blank())
  return(p1_1)
  
}

p_PC1_PC2<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC1",y_axis="PC2")

p_PC3_PC4<- plot_PCA(PCA_input=pca_TAC_ld0.9_no_rm,
                     tracy_for_plot_input = tracy_for_plot,
                     x_axis="PC3",y_axis="PC4")

#### output PCA table ###
Table_PCA<-pca_TAC_ld0.9_no_rm %>% 
  dplyr::filter(PC1 < -0.05) %>% 
  dplyr::select(-strain, -lat, -long, -PC5, -PC6)

write.table(Table_PCA,
            "../../tables/TableS3_PCA.tsv",
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

Ct_all_isotypes_PCA_plot_PC1_4 <- (
  p_PC1_PC2 | p_PC3_PC4  
)

saveRDS(Ct_all_isotypes_PCA_plot_PC1_4, 
        file = "../../processed_data/Ct_pruned_VCF_and_PCA/Ct_all_isotypes_PCA_plot_PC1_4.rds")

####################################
########## some details ###########
####################################
big_cluster<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC1 > -0.1) 
nrow(big_cluster)
#584

left_up_cluster<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC1 < -0.1) %>% 
  filter(PC2 > 0) 
nrow(left_up_cluster)
# 24

left_down_cluster<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC1 < -0.1) %>% 
  filter(PC2 < 0) 
nrow(left_down_cluster)
# 14
# 1 Taiwan, 13 Hawaii

left_down_PC34<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC4 < -0.1) 
nrow(left_down_PC34)
# 10

left_up_PC34<-pca_TAC_ld0.9_no_rm %>% 
  filter(PC4 > 0.1) 
nrow(left_up_PC34)
# 9

