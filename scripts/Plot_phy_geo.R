# The input .treefile were generated with iqtree

rm(list=ls())


library(ggplot2)
library(treeio)
library(dplyr)
library(reshape2)
library(stringr)
library(ape)
library(geosphere)
library(ggpubr)
library(grid)



# build a plot function
plot_phy_geo_cor<- function(tree_file,plot_title=NULL,x_lab=NULL,y_lab=NULL)
  {


# Calculate the matrix.
## cophenetic.phylo computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
phy_matrix<-ape::cophenetic.phylo(tree_file)

# sort according to names
phy_row_names <- rownames(phy_matrix)
phy_col_names <- colnames(phy_matrix)
phy_matrix_sorted <- phy_matrix[order(phy_row_names), order(phy_col_names)]


phy_matrix_sorted[upper.tri(phy_matrix_sorted,diag=TRUE)]<-NA
phy_distance<-phy_matrix_sorted %>% 
  reshape2::melt() %>% 
  dplyr::filter(value!="NA") %>% 
  dplyr::mutate(paired = stringr::str_c(Var1,Var2,sep = "_") ) %>%
  dplyr::select(paired,value) %>%
  # dplyr::mutate(value=value/100)%>% # change the unite into kilometer (km)
  dplyr::rename(phy_distance=value)




# Load geo data
geo_info<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE,row.names = 1)
colnames(geo_info)
geo_info<-geo_info %>%
  dplyr::select(long,lat)

# Calculate the matrix. Default unit in the matrix is meter(m)
site_dis_matrix <- geosphere::distm(geo_info, fun = distVincentyEllipsoid)
site_dis_matrix[upper.tri(site_dis_matrix,diag=TRUE)]<-NA
row.names(site_dis_matrix)<-row.names(geo_info)
colnames(site_dis_matrix)<-row.names(geo_info)

# Melt the matrix 
geo_distance<-site_dis_matrix %>% 
  reshape2::melt() %>% 
  dplyr::filter(value!="NA") %>% 
  dplyr::mutate(paired = stringr::str_c(Var1,Var2,sep = "_") ) %>%
  dplyr::select(paired,value) %>%
  dplyr::mutate(value=value/1000)%>% # change the unite into kilometer (km)
  dplyr::rename(geo_distance=value)





merged_data<-merge(phy_distance,geo_distance,by="paired")






p<-ggplot2::ggplot(merged_data, aes(x = geo_distance, y = phy_distance)) +
  geom_point(size = 0.5, alpha = 0.2, color = "lightgrey",stroke=0) +  # 通过颜色表示点的密度
  geom_smooth(method = "lm", se = FALSE, color="red",alpha = 0.7) +  
  # labs(x = "Geographic distance (km)", y = "Phylogenetic distance") +
  theme_classic()+
  ggpubr::stat_cor(method = "pearson")
p

if (!is.null(plot_title)) {
  p <- p + labs(title = plot_title)
}

if (!is.null(x_lab)) {
  p <- p + labs(x = x_lab)
} else{p<-p + theme(axis.title.x = element_blank())}

if (!is.null(y_lab)) {
  p <- p + labs(y = y_lab)
} else{p<-p + theme(axis.title.y = element_blank())}


print(p)
}



# load tree data
tree_LD_0.6<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.6_singletons_eiganstrat_input.min4.phy.treefile")
tree_LD_0.7<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.7_singletons_eiganstrat_input.min4.phy.treefile")
tree_LD_0.8<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.8_singletons_eiganstrat_input.min4.phy.treefile")
tree_LD_0.9<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.9_singletons_eiganstrat_input.min4.phy.treefile")


#plot
p_06<-plot_phy_geo_cor(tree_LD_0.6,plot_title=expression(paste("LD pruned (",r^{2}," threshold = 0.6)")))
p_07<-plot_phy_geo_cor(tree_LD_0.7,plot_title=expression(paste("LD pruned (",r^{2}," threshold = 0.7)")))
p_08<-plot_phy_geo_cor(tree_LD_0.8,plot_title=expression(paste("LD pruned (",r^{2}," threshold = 0.8)")))
p_09<-plot_phy_geo_cor(tree_LD_0.9,plot_title=expression(paste("LD pruned (",r^{2}," threshold = 0.9)")))

p_09_main<-plot_phy_geo_cor(tree_LD_0.9,
                            plot_title="Geo. vs. Phylo. Distance Correlation",
                            x_lab="Geographic distance (km)",
                            y_lab = "Phylogenetic distance")
p_09_main
# ggsave("../plots/phy_geo_LD_0.9.pdf", plot = p_09_main, width = 3.75, height = 3, units = "in")




p_all<-ggpubr::ggarrange(p_06,p_07,p_08,p_09,
                         ncol=2,
                         nrow = 2,
                         common.legend = TRUE)
p_all<-ggpubr::annotate_figure(p_all,
                               left = grid::textGrob(("Phylogenetic distance"), rot = 90, vjust = 0.5),
                               bottom = textGrob("Geographic distance (km)"))
p_all

ggsave("../plots/phy_geo_all.pdf", plot = p_all, width = 7.5, height = 6, units = "in")




