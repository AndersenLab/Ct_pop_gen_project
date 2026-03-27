rm(list=ls())

library(ggtree)
library(ggplot2)
library(treeio)
library(ape)
library(dplyr)
library(ggnewscale)
library(phytools)
library(RColorBrewer)
library(cowplot)

source("../utilities.R")

isotype_map_0.9_raw<-treeio::read.tree("../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy.contree")

### set midpoint rooting of phylogenetic tree
isotype_map_0.9<-phytools::midpoint_root(isotype_map_0.9_raw)

## write rooted tree
write.tree(isotype_map_0.9,
           file = "../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy.contree.rooted")

##### annotation_maps
annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)

############ plot the tree with heatmaps #########
plot_tree_heatmaps <- function(tree_file, 
                                    annotation_map_file,
                                    annotation_temperature_file,
                                    tree_layout="circular",
                                    add_tiplab=TRUE,
                                    heatmap_lat_offset = 0.0013,
                                    heatmap_geo_offset = 0.0024
                               ){
  data_tree_file<-ggtree::fortify(tree_file)
  annotation_maps<- annotation_map_file
 
  heatmap<- annotation_maps%>%
    select(lat)
  row.names(heatmap)<-annotation_maps$isotype
  heatmap$lat<-abs(heatmap$lat)
 
  heatmap_3<- annotation_maps%>%
    select(geo)
  row.names(heatmap_3)<-annotation_maps$isotype
 
  # plot tree
  p0<-ggtree::ggtree(data_tree_file, aes(col=geo), layout=tree_layout, size=0.2) %<+% annotation_maps +
    geom_tippoint(aes(color=geo), size=0.1, alpha =0.6, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours)
    )+
    theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  if (add_tiplab==TRUE) {
    p0<- p0+geom_tiplab(aes(label=label),
                        color = "grey50",
                        linetype = "dashed",
                        size = 0.6, hjust = -0.5,
                        align = TRUE, linesize = 0.05,
                        alpha = 0.6)
  }
  p1<- ggtree::open_tree(p0, 1) %>% ggtree::rotate_tree(90)
  p2<-p1+new_scale_fill()+new_scale_color()
  p3<-gheatmap(p2, heatmap, offset = heatmap_lat_offset, width=0.05, font.size=2.5, 
               color = NULL, hjust = -0.1, colnames_level=colnames(heatmap), 
               colnames_angle=0, 
               legend_title=" Absolute latitude",
               custom_column_labels  = rep("", ncol(heatmap))
  ) + 
    scale_fill_gradient2(
      low = "#D73027", mid = "#FFFFBF",high = "#4575B4",
      midpoint = mean(c(max(heatmap$lat, na.rm = TRUE), min(heatmap$lat, na.rm = TRUE))),
      na.value = "white",
      name = "Absolute latitude"
    )
  p6<-p3 + new_scale_fill() + new_scale_color()
  p7<-gheatmap(p6, heatmap_3, offset = heatmap_geo_offset, width=0.05, font.size=2.5, 
               color = NULL,
               hjust = -0.1, colnames_level=colnames(heatmap_3), 
               colnames_angle=0, legend_title="",
               custom_column_labels  = rep("", ncol(heatmap_3))
  ) + 
    scale_fill_manual(values = geo.colours,
                      name = "Geo")+
    theme(legend.position = c(1.1,0.92),
          legend.background = element_rect(fill = "transparent",
                                           colour = NA))
  return(p7)
}

plot_0.9_tree<-plot_tree_heatmaps(tree_file = isotype_map_0.9,
                                  annotation_map_file = annotation_maps)
ggsave("../../figures/raw_FigureS4_Ct_tree.pdf", plot = plot_0.9_tree, width = 7.5, height = 7.5, units = "in")


plot_equal_angle_tree <- function(tree_file, 
                                  annotation_map_file){
  
  annotation_maps<-annotation_map_file
  
  p_ea<-ggtree::ggtree(tree_file, 
                       layout="equal_angle", size=0.15) %<+% annotation_maps +
    geom_tippoint(aes(color=geo), size=1, alpha =0.8, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours)
    )+
    theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  return(p_ea)
}

annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)
plot_0.9_tree_equal_angle<-plot_equal_angle_tree(isotype_map_0.9_raw, 
                                                 annotation_map_file = annotation_maps)
ggsave("../../figures/FigureS5_tree_unrooted.pdf", 
       plot = plot_0.9_tree_equal_angle, 
       width = 7, height = 7, units = "in")




