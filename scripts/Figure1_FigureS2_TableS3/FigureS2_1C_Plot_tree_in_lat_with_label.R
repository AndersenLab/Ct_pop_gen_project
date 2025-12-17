# The input .treefile were generated with iqtree
# Annotation file "indep_isotype_info_geo.csv" is from Isotype_tropicalis_global_map.R 

rm(list=ls())

library(ggtree)
library(ggplot2)
library(treeio)
library(ape)
library(dplyr)
library(ggnewscale)
library(phytools)

library(RColorBrewer)

source("../utilities.R")



display.brewer.all()
display.brewer.pal(11,"RdYlBu")
brewer.pal(11,"RdYlBu")
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")

isotype_map_0.9_raw<-treeio::read.tree("../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy.contree")

### set midpoint rooting of phylogenetic tree

isotype_map_0.9<-phytools::midpoint_root(isotype_map_0.9_raw)

## write rooted tree
write.tree(isotype_map_0.9,
           file = "../../processed_data/LD_pruned_trees/LD_0.9/phy_file_LD_0.9.phy.contree.rooted")



##### annotation_maps
annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)





# 
# # function: plot the tree
# plot_tree <- function(tree_file){
#   
# data_tree_file<-ggtree::fortify(tree_file)
# 
# annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)
# heatmap<- annotation_maps%>%
#   select(lat)
# row.names(heatmap)<-annotation_maps$isotype
# heatmap$lat<-abs(heatmap$lat)
# 
# 
# ## geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
# ##                  "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
# ##                  "Taiwan" = "#E5C494")
# 
# 
# # plot tree
# p1<-ggtree::ggtree(data_tree_file, aes(col=geo), layout="circular", size=0.3) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
#   geom_tippoint(aes(color=geo), size=0.1, alpha =0.6)+
#   scale_color_manual(values=geo.colours)+
#   geom_tiplab(aes(label=label, col=geo), size = 0.8,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
#   geom_tiplab(aes(label=NA, col=NA), size=0.8)+
#   theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
#  
# 
# p1<- ggtree::open_tree(p1, 1) %>% ggtree::rotate_tree(90) 
# 
# 
# lat_color<-brewer.pal(11,"RdYlBu")[2:10]
# p2<-p1+new_scale_fill()+new_scale_color()
# p3<-gheatmap(p2, heatmap, offset = 0.0015, width=0.1, font.size=6, 
#              color = NULL, hjust = -0.1, colnames_level=colnames(heatmap), 
#              colnames_angle=0, legend_title=" Absolute latitude") + 
#   scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF",high = "#4575B4",midpoint = mean(c(max(heatmap$lat), min(heatmap$lat))))+
#     #### scale_fill_viridis_c(option = "H")+
#   theme(legend.position = c(0.1,0.1),
#         legend.background = element_rect(fill = "transparent",  
#                                          colour = NA))
# 
# 
# 
# p3
# 
# return(p3)
# }
# 








########### new function ###########
############ plot the tree with  heatmaps #########

# function: plot the tree
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
                       breaks = names(geo.colours) # remove NA from Geo legend
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
  
  
  # p1<- ggtree::open_tree(p0, 90) %>% ggtree::rotate_tree(90) 
  p1<- ggtree::open_tree(p0, 1) %>% ggtree::rotate_tree(90)
  
  
  # lat_color<-brewer.pal(11,"RdYlBu")[2:10]
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
      na.value = "white",  # don't show NA in the heatmap
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
plot_0.9_tree

ggsave("../../plots/FigureS2_raw_Ct_tree_0.9_unrooted.pdf", plot = plot_0.9_tree, width = 7.5, height = 7.5, units = "in")









##########  function: plot the equal_angle tree
plot_equal_angle_tree <- function(tree_file, 
                                  annotation_map_file){
  
  annotation_maps<-annotation_map_file
  
  p_ea<-ggtree::ggtree(tree_file, 
                       # aes(col=geo), 
                       layout="equal_angle", size=0.15) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
    geom_tippoint(aes(color=geo), size=1, alpha =0.5, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
    # geom_tiplab(aes(label=label, col=geo), size = 0.5, hjust=-0.5, align=TRUE, linesize=0.1, alpha=0.6)+
    # geom_tiplab(aes(label=NA, col=NA), size=0.8)+
    theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  return(p_ea)
}

annotation_maps<- read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv", header = TRUE)
plot_0.9_tree_equal_angle<-plot_equal_angle_tree(isotype_map_0.9_raw, 
                                                 annotation_map_file = annotation_maps)
plot_0.9_tree_equal_angle
# ggsave("Ct_tree_0.9_unrooted.pdf", plot = plot_0.9_tree_equal_angle, width = 7.5, height = 7.5, units = "in")











##########  function: plot the equal_angle tree for Figure 1 - no legend ###
library(cowplot)
plot_equal_angle_tree_rotated <- function(tree_file, 
                                          annotation_map_file){
  
  annotation_maps<-annotation_map_file
  
  p_ea<-ggtree::ggtree(tree_file, 
                       # aes(col=geo), 
                       layout="equal_angle", size=0.15) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
    geom_tippoint(aes(color=geo), size=1, alpha =0.5, shape = 16)+
    scale_color_manual(values = geo.colours,
                       name = "Geo",
                       breaks = names(geo.colours) # remove NA from Geo legend
    )+
    # theme_tree2() +
    # geom_tiplab(aes(label=label, col=geo), size = 0.5, hjust=-0.5, align=TRUE, linesize=0.1, alpha=0.6)+
    # geom_tiplab(aes(label=NA, col=NA), size=0.8)+
    theme(legend.title=element_text(face="bold"), 
          legend.position="none",
          legend.box="horizontal", 
          legend.text=element_text(size=rel(0.7))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
  
  
  # convert ggplot abject into grob
  g <- ggplotGrob(p_ea)
  rotated_grob <- grid::grobTree(
    g,
    vp = grid::viewport(
      angle = 130, 
      x = unit(0.5, "npc"), 
      y = unit(0.5, "npc"),
      width = unit(1, "npc"),
      height = unit(1, "npc")
    )
  )
  
  # use cowplot to draw
  p_rotated <- cowplot::ggdraw() + draw_grob(rotated_grob)
  
  return(p_rotated)
  
  
  # return(p_ea)
}


plot_0.9_tree_equal_angle_rotated<-plot_equal_angle_tree_rotated(isotype_map_0.9_raw, 
                                                                 annotation_map_file = annotation_maps)

plot_0.9_tree_equal_angle_rotated
saveRDS(plot_0.9_tree_equal_angle_rotated, 
        file = "../../processed_data/assemble_figure_1/Ct_plot_0.9_tree_equal_angle_rotated.rds")

# ggsave("Ct_tree_0.9_unrooted_rotated.pdf", plot = plot_0.9_tree_equal_angle_rotated, width = 7.5, height = 7.5, units = "in")
#  









