# The input .treefile were generated with iqtree
# Annotation file "indep_isotype_info_geo.csv" is from Isotype_tropicalis_global_map.R 

rm(list=ls())

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ggnewscale)
library(ggpubr)
library(grid)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(11,"RdYlBu")
brewer.pal(11,"RdYlBu")
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")

isotype_map_0.6<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.6_singletons_eiganstrat_input.min4.phy.treefile")
isotype_map_0.7<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.7_singletons_eiganstrat_input.min4.phy.treefile")
isotype_map_0.8<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.8_singletons_eiganstrat_input.min4.phy.treefile")
isotype_map_0.9<-treeio::read.tree("../processed_data/LD_pruned_trees/LD_0.9_singletons_eiganstrat_input.min4.phy.treefile")


annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)
heatmap<- annotation_maps%>%
  select(lat)
row.names(heatmap)<-annotation_maps$isotype
heatmap$lat<-abs(heatmap$lat)

geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")



plot_tree <- function(tree_file,tree_title){
  
p<-ggtree::ggtree(tree_file, layout="equal_angle") %<+% annotation_maps +
  geom_tippoint(aes(color = geo), size = 0.5,alpha =0.7)+
  scale_color_manual(values=geo.colours) +
  labs(title = tree_title)+
  theme(plot.title = element_text(size=10))
p

return(p)
}




plot_0.6_tree<-plot_tree(isotype_map_0.6,expression(paste("LD pruned (",r^{2}," threshold = 0.6)")))
# ggsave("Tree_0.6.pdf", plot = plot_0.6_tree, width = 7.5, height = 7.5, units = "in")

plot_0.7_tree<-plot_tree(isotype_map_0.7,expression(paste("LD pruned (",r^{2}," threshold = 0.7)")))
# ggsave("Tree_0.7.pdf", plot = plot_0.7_tree, width = 7.5, height = 7.5, units = "in")

plot_0.8_tree<-plot_tree(isotype_map_0.8,expression(paste("LD pruned (",r^{2}," threshold = 0.8)")))
# ggsave("Tree_0.8.pdf", plot = plot_0.8_tree, width = 7.5, height = 7.5, units = "in")

plot_0.9_tree<-plot_tree(isotype_map_0.9,expression(paste("LD pruned (",r^{2}," threshold = 0.9)")))
# ggsave("Tree_0.9.pdf", plot = plot_0.9_tree, width = 7.5, height = 7.5, units = "in")

plot_0.6_tree
plot_0.7_tree
plot_0.8_tree
plot_0.9_tree


p1 <- ggpubr::ggarrange(plot_0.6_tree,
                       plot_0.7_tree,
                       plot_0.8_tree,
                       plot_0.9_tree,
                       labels = c("a","b", "c", "d"),
                       vjust = 3,
                       ncol = 2, nrow = 2,
                       align = "hv",
                       common.legend = TRUE, legend = "top",
                       font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
p1


ggsave("../plots/unrooted_tree_equal_angle.pdf", plot = p1, width = 7.5, height = 7.5, units = "in")



#save image 7.5 inches*inches
#rearrange legend positions and titles with Adobe Illustrator


