# Dendrogram and geographical frequency of Medea haplotypes on chromosome III and V


# BiocManager::install("fastreeR")

rm(list = ls())

library(fastreeR)
library(ggtree)
library(ape)
library(ggplot2)
library(TreeTools)
library(dplyr)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")



####### 1. Chr3 - JU1373 Clade #######

myVcfDist_Chr3 <- fastreeR::vcf2dist(inputFile = "../processed_data/Medea/Medea_Chr3_127_151.recode.vcf", threads = 2)
myVcfTree_Chr3 <- fastreeR::dist2tree(inputDist = myVcfDist_Chr3)
myVcfTree_Chr3 <- gsub('"', '', myVcfTree_Chr3)

# Convert the tree in string format to a phylo object.
phylo_Chr3 <- ape::read.tree(text = myVcfTree_Chr3)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p1<-ggtree::ggtree(phylo_Chr3, aes(col=geo), 
                   # layout="circular",
                   size=0.3) %<+% annotation_maps + # %<+% annotation_maps 

  geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+

  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+

  geom_tiplab(aes(label=NA, col=NA), size=0.4)+

  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p1

p2<-p1+geom_strip("JU3170","JU1641",label = "JU1373\nClade",offset = 0.02, offset.text = 0.1,
                  barsize = 1, extend = 0.2, fontsize = 4, angle = 0, hjust = 0.9,
                  color = "#7B5AA3")+
  theme(legend.position = "none")
p2<-p2+geom_treescale(x = 0.05, y = 60,width = 0.05,offset = 3)+
  ggtitle("Chromosome III Medea region")+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.1, color = "black"))
p2
# ggsave("Chr3.pdf", plot = p2, width = 7.5, height = 7.5, units = "in")
# 
# ggsave("Chr3_small.pdf", plot = p2, width = 3.75, height = 9, units = "in")




######## Chr 3 end ##########





####### 2. Chr5 - JU1373 and NIC58 Clade #######


myVcfDist_Chr5 <- fastreeR::vcf2dist(inputFile = "../processed_data/Medea/Medea_Chr5_1294_1332.recode.vcf", threads = 2)
myVcfTree_Chr5 <- fastreeR::dist2tree(inputDist = myVcfDist_Chr5)
myVcfTree_Chr5 <- gsub('"', '', myVcfTree_Chr5)

# Convert the tree in string format to a phylo object.
phylo_Chr5 <- ape::read.tree(text = myVcfTree_Chr5)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p3<-ggtree::ggtree(phylo_Chr5, aes(col=geo), 
                   # layout="circular",
                   size=0.3) %<+% annotation_maps + # %<+% annotation_maps 

  ggtree::geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+

  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+

  geom_tiplab(aes(label=NA, col=NA), size=0.4)+

  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p3

p4<-p3 + geom_strip("NIC1534","NIC1531",label = "JU1373\nClade",offset = 0.02, offset.text = 0.1,
                  barsize = 1, extend = 0.2, fontsize = 4, angle = 0, hjust = 0.9,
                  color = "#7B5AA3")+
  theme(legend.position = "none")
p5<-p4 +geom_strip("NIC1382","NIC1180",label = "NIC58\nClade",offset = 0.02, offset.text = 0.1,
                   barsize = 1, extend = 0.2, fontsize = 4, angle = 0, hjust = 0.9,
                   color = "#EF3751")
p5
p6<-p5+geom_treescale(x = 0.05, y = 60,width = 0.05,offset = 3)+
  ggtitle("Chromosome V Medea region")+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.1, color = "black"))
p6
# ggsave("Chr5.pdf", plot = p6, width = 7.5, height = 7.5, units = "in")
# 
# ggsave("Chr5_small.pdf", plot = p6, width = 3.75, height = 10, units = "in")


##### Chr5 - end #####




##### 3. merge figures #####

library(ggpubr)
plot_all<-ggpubr::ggarrange(p2, p6, 
                            ncol = 2, 
                            widths = c(1, 1),
                            labels = c("a","b"))
# ggsave("Medea_with_out_pi.pdf", plot = plot_all, width = 7.5, height = 9, units = "in")

#### merge - end ######






#### 4. frequency of alternative haplotypes & add pie Chart#####


library(BAMMtools)
library(phytools)


# function 1. filter_strain_list from specific clade of the tree

filter_strain_list<-function(phylo_file, start, stop){
# Find the most recent common ancestor
tmp_a<-BAMMtools::getmrca(phylo_file, start, stop)
tmp_b<-phytools::getDescendants(phylo_file,tmp_a)
# print the label
stain_list<-phylo_file[["tip.label"]][tmp_b]
stain_list<-stain_list[!is.na(stain_list)]
return(stain_list)
}


# calculate the strain lists

Chr3_JU<- filter_strain_list(phylo_Chr3,"JU1641","JU3170")
Chr5_NIC<- filter_strain_list(phylo_Chr5,"NIC1180","NIC1382")
Chr5_JU<- filter_strain_list(phylo_Chr5,"NIC1531","NIC1534")




# function 2: plot_Pie_chart

plot_pie_chart<-function(strain_list,pie_title=NULL,title_color = "black"){
indep_isotype_info<-read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv")
geo_freq <- indep_isotype_info %>%
  dplyr::filter(isotype %in% strain_list) %>% 
  dplyr::select(isotype,geo) %>% 
  dplyr::group_by(geo) %>%
  dplyr::summarize(frequency = dplyr::n()) %>%
  dplyr::arrange(desc(frequency))


pie_plot<-ggplot(geo_freq, aes(x = "", y = frequency, fill = geo)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  # geom_text(aes(label = frequency), position = position_stack(vjust = 0.5)) +
  geom_text(aes(x = 1.7, label = frequency, color = geo),
            position = position_stack(vjust = 0.5),
            size = 5,
            hjust = 0.5) +
  scale_fill_manual(values = geo.colours)+
  scale_color_manual(values = geo.colours)
# save image 4*4 inches

if (!is.null(title)) {
  pie_plot <- pie_plot +
    ggtitle(pie_title)+
    theme(plot.title = element_text(hjust = 0.5,color = title_color))

  
  return(pie_plot)
  
}



return(pie_plot)

}


pie_Chr3_JU<-plot_pie_chart(Chr3_JU,
                             pie_title = "JU1373 Clade",
                             title_color = "#7B5AA3")
pie_Chr3_JU

pie_Chr5_NIC<-plot_pie_chart(Chr5_NIC,
                             pie_title = "NIC58 Clade",
                             title_color = "#EF3751")
pie_Chr5_NIC

pie_Chr5_JU<-plot_pie_chart(Chr5_JU,
                            pie_title = "JU1373 Clade",
                            title_color = "#7B5AA3")
pie_Chr5_JU











##### 5. merge all plots ####





# Chr5 x range 
layer_scales(p2)$x$get_limits()
# 0.0000000 0.4656064

# Chr5 y range 
layer_scales(p2)$y$get_limits()
# 1 518

# Merge Chr5
p2_1<-p2 + 
  ggplot2::annotation_custom(ggplotGrob(pie_Chr3_JU), 
                             xmin = 0, xmax = 0.3, 
                             ymin = 140, ymax = 360)

p2_1
# ggsave("test_p2_1.pdf", plot = p2_1, width = 3.75, height = 9, units = "in")





# Chr5 x range 
layer_scales(p6)$x$get_limits()
# 0.0000000 0.53391

# Chr5 y range 
layer_scales(p6)$y$get_limits()
# 0.8 518.2

# Merge Chr5
p6_1<-p6 + 
  ggplot2::annotation_custom(ggplotGrob(pie_Chr5_JU), 
                             xmin = 0, xmax = 0.35, 
                             ymin = 140, ymax = 360)+
  ggplot2::annotation_custom(ggplotGrob(pie_Chr5_NIC), 
                             xmin = 0, xmax = 0.35, 
                             ymin = 290, ymax = 510)
  
p6_1
# ggsave("test_p6_1.pdf", plot = p6_1, width = 3.75, height = 9, units = "in")


# merge all
library(ggpubr)
plot_all_1<-ggpubr::ggarrange(p2_1, p6_1, 
                            ncol = 2, 
                            widths = c(1, 1),
                            labels = c("a","b"))
ggsave("../plots/Medea_all_raw.pdf", plot = plot_all_1, width = 7.5, height = 9, units = "in")



#### end #####




# 
# 
# #### supplement plot legend ####
# 
# 
# p_leg<-ggtree::ggtree(phylo_Chr3, aes(col=geo),
#                    # layout="circular",
#                    size=0.3) %<+% annotation_maps + # %<+% annotation_maps # 引入注释文件
#   # tip color and size
#   geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
#   scale_color_manual(values=geo.colours)+
#   # line between tree and outer ring
#   geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
#   # tip color
#   geom_tiplab(aes(label=NA, col=NA), size=0.4)+
#   # annotation
#   ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
#   ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))+
#   geom_strip("JU3170","JU1641",label = "JU1373\nClade",offset = 0.01, offset.text = 0.025,
#                   barsize = 1, extend = 0.2, fontsize = 1.8, angle = 0, hjust = 0.5,
#                   color = "#7B5AA3")
# 
# ggsave("p_leg_raw.pdf", plot = p_leg, width = 10, height = 10, units = "in")
# 
# 



