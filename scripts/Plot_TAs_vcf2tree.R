# Dendrogram and geographical frequency of TA, slow-1 on chromosome III


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



####### 1. Chr3 - slow1 #######

myVcfDist_slow1 <- fastreeR::vcf2dist(inputFile = "../processed_data/TAs/TAs_slow_1.recode.vcf", threads = 2)
myVcfTree_slow1 <- fastreeR::dist2tree(inputDist = myVcfDist_slow1)
myVcfTree_slow1 <- gsub('"', '', myVcfTree_slow1)

# Convert the tree in string format to a phylo object.
phylo_slow1 <- ape::read.tree(text = myVcfTree_slow1)



geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p<-ggtree::ggtree(phylo_slow1, aes(col=geo), 
                   # layout="circular",
                   size=0.3) %<+% annotation_maps + # %<+% annotation_maps 
  
  geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+
  
  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+

  geom_tiplab(aes(label=NA, col=NA), size=0.4)+
  
  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p
# ggsave("slow1.pdf", plot = p, width = 15, height = 15, units = "in")

# 
# p2<-p1+geom_strip("JU3170","JU1641",label = "JU1373\nClade",offset = 0.01, offset.text = 0.025,
#                   barsize = 1, extend = 0.2, fontsize = 1.8, angle = 0, hjust = 0.5,
#                   color = "#E41A1C")+
#   theme(legend.position = "none")
# 
# p2<-p2+geom_treescale(x = 0.05, y = 60,width = 0.05,offset = 3)+
#   ggtitle("Chromosome III Medea region")+
#   theme(plot.title = element_text(hjust = 0.5,color = "black"))
# p2
# ggsave("Chr3.pdf", plot = p2, width = 7.5, height = 7.5, units = "in")
# 
# ggsave("Chr3_small.pdf", plot = p2, width = 3.75, height = 9, units = "in")




######## Chr 3 end ##########






####### 2. NIL_NIC_ChrII #######

myVcfDist_NIL_NIC_ChrII <- fastreeR::vcf2dist(inputFile = "../processed_data/TAs/TAs_NIL_NIC_ChrII.recode.vcf", threads = 2)
myVcfTree_NIL_NIC_ChrII <- fastreeR::dist2tree(inputDist = myVcfDist_NIL_NIC_ChrII)
myVcfTree_NIL_NIC_ChrII <- gsub('"', '', myVcfTree_NIL_NIC_ChrII)

# Convert the tree in string format to a phylo object.
phylo_NIL_NIC_ChrII <- ape::read.tree(text = myVcfTree_NIL_NIC_ChrII)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p1<-ggtree::ggtree(phylo_NIL_NIC_ChrII, aes(col=geo), 
                   # layout="circular",
                   size=0.3) %<+% annotation_maps + # %<+% annotation_maps
  geom_tippoint(aes(color=geo), size=0.1, alpha =0.5)+
  scale_color_manual(values=geo.colours)+
  geom_tiplab(aes(label=label, col=geo), size = 0.22,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
  # geom_tiplab(aes(label=NA, col=NA), size=0.4)+
  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p1
# ggsave("NIL_NIC_ChrII.pdf", plot = p1, width = 10, height = 10, units = "in")

p1_1<-p1+theme(legend.position = c(0.2, 0.8))+
  ggplot2::theme(plot.margin = unit(c(0, 0.1, 0, 0), "in"))  +
  # geom_strip("EG6180","EG6180",label = "EG6180",offset = 0.002, offset.text = 0.007,
  #                               barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
  #                               color = "#E41A1C")+
  geom_strip("NIC203","NIC203",label = "NIC203",offset = 0.002, offset.text = 0.007,
             barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
             color = "#E41A1C")+
  ggtitle("Phylogenetic Tree of NIL: NIC203 Chr. II") +  
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  vjust = -2,
                                  face = "bold"))  

  

p1_1
ggsave("../plots/NIL_NIC_ChrII_plot.pdf", plot = p1_1, width = 7.5, height = 7.5, units = "in")



# 
# p22<-p11+geom_strip("NIC203","NIC203",label = "NIC203\nMaternal",offset = 0.01, offset.text = 0.025,
#                   barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
#                   color = "#E41A1C")+
#   theme(legend.position = c(0.2, 0.8))
# p22
# ggsave("NIL_NIC_ChrII_test.pdf", plot = p22, width = 7.5, height = 9, units = "in")


# ------------
# 
# 
# p22<-p22+geom_treescale(x = 0.05, y = 60,width = 0.05,offset = 3)+
#   ggtitle("Chromosome III Medea region")+
#   theme(plot.title = element_text(hjust = 0.5,color = "black"))
# p22
# ggsave("Chr3.pdf", plot = p22, width = 7.5, height = 7.5, units = "in")
# 
# ggsave("Chr3_small.pdf", plot = p22, width = 3.75, height = 9, units = "in")




######## NIL_NIC_ChrII end ##########







####### 3. TAs_NIL_NIC_ChrIII #######

myVcfDist_NIL_NIC_ChrIII <- fastreeR::vcf2dist(inputFile = "../processed_data/TAs/TAs_NIL_NIC_ChrIII.recode.vcf", threads = 2)
myVcfTree_NIL_NIC_ChrIII <- fastreeR::dist2tree(inputDist = myVcfDist_NIL_NIC_ChrIII)
myVcfTree_NIL_NIC_ChrIII <- gsub('"', '', myVcfTree_NIL_NIC_ChrIII)

# Convert the tree in string format to a phylo object.
phylo_NIL_NIC_ChrIII <- ape::read.tree(text = myVcfTree_NIL_NIC_ChrIII)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p2<-ggtree::ggtree(phylo_NIL_NIC_ChrIII, aes(col=geo), 
                    # layout="circular",
                    size=0.3) %<+% annotation_maps + # %<+% annotation_maps
  geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+
  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
  geom_tiplab(aes(label=NA, col=NA), size=0.4)+
  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p2
# ggsave("NIL_NIC_ChrIII.pdf", plot = p2, width = 10, height = 10, units = "in")


p2_1<-p2+theme(legend.position = c(0.2, 0.8))+
  ggplot2::theme(plot.margin = unit(c(0, 0.1, 0, 0), "in"))  +
  # geom_strip("EG6180","EG6180",label = "EG6180",offset = 0.002, offset.text = 0.007,
  #                               barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
  #                               color = "#E41A1C")+
  geom_strip("NIC203","NIC203",label = "NIC203",offset = 0.01, offset.text = 0.02,
             barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
             color = "#E41A1C")+
  ggtitle("Phylogenetic Tree of NIL: NIC203 Chr. III") +  
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  vjust = -2,
                                  face = "bold"))  



p2_1
ggsave("../plots/NIL_NIC_ChrIII_plot.pdf", plot = p2_1, width = 7.5, height = 7.5, units = "in")




######## NIL_NIC_ChrIII end ##########







####### 4. TAs_NIL_NIC_ChrV #######

myVcfDist_NIL_NIC_ChrV <- fastreeR::vcf2dist(inputFile = "../processed_data/TAs/TAs_NIL_NIC_ChrV.recode.vcf", threads = 2)
myVcfTree_NIL_NIC_ChrV <- fastreeR::dist2tree(inputDist = myVcfDist_NIL_NIC_ChrV)
myVcfTree_NIL_NIC_ChrV <- gsub('"', '', myVcfTree_NIL_NIC_ChrV)

# Convert the tree in string format to a phylo object.
phylo_NIL_NIC_ChrV <- ape::read.tree(text = myVcfTree_NIL_NIC_ChrV)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p3<-ggtree::ggtree(phylo_NIL_NIC_ChrV, aes(col=geo), 
                      # layout="circular",
                      size=0.3) %<+% annotation_maps + # %<+% annotation_maps
  geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+
  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
  geom_tiplab(aes(label=NA, col=NA), size=0.4)+
  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p3
# ggsave("NIL_NIC_ChrV.pdf", plot = p3, width = 10, height = 10, units = "in")


p3_1<-p3+theme(legend.position = c(0.2, 0.8))+
  ggplot2::theme(plot.margin = unit(c(0, 0.1, 0, 0), "in"))  +
  # geom_strip("EG6180","EG6180",label = "EG6180",offset = 0.002, offset.text = 0.007,
  #                               barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
  #                               color = "#E41A1C")+
  geom_strip("NIC203","NIC203",label = "NIC203",offset = 0.01, offset.text = 0.03,
             barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
             color = "#E41A1C")+
  ggtitle("Phylogenetic Tree of NIL: NIC203 Chr. V") +  
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  vjust = -2,
                                  face = "bold"))  



p3_1
ggsave("../plots/NIL_NIC_ChrV_plot.pdf", plot = p3_1, width = 7.5, height = 7.5, units = "in")



######## NIL_NIC_ChrV end ##########











####### 5. TAs_NIL_EG_ChrII #######

myVcfDist_NIL_EG_ChrII <- fastreeR::vcf2dist(inputFile = "../processed_data/TAs/TAs_NIL_EG_ChrII.recode.vcf", threads = 2)
myVcfTree_NIL_EG_ChrII <- fastreeR::dist2tree(inputDist = myVcfDist_NIL_EG_ChrII)
myVcfTree_NIL_EG_ChrII <- gsub('"', '', myVcfTree_NIL_EG_ChrII)

# Convert the tree in string format to a phylo object.
phylo_NIL_EG_ChrII <- ape::read.tree(text = myVcfTree_NIL_EG_ChrII)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p4<-ggtree::ggtree(phylo_NIL_EG_ChrII, aes(col=geo), 
                      # layout="circular",
                      size=0.3) %<+% annotation_maps + # %<+% annotation_maps
  geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+
  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
  geom_tiplab(aes(label=NA, col=NA), size=0.4)+
  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p4
# ggsave("NIL_EG_ChrII.pdf", plot = p4, width = 10, height = 10, units = "in")

p4_1<-p4+theme(legend.position = c(0.2, 0.8))+
  ggplot2::theme(plot.margin = unit(c(0, 0.1, 0, 0), "in"))  +

  geom_strip("EG6180","EG6180",label = "EG6180",offset = 0.01, offset.text = 0.03,
             barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
             color = "#E41A1C")+
  ggtitle("Phylogenetic Tree of NIL: EG6180 Chr. II") +  
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  vjust = -2,
                                  face = "bold"))  



p4_1
ggsave("../plots/NIL_EG_ChrII_plot.pdf", plot = p4_1, width = 7.5, height = 7.5, units = "in")



######## NIL_EG_ChrII end ##########












####### 6. TAs_NIL_EG_ChrV #######

myVcfDist_NIL_EG_ChrV <- fastreeR::vcf2dist(inputFile = "../processed_data/TAs/TAs_NIL_EG_ChrV.recode.vcf", threads = 2)
myVcfTree_NIL_EG_ChrV <- fastreeR::dist2tree(inputDist = myVcfDist_NIL_EG_ChrV)
myVcfTree_NIL_EG_ChrV <- gsub('"', '', myVcfTree_NIL_EG_ChrV)

# Convert the tree in string format to a phylo object.
phylo_NIL_EG_ChrV <- ape::read.tree(text = myVcfTree_NIL_EG_ChrV)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

annotation_maps<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE)

p5<-ggtree::ggtree(phylo_NIL_EG_ChrV, aes(col=geo), 
                      # layout="circular",
                      size=0.3) %<+% annotation_maps + # %<+% annotation_maps
  geom_tippoint(aes(color=geo), size=0.4, alpha =0.8)+
  scale_color_manual(values=geo.colours)+
  geom_tiplab(aes(label=label, col=geo), size = 0.52,hjust=-0.5, align=TRUE, linesize=0.5, alpha=0.8)+
  geom_tiplab(aes(label=NA, col=NA), size=0.4)+
  ggplot2::theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) +
  ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "in"))

p5
# ggsave("../plots/NIL_EG_ChrV.pdf", plot = p5, width = 10, height = 10, units = "in")




p5_1<-p5+theme(legend.position = c(0.2, 0.8))+
  ggplot2::theme(plot.margin = unit(c(0, 0.1, 0, 0), "in"))  +
  
  geom_strip("EG6180","EG6180",label = "EG6180",offset = 0.01, offset.text = 0.03,
             barsize = 2, extend = 0.2, fontsize = 3, angle = 0, hjust = 0.5,
             color = "#E41A1C")+
  ggtitle("Phylogenetic Tree of NIL: EG6180 Chr. V") +  
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14, 
                                  vjust = -2,
                                  face = "bold"))  



p5_1
ggsave("../plots/NIL_EG_ChrV_plot.pdf", plot = p5_1, width = 7.5, height = 7.5, units = "in")


######## NIL_EG_ChrV end ##########



