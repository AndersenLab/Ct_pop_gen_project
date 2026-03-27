rm(list = ls())

source("../utilities.R")

library(dplyr)
library(ggplot2)
library(RColorBrewer)

LD_all_summary_Ct<-read.table("../../processed_data/LD_three_species/C_tro_all.summary")
LD_all_summary_Cb<-read.table("../../processed_data/LD_three_species/C_bri_all.summary")
LD_all_summary_Ce<-read.table("../../processed_data/LD_three_species/C_ele_all.summary")

rearrange_data <- function(summary_data){
  summary_data<-summary_data%>% 
  dplyr::filter(V1!="CHR_A"& V1!="MtDNA")
colnames(summary_data)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
summary_data$chr_A<-gsub("23","X",summary_data$chr_A)
summary_data$chr_B<-gsub("23","X",summary_data$chr_B)
summary_data$Distance<-(summary_data$pos_B-summary_data$pos_A)/10
return(summary_data)
}

LD_Ct<-rearrange_data(LD_all_summary_Ct)
LD_Cb<-rearrange_data(LD_all_summary_Cb)
LD_Ce<-rearrange_data(LD_all_summary_Ce)

LD_Ct$Species<-c("C. tropicalis")
LD_Cb$Species<-c("C. briggsae")
LD_Ce$Species<-c("C. elegans")

LD_all_summary<-rbind(LD_Ct,
                      LD_Cb,
                      LD_Ce)

LD_all_summary$Species<-factor(LD_all_summary$Species,
                               levels = c("C. elegans",
                                          "C. briggsae",
                                          "C. tropicalis"))

plot_figure<-function(LD_summary_data){
p1 <- ggplot(LD_summary_data, aes(x = Distance, y = Mean, color= Species, fill=Species, ordered = FALSE))  +
  labs(title ="", x = "Distance (Mb)", y = expression(paste("Linkage disequilibrium (",r^{2},")"))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),limits=c(0,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right")+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  stat_summary(fun = median, geom="line",linewidth=0.75)+
  theme(axis.title.y = element_text(size = 9))+
  theme(axis.title.x = element_text(size = 9))+
  scale_color_manual(values = c(
    "C. elegans"    = "#DB6333",
    "C. briggsae"   = "#53886C",
    "C. tropicalis" = "#0719BC"
  )) +
  scale_fill_manual(values = c(
    "C. elegans"    = "#DB6333",
    "C. briggsae"   = "#53886C",
    "C. tropicalis" = "#0719BC"
  )) 

return(p1)
}

plot_LD_all_summary<-plot_figure(LD_all_summary)
ggsave("../../figures/FigureS3_LD_decay.pdf", 
       plot = plot_LD_all_summary, 
       width = 7, height = 4, units = "in")



