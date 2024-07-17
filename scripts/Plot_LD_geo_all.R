#plot LD_all

rm(list = ls())

library(dplyr)
library(ggplot2)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")

#imput raw data
LD_all_summary_All_Raw<-read.table("../processed_data/LD_all_isotype/C_tro_all.summary")
LD_all_summary_Africa_Raw<-read.table("../processed_data/LD_geo/C_tro_Africa_LD.ld_LD.summary",sep= ' ')
LD_all_summary_Australia_Raw<-read.table("../processed_data/LD_geo/C_tro_Australia_LD.ld_LD.summary",sep= ' ')
LD_all_summary_Caribbean_Raw<-read.table("../processed_data/LD_geo/C_tro_Caribbean_LD.ld_LD.summary",sep= ' ')
LD_all_summary_Central_America_Raw<-read.table("../processed_data/LD_geo/C_tro_Central America_LD.ld_LD.summary",sep= ' ')
LD_all_summary_Hawaii_Raw<-read.table("../processed_data/LD_geo/C_tro_Hawaii_LD.ld_LD.summary",sep= ' ')
LD_all_summary_South_America_Raw<-read.table("../processed_data/LD_geo/C_tro_South America_LD.ld_LD.summary",sep= ' ')
LD_all_summary_Taiwan_Raw<-read.table("../processed_data/LD_geo/C_tro_Taiwan_LD.ld_LD.summary",sep= ' ')




# rearrange the data set
rearrange_data <- function(summary_data){

# remove the wrong header
  summary_data<-summary_data%>% 
  dplyr::filter(V1!="CHR_A"& V1!="MtDNA")
# reset the header names
colnames(summary_data)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
#change the "23" into "X"
summary_data$chr_A<-gsub("23","X",summary_data$chr_A)
summary_data$chr_B<-gsub("23","X",summary_data$chr_B)


summary_data$Distance<-(summary_data$pos_B-summary_data$pos_A)/10
# summary_data[summary_data$chr_A != summary_data$chr_B, ]$Distance<-"NA" #remove intrachromosomal LD

return(summary_data)
}



LD_all_summary_All<-rearrange_data(LD_all_summary_All_Raw)
LD_all_summary_Africa<-rearrange_data(LD_all_summary_Africa_Raw)
LD_all_summary_Australia<-rearrange_data(LD_all_summary_Australia_Raw)
LD_all_summary_Caribbean<-rearrange_data(LD_all_summary_Caribbean_Raw)
LD_all_summary_Central_America<-rearrange_data(LD_all_summary_Central_America_Raw)
LD_all_summary_Hawaii<-rearrange_data(LD_all_summary_Hawaii_Raw)
LD_all_summary_South_America<-rearrange_data(LD_all_summary_South_America_Raw)
LD_all_summary_Taiwan<-rearrange_data(LD_all_summary_Taiwan_Raw)


LD_all_summary_All$Geo<-c("All")
LD_all_summary_Africa$Geo<-c("Africa")
LD_all_summary_Australia$Geo<-c("Australia")
LD_all_summary_Caribbean$Geo<-c("Caribbean")
LD_all_summary_Central_America$Geo<-c("Central America")
LD_all_summary_Hawaii$Geo<-c("Hawaii")
LD_all_summary_South_America$Geo<-c("South America")
LD_all_summary_Taiwan$Geo<-c("Taiwan")

LD_all_summary<-rbind(LD_all_summary_All,
                      LD_all_summary_Africa,
                      LD_all_summary_Australia,
                      LD_all_summary_Caribbean,
                      LD_all_summary_Central_America,
                      LD_all_summary_Hawaii,
                      LD_all_summary_South_America,
                      LD_all_summary_Taiwan)

LD_all_summary$Geo<-factor(LD_all_summary$Geo,levels = c("All",
                                                         "Africa",
                                                         "Australia",
                                                         "Caribbean",
                                                         "Central America",
                                                         "Hawaii",
                                                         "South America",
                                                         "Taiwan"))
# generate data set for the main figure
LD_main_fig_summary<-LD_all_summary %>% 
  dplyr::filter(Geo != "Africa" & Geo != "Australia")


# plot the figure
plot_figure<-function(LD_summary_data){

  geo.colours <- c("All" ="lightgrey","Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                   "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                   "Taiwan" = "#E5C494")

p1 <- ggplot(LD_summary_data, aes(x = Distance, y = Mean, color= Geo, fill=Geo, ordered = FALSE))  +
  labs(title ="", x = "Distance (Mb)", y = expression(paste("Linkage disequilibrium (",r^{2},")"))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),limits=c(0,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right")+
  # theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(),
        # strip.text = element_text(size = rel(1)),
        panel.grid.minor = element_blank()) +
  stat_summary(fun = median, geom="line",linewidth=0.75)+
  theme(axis.title.y = element_text(size = 9))+ # face = "bold",
  theme(axis.title.x = element_text(size = 9))+ #face = "bold",
  # theme(axis.text.x = element_text(size = 9))+ # face = "bold",
  # theme(axis.text.y = element_text(size = 9))+ # face = "bold",
  # theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) 
  scale_color_manual(values = geo.colours,name="") +
  scale_fill_manual(values = geo.colours,name="")
  # theme(legend.title = c("Chromsom"))
  # facet_grid(species~., space="free")

return(p1)
}





plot_LD_all_summary<-plot_figure(LD_all_summary)
plot_LD_all_summary
# save image 7.5 * 4 inches

# ggsave("../plots/Plot_LD_geo_all.pdf", plot = plot_LD_all_summary, width = 7.5, height = 4, units = "in")


plot_LD_main_fig_summary<-plot_figure(LD_main_fig_summary)
plot_LD_main_fig_summary
ggsave("../plots/Plot_LD_geo_all_main_fig.pdf", plot = plot_LD_main_fig_summary, width = 7.5, height = 4, units = "in")




