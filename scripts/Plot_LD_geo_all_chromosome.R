#plot LD_all

rm(list = ls())

library(dplyr)
library(ggplot2)
library(grid)
library(ggpubr)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,"Set1")
brewer.pal(9,"Set1")

#input raw data
LD_all_summary_All_Raw<- read.table("../processed_data/LD_all_isotype/C_tro_all.summary",sep = ' ')
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
  dplyr::filter(V1!="CHR_A" & V1!="MtDNA")
# reset the header names
colnames(summary_data)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
#change the "23" into "X"
summary_data$chr_A<-gsub("23","X",summary_data$chr_A)
summary_data$chr_B<-gsub("23","X",summary_data$chr_B)


summary_data$Distance<-(summary_data$pos_B-summary_data$pos_A)/10
# summary_data[summary_data$chr_A != summary_data$chr_B, ]$Distance<-"NA" #remove intrachromosomal LD

return(summary_data)
}










# plot the figure
plot_figure<-function(LD_summary_data,y_range = NULL){

chr_colors<-c("I"="#E41A1C",
              "II"="#FF7F00",
              "III"="#FFFF33",
              "IV"="#4DAF4A",
              "V"="#377EB8",
              "X"="#984EA3")

p1 <- ggplot2::ggplot(LD_summary_data, aes(x = Distance, y = Mean, color= chr_A, fill=chr_A, ordered = FALSE))  +
  labs(title ="", x = "Distance (Mb)", y = expression(paste("Linkage Disequilibrium (",r^{2},")"))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),limits=c(0,1)) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0, face = "bold"),strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  theme(legend.position = c(0.85,0.7))+
  # theme(strip.text.y = element_text(face = "italic"))+
  theme(strip.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = rel(1))) +
  stat_summary(fun = median, geom="line",linewidth=0.75)+
  # theme(plot.margin = unit(c(1, 3, 1.5, 3), "pt")) 
  scale_color_manual(values = chr_colors,name="") +
  scale_fill_manual(values = chr_colors,name="")
  # theme(legend.title = c("Chromsom"))
  # facet_grid(species~., space="free")

# custom y axis range
if (!is.null(y_range)) {
  p1 <- p1 + ylim(y_range)
}

return(p1)
}




LD_all_summary_All<-rearrange_data(LD_all_summary_All_Raw)
LD_all_summary_Africa<-rearrange_data(LD_all_summary_Africa_Raw)
LD_all_summary_Australia<-rearrange_data(LD_all_summary_Australia_Raw)
LD_all_summary_Caribbean<-rearrange_data(LD_all_summary_Caribbean_Raw)
LD_all_summary_Central_America<-rearrange_data(LD_all_summary_Central_America_Raw)
LD_all_summary_Hawaii<-rearrange_data(LD_all_summary_Hawaii_Raw)
LD_all_summary_South_America<-rearrange_data(LD_all_summary_South_America_Raw)
LD_all_summary_Taiwan<-rearrange_data(LD_all_summary_Taiwan_Raw)


plot_All<-plot_figure(LD_all_summary_All)
plot_Africa<-plot_figure(LD_all_summary_Africa)
plot_Australia<-plot_figure(LD_all_summary_Australia)
plot_Caribbean<-plot_figure(LD_all_summary_Caribbean)
plot_Central_America<-plot_figure(LD_all_summary_Central_America)
plot_Hawaii<-plot_figure(LD_all_summary_Hawaii)
plot_South_America<-plot_figure(LD_all_summary_South_America)
plot_Taiwan<-plot_figure(LD_all_summary_Taiwan)




plot_All
plot_Africa
plot_Australia
plot_Caribbean
plot_Central_America
plot_Hawaii
plot_South_America
plot_Taiwan





# Assemble the 8 plots
p <- ggpubr::ggarrange(plot_All+ rremove("ylab") + rremove("xlab"),
                       plot_Africa+ rremove("ylab") + rremove("xlab"),
                       plot_Australia+ rremove("ylab") + rremove("xlab"),
                       plot_Caribbean+ rremove("ylab") + rremove("xlab"),
                       plot_Central_America+ rremove("ylab") + rremove("xlab"),
                       plot_Hawaii+ rremove("ylab") + rremove("xlab"), 
                       plot_South_America+ rremove("ylab") + rremove("xlab"),
                       plot_Taiwan+ rremove("ylab") + rremove("xlab"),
                       labels = c("All","Africa", "Australia", "Caribbean", "Central America", "Hawaii", "South America", "Taiwan"),
                       vjust = 3,
                       ncol = 2, nrow = 4,
                       align = "hv",
                       common.legend = TRUE, legend = "top",
                       font.label = list(size = 9, color = "black", face = "bold", family = NULL, position = "top"))
p1<-annotate_figure(p,left = grid::textGrob(expression(paste("Linkage disequilibrium (",r^{2},")")), rot = 90, vjust = 0.5, gp = gpar(cex = 1)),
                bottom = textGrob("Distance (Mb)", gp = gpar(cex = 1)))

p1

ggsave("../plots/Plot_LD_geo_all_chromosome.pdf", plot = p1, width = 7.5, height = 9, units = "in")






