# raw data from https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint/data-download#close
# https://www.nature.com/articles/s41597-022-01284-8#Sec12


### This is the downstream analysis of the R code human_footprint_index.R
### /home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/scripts
# install.packages("ggbeeswarm")

rm(list=ls())


library(dplyr)


#### read C.tropicalis data
# read hfi data
tro_hfi_mean <- read.csv("../processed_data/tif_results/C_tro_samples_with_hfi_mean.csv")
# read geo
tro_geo_data_raw <- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv")
# merge
tro_merged_data <- merge(tro_geo_data_raw, tro_hfi_mean, by.x = "isotype", by.y = "SampleID")
tro_merged_data<-tro_merged_data %>% 
  dplyr::select(isotype,Longitude,Latitude,geo,HFI) %>% 
  dplyr::mutate(species="C. tropicalis")



#### read C.elegans data
# read hfi data
ele_hfi_mean <- read.csv("../processed_data/tif_results/C_ele_samples_with_hfi_mean.csv")
# read 
ele_geo_data_raw <- read.csv("../processed_data/Geo_info/20231213_c_elegans_strain_data.csv")
ele_geo_data<-ele_geo_data_raw %>% 
  dplyr::filter(strain==isotype) %>% 
  dplyr::select(isotype,longitude,latitude) %>% 
  na.omit()
# merge
ele_merged_data <- merge(ele_geo_data, ele_hfi_mean, by.x = "isotype", by.y = "SampleID")
ele_merged_data<-ele_merged_data %>% 
  dplyr::select(isotype,Longitude,Latitude,HFI) %>% 
  dplyr::mutate(species="C. elegans")
ele_merged_data<-ele_merged_data %>% 
  dplyr::mutate(geo = ifelse(Latitude > 18.450477 & Latitude < 22.727595 &
                        Longitude > -161.745723 & Longitude < -154.033320,
                      "Hawaii", "non_Hawaii"))

merged_all<-rbind(tro_merged_data, ele_merged_data)
merged_all<-na.omit(merged_all)
merged_all<-merged_all %>% 
  dplyr::mutate(HFI_type = ifelse(merged_all$HFI < 1, "wilderness", 
                           ifelse(merged_all$HFI < 4, "intact", "highly_modified")))
merged_all<-merged_all %>% 
  dplyr::mutate(geo_type = ifelse(merged_all$geo == "Hawaii", "Hawaii", "Others"))




##### plot figures
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,"Set3")
brewer.pal(12,"Set3")

library(ggplot2)
library(ggpubr)
# library(gghalves)
library(ggbeeswarm)


geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")

geo_type.colours <- c("Hawaii"="#66C2A5","Others" = "darkgrey")



#### 1. plot comparison ele vs tro
compaired_ele_tro <- list(c('C. elegans', 'C. tropicalis'))
p1_1<- ggplot2::ggplot(merged_all, aes(x=species, y=HFI)) + 
  ggbeeswarm::geom_beeswarm(aes(fill = geo_type), 
              shape=21, 
              size = 0.8, 
              cex=0.8,
              stroke = 0.01,
              alpha = 0.8)+ 
  ggplot2::stat_boxplot(geom = 'errorbar', width = 0.2, 
               linewidth=0.2,
               color="black",) + # error bar
  
  ggplot2::geom_boxplot(width = 0.25,
  size=0.3, alpha=0.1,outlier.shape = NA) + # outliers

ggpubr::stat_compare_means(comparisons = compaired_ele_tro, label = "p.signif", method = "wilcox.test", vjust = 0.7) + 
  ggplot2::theme_bw() + 
  ggplot2::theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  ggplot2::scale_fill_manual(values=geo_type.colours)+
  
  ggplot2::ylab("Human Footprint Index (HFI)")+xlab(NULL)+ 
  # ylim(0,50)+
  ggplot2::theme(legend.position="none") +
  ggplot2::labs(title = expression(~italic("C. elegans ")~ "vs." ~italic(" C. tropicalis")~"")) +
  ggplot2::theme(axis.text.x = element_text(face = "italic"))+
  ggplot2::theme(plot.title = element_text(hjust = 0.5,vjust = -1.5))

p1_1

# ggsave("HFI_species.pdf", plot = p1_1, width = 3.75, height = 3, units = "in")






#### 2. plot C.tropicalis geo regions
merged_all_plot_p1_2<-merged_all %>% 
  filter(species=="C. tropicalis")
sorted_geo <- merged_all_plot_p1_2 %>%
  dplyr::group_by(geo) %>%
  dplyr::summarize(mean_HFI = median(HFI)) %>%
  dplyr::arrange(-mean_HFI) %>%
  dplyr::pull(geo)

merged_all_plot_p1_2$geo <- factor(merged_all_plot_p1_2$geo, levels = rev(sorted_geo))


# plot
p1_2_raw<- ggplot(merged_all_plot_p1_2, aes(x=geo, y=HFI, color=geo,fill=geo)) + 
  ggbeeswarm::geom_beeswarm(aes(fill = geo), 
                            shape=21, 
                            cex=0.5,
                            size = 0.8, 
                            stroke = 0.01,
                            alpha = 0.8)+ 
  # ggplot2::geom_violin(alpha=0.5,linewidth = 0.3,width = 1)+
  ggplot2::stat_boxplot(geom = 'errorbar', width = 0.2, linewidth=0.2) +
  ggplot2::geom_boxplot(width = 0.25,

               size=0.3, alpha=0.1,outlier.shape = NA) + # outliers

  
ggplot2::theme_bw() + 
  ggplot2::theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  
  ggplot2::scale_fill_manual(values=geo.colours)+
  ggplot2::scale_color_manual(values=geo.colours)+
  
  ggplot2::ylab("Human Footprint Index (HFI)")+xlab(NULL)+ 
  # ylim(0,50)+
  ggplot2::theme(legend.position="none") +
  ggplot2::labs(title = expression("Geographic HFI - "~italic("C. tropicalis")~"")) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5,vjust = -1.5))
  # coord_flip()

p1_2_raw


# add subset correlations

Ne_raw<-read.csv("../tables/table_Ne_resultes.csv")
Ne_raw$Ne_value<-Ne_raw$Ne_value/1000
colnames(Ne_raw)[1] <- "geo"

geo_avg_HFI<-merged_all_plot_p1_2 %>% 
  dplyr::group_by(geo) %>% 
  dplyr::summarise(mean_HFI = mean(HFI, na.rm = TRUE))

merged_all_plot_p1_2_subset<- merge(geo_avg_HFI, Ne_raw, by = "geo", all.x = TRUE)



p1_2_subset <- ggplot(merged_all_plot_p1_2_subset, aes(x = mean_HFI, y = Ne_value,color = geo)) +
  ggplot2::geom_point(shape = 16, size = 3, alpha = 0.9, stroke = 0.2) +  
  ggplot2::scale_color_manual(values = geo.colours)+
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "darkgrey",alpha = 0.5) +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "mean HFI", y = expression(Ne ~ "x" ~ 10^3)) +
    ggpubr::stat_cor(method = "pearson", color = "black") +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(plot.background = element_rect(fill = "transparent", color = "transparent"))+  # transparent background and frame
      # labs(title = "Pearson cor.")+
  ggplot2::theme(plot.title = element_text(vjust = -0.5,hjust = 0.5),
        axis.title.x = element_text(vjust = 0.9),
        axis.title.y = element_text(vjust = -0.9))+
  ggplot2::theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  
p1_2_subset



## addannotation figure
p1_2<-p1_2_raw + 
  ggplot2::annotation_custom(ggplotGrob(p1_2_subset), 
                             xmin = -0.1, xmax = 2.5, 
                             ymin = 19, ymax = 48)

p1_2

# ggsave("HFI_tro_geo.pdf", plot = p1_2, width = 7.5, height = 3.75, units = "in")





#### 3. plot C.elegans Hawaii vs. non-Hawaii
merged_all_plot_p1_3<-merged_all %>% 
  dplyr::filter(species=="C. elegans")
compaired_ele_Hawaii <- list(c('Hawaii', 'non_Hawaii'))

p1_3<- ggplot(merged_all_plot_p1_3, aes(x=geo, y=HFI,fill=geo_type)) + 
  ggbeeswarm::geom_beeswarm(aes(fill = geo_type), 
                            cex=0.8,
                            shape=21, 
                            size = 0.8, 
                            stroke = 0.01,
                            alpha = 0.8)+ 
    ggplot2::stat_boxplot(geom = 'errorbar', width = 0.2, 
               linewidth=0.2) + # error bar
  ggplot2::geom_boxplot(width = 0.25,
               # fill= HFI_type,
               size=0.3, alpha=0.1,outlier.shape = NA) + # outliers
   stat_compare_means(comparisons = compaired_ele_Hawaii, label = "p.signif", method = "wilcox.test", vjust = 0.7) + #vjust是调节星标的垂直位置
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + #去掉各种背景

  ggplot2::scale_fill_manual(values=geo_type.colours)+
  ggplot2::ylab("Human Footprint Index (HFI)")+xlab(NULL)+ 
  # ylim(0,50)+
  ggplot2::theme(legend.position="none") +
  ggplot2::labs(title = expression(~italic("C. elegans")~" Hawaii vs. non-Hawaii")) +
  ggplot2::theme(plot.title = element_text(hjust = 0.5,vjust = -1.5))

p1_3
# ggsave("HFI_C_elegans_Hawaii.pdf", plot = p1_3, width = 3.75, height = 3, units = "in")




#####  assemble pics
p_all <- ggpubr::ggarrange(
  ggpubr::ggarrange(p1_3, p1_1, labels = c("A", "B"), ncol = 2),
  p1_2, 
  labels = c("", "C"), 
  nrow = 2
)
p_all
ggsave("../plots/HFI_all.pdf", plot = p_all, width = 7.5, height = 7.5, units = "in")






