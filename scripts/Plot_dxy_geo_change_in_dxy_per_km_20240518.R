# 

rm(list=ls())


library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(geosphere)
library(ggpubr)
library(grid)


# load all_518_paired_dxy_summary
## note that the all_518_paired_dxy_summary.txt file was merged from all chunks during the dxy calculation
## Hence, we still need to calculate the average dxy for each pair
all_paired_dxy_raw<-read.table(file = "../processed_data/fst_dxy/dxy/chunks/C_tro_all_518_paired_dxy_summary.txt",
                           sep = ' ',
                           header = FALSE)
colnames(all_paired_dxy_raw)<-c("isotype1","isotype2","dxy")

all_paired_dxy<-all_paired_dxy_raw %>%
  dplyr::group_by(isotype1,isotype2) %>% 
  dplyr::summarize(ave_dxy=mean(dxy))




###### now, load geo data ######
# Load geo data
geo_info_raw<- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv", header = TRUE,row.names = 1)
colnames(geo_info_raw)
geo_info<-geo_info_raw %>%
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
  # dplyr::mutate(paired = stringr::str_c(Var1,Var2,sep = "_") ) %>%
  # dplyr::select(paired,value) %>%
  dplyr::mutate(value=value/1000) %>% # change the unite into kilometer (km)
  dplyr::rename(geo_distance=value)
colnames(geo_distance)<-c("isotype1","isotype2","geo_distance")




# merge the two data sets all_paired_dxy & geo_distance

merged_df_tmp_1 <- left_join(geo_distance, all_paired_dxy, by = c("isotype1", "isotype2"))
avg_dxy_mean_tmp<-all_paired_dxy %>%
  rename(isotype2=isotype1, isotype1=isotype2)

merged_df_tmp_2 <- left_join(geo_distance, avg_dxy_mean_tmp, by = c("isotype1", "isotype2"))

merged_df_tmp_3<- merge(merged_df_tmp_1,merged_df_tmp_2, by= c("isotype1","isotype2"))
colnames(merged_df_tmp_3)
###mean_avg_dxy is the previously calculated average number
merged_df_tmp_4 <- merged_df_tmp_3 %>%
  mutate(merged_avg_dxy = paste(ifelse(is.na(ave_dxy.x), "", ave_dxy.x),
                                ifelse(is.na(ave_dxy.y), "", ave_dxy.y), sep = ""))
merged_df_tmp_4$merged_avg_dxy<-as.numeric(merged_df_tmp_4$merged_avg_dxy)
colnames(merged_df_tmp_4)

merged_df<- merged_df_tmp_4 %>% 
  select("isotype1", "isotype2", "geo_distance.x", "merged_avg_dxy") %>%
  rename(geo_distance = geo_distance.x) 


# add geo info into the data set 

geo_info_raw_tmp_1<- geo_info_raw
geo_info_raw_tmp_1$isotype<-row.names(geo_info_raw_tmp_1)

merged_df_plot_tmp_1 <- left_join(merged_df, geo_info_raw_tmp_1, by = c("isotype1" = "isotype")) %>%
  rename(geo1 = geo)  

merged_df_plot_tmp_2 <- left_join(merged_df_plot_tmp_1, geo_info_raw_tmp_1, by = c("isotype2" = "isotype")) %>%
  rename(geo2 = geo)  

# use mutate() to add new col "geo"
merged_df_plot <- merged_df_plot_tmp_2 %>%
  mutate(
    geo = ifelse(geo1 == geo2, geo1, "others")
  ) %>%
  select(-geo1, -geo2)

# View(merged_df_plot)



### we use merged_df_plot data.frame to plot all pairs ###

p1 <- ggplot2::ggplot(merged_df_plot, aes(x = geo_distance, y = merged_avg_dxy)) +
  geom_point(shape = 19, size = 0.5, alpha = 0.4, stroke = 0.2, color= "lightgrey") +  
  geom_smooth(method = "lm", se = FALSE, color = "black",linewidth = 0.5) +
  theme_classic()+
  labs(x = "Geographic distance (km)", y = "Average dxy (%)") +
  # labs(title = "Geo. Distance vs. Dxy Correlation")+
  # ggpubr::stat_cor(method = "pearson", color = "black") +
  ggpubr::stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson",label.sep = "\n", size = 2)+
  theme(text = element_text(size = 8)) +
  labs(title = "All")+
  theme(plot.title = element_text(size = 7))

p1

# ggsave("p.pdf", plot = p, width = 3.75, height = 2.5, units = "in")








### Then, we use data_without_others data.frame to plot all pairs ###
data_without_others <- merged_df_plot %>%
  filter(geo != "others")

geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")


# write a function
plot_custom_geo <- function(geo_value, geo_color) {
  library(dplyr)
  p <- ggplot(data_without_others %>%
                 filter(geo == geo_value), aes(x = geo_distance, y = merged_avg_dxy)) +
    geom_point(shape = 19, size = 0.5, alpha = 0.4, stroke = 0.2, color = geo_color) +  
    geom_smooth(method = "lm", se = FALSE, color = "black",linewidth = 0.5) +
    theme_classic() +
    labs(x = "Geographic distance (km)", y = "Average dxy (%)") +
    # ggpubr::stat_cor(method = "pearson", color = "black") +
    ggpubr::stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson",label.sep = "\n", size = 2)+
    theme(text = element_text(size = 8)) +
    theme(legend.position = "none") +
    labs(title = geo_value)+
    theme(plot.title = element_text(size = 7))
    
    
return(p)
  
}

p2<-plot_custom_geo(geo_value = "Africa",geo_color = "#A6D854")
p3<-plot_custom_geo(geo_value = "Australia",geo_color = "#FC8D62")
p4<-plot_custom_geo(geo_value = "Caribbean",geo_color = "#FFD92F")
p5<-plot_custom_geo(geo_value = "Central America",geo_color = "#8DA0CB")
p6<-plot_custom_geo(geo_value = "Hawaii",geo_color = "#66C2A5")
p7<-plot_custom_geo(geo_value = "South America",geo_color = "#E78AC3")
p8<-plot_custom_geo(geo_value = "Taiwan" ,geo_color = "#E5C494")
p2
p3
p4
p5
p6
p7
p8


# Assemble the 8 plots
p <- ggpubr::ggarrange(p1+ rremove("ylab") + rremove("xlab"),
                       p2+ rremove("ylab") + rremove("xlab"),
                       p3+ rremove("ylab") + rremove("xlab"),
                       p4+ rremove("ylab") + rremove("xlab"),
                       p5+ rremove("ylab") + rremove("xlab"),
                       p6+ rremove("ylab") + rremove("xlab"), 
                       p7+ rremove("ylab") + rremove("xlab"),
                       p8+ rremove("ylab") + rremove("xlab"),
                       # labels = c("All","Africa", "Australia", "Caribbean", "Central America", "Hawaii", "South America", "Taiwan"),
                       vjust = 3,
                       ncol = 4, nrow = 2,
                       align = "hv",
                       common.legend = TRUE, legend = "top",
                       font.label = list(size = 8, color = "black", face = "bold", family = NULL, position = "top"))
pp<-annotate_figure(p,left = grid::textGrob("Dxy (%)", rot = 90, vjust = 0.5, gp = gpar(cex = 1.3,fontsize = 8)), # ,fontface="bold"
                    bottom = textGrob("Geographic distance (km)", gp = gpar(cex = 1.3,fontsize = 8))) # ,fontface="bold" 


pp

ggsave("plot_dxy_geo.pdf", plot = pp, width = 7.5, height = 3, units = "in")








##### plot change in dxy per km #####

# calculate change in dxy per km
plot_dxy_change_per_km<-data_without_others
plot_dxy_change_per_km$dxy_change_per_km <- plot_dxy_change_per_km$merged_avg_dxy / plot_dxy_change_per_km$geo_distance

# Some pairs has geo_distance == 0
# So we need to remove the ones "dxy_change_per_km == Inf" 

plot_dxy_change_per_km<- plot_dxy_change_per_km %>%
  dplyr::filter(dxy_change_per_km != Inf) %>%
  dplyr::select(geo,dxy_change_per_km)

# differences are too big, log10 the data

plot_dxy_change_per_km$dxy_change_per_km<-log10(plot_dxy_change_per_km$dxy_change_per_km)


avg_dxy <- plot_dxy_change_per_km %>%
  group_by(geo) %>%
  summarize(avg_dxy = mean(dxy_change_per_km)) %>%
  arrange(desc(avg_dxy))

plot_dxy_change_per_km$geo <- factor(plot_dxy_change_per_km$geo, levels = avg_dxy$geo)


library(gghalves)
# plot 
p_box<-ggplot2::ggplot(plot_dxy_change_per_km,aes(x= geo, y= dxy_change_per_km,fill=geo,color=geo))+
  ggplot2::scale_fill_manual(values = geo.colours)+
  ggplot2::scale_colour_manual(values = geo.colours)+
  ggplot2::geom_violin(alpha=0.2,linewidth = 0.5)+
  ggplot2::geom_boxplot(width =0.2,alpha=0.6, outlier.size = 0.5)+
  ggplot2::theme_bw()+
  ggplot2::theme(panel.grid=element_blank())+
  labs(x="Geographic regions",
       y = expression(Log[10]("Dxy change per km")),
       color="Geographic regions",
       fill="Geographic regions")+
  theme(axis.title.y = element_text(size = 9,face = "bold"),
        axis.title.x = element_text(size = 9,face = "bold"),
        # axis.title.y = element_text(face = "bold", size = 9),
        # axis.title.x = element_text(face = "bold", size = 9),
        axis.text = element_text(size = 8))+
  theme(legend.position = "top",
        legend.title = element_text(size = 8,face = "bold"),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(10, "pt"),
        legend.margin = margin(t = 10))+
  # legend.margin = margin(0, 0, -8, 0))+
  guides(fill = guide_legend(ncol = 7))  
# legend.key.width = unit(5, "pt"),
# legend.key.height = unit(5, "pt"))
  


p_box

# ggsave("plot_dxy_change_per_km.pdf", plot = p_box, width = 7.5, height = 3, units = "in")







plot_dxy_change_per_km_2<-plot_dxy_change_per_km %>%
  dplyr::filter(geo != "Africa" & geo != "Australia")



p_box_2<-ggplot2::ggplot(plot_dxy_change_per_km_2,aes(x= geo, y= dxy_change_per_km,fill=geo,color=geo))+
  ggplot2::scale_fill_manual(values = geo.colours)+
  ggplot2::scale_colour_manual(values = geo.colours)+
  ggplot2::geom_violin(alpha=0.2,linewidth = 0.5)+
  ggplot2::geom_boxplot(width =0.2,alpha=0.6, outlier.size = 0.5)+
  ggplot2::theme_bw()+
  ggplot2::theme(panel.grid=element_blank())+
  labs(x="Geographic regions",
       y = expression(Log[10]("Dxy change per km")),
       color="Geographic regions",
       fill="Geographic regions")+
  theme(axis.title.y = element_text(size = 9), # ,face = "bold"
        axis.title.x = element_text(size = 9), # ,face = "bold"
        # axis.title.y = element_text(face = "bold", size = 9),
        # axis.title.x = element_text(face = "bold", size = 9),
        axis.text = element_text(size = 8))+
  theme(legend.position = "top",
        legend.title = element_text(size = 8,face = "bold"),
        legend.text  = element_text(size = 6),
        legend.key.size = unit(10, "pt"),
        legend.margin = margin(t = 10))+
  # legend.margin = margin(0, 0, -8, 0))+
  guides(fill = guide_legend(ncol = 7))  
# legend.key.width = unit(5, "pt"),
# legend.key.height = unit(5, "pt"))


p_box_2

# ggsave("plot_dxy_change_per_km_2.pdf", plot = p_box_2, width = 7.5, height = 3, units = "in")





p_all <- ggpubr::ggarrange(pp + theme(plot.margin = margin(t = 10, l = 10)), # add margin to make sure label have enough space in the top left (top margin & left margin)
                           p_box_2 + theme(plot.margin = margin(t = 10, l = 10)), # same way to add margin
                           ncol = 1, nrow = 2,
                           labels = c("a","b"),
                           label.x = 0.01, label.y = 0.95, # readjust label.x and label.y to change the labels' position
                           font.label = list(size = 12, face = "bold")) # font style
p_all
# ggsave("plot_dxy_all.pdf", plot = p_all, width = 7.5, height = 6, units = "in")


#####





##### assemble figures #####

p_all <- ggpubr::ggarrange(pp,p_box_2,
                       vjust = 20,
                       ncol = 1, nrow = 2,
                       labels = c("a","b"))


p_all

ggsave("../plots/plot_dxy_all.pdf", plot = p_all, width = 7.5, height = 6, units = "in")


##### end assemble figures #####









