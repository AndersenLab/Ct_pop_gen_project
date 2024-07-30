# 

rm(list=ls())


library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(pheatmap)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(8,"Set3")
brewer.pal(8,"Set3")



# load species_wide_dxy_summary
## note that the all_518_paired_dxy_summary.txt file was merged from all chunks during the dxy calculation

species_wide_dxy<-read.table(file = "../processed_data/fst_dxy/dxy/chunks/C_tro_species_wide_dxy_summary.txt",
                           sep = ' ',
                           header = FALSE)
species_wide_dxy<- species_wide_dxy %>%
  dplyr::filter(V1 != "chromosome")

colnames(species_wide_dxy)<-c("chrom","window_start_dxy","window_stop_dxy","dxy")
species_wide_dxy$window_start_dxy<-as.numeric(species_wide_dxy$window_start_dxy)
species_wide_dxy$window_stop_dxy<-as.numeric(species_wide_dxy$window_stop_dxy)

species_wide_dxy <- species_wide_dxy %>%
  dplyr::mutate(window_start_kb = round(window_start_dxy / 1000),
                window_stop_kb=round(window_stop_dxy / 1000))





##### p1 - plot genome wide dxy


plot_species_wide_dxy<-species_wide_dxy
plot_species_wide_dxy$x <- (plot_species_wide_dxy$window_start_dxy + plot_species_wide_dxy$window_stop_dxy) / 2
plot_species_wide_dxy<-plot_species_wide_dxy%>%
  dplyr::group_by(chrom) %>%
  dplyr::mutate(x_normalized = (x - min(x)) / (max(x) - min(x)))


p1 <- ggplot2::ggplot(data = plot_species_wide_dxy, mapping = aes(x = x_normalized, y = dxy))+
  geom_point(color = "lightgray", size = 0.02, alpha = 0.8)+
  facet_grid(. ~ chrom, scales = 'free')+
  xlab("Normalized genome position")+
  ylab("Dxy (%)")+
  geom_smooth(method = "loess", se = FALSE,alpha = 0.8,color = "#BEBADA") +
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(strip.background = element_blank())+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.title.y = element_text(face = "bold", size = 9))+
  theme(axis.title.x = element_text(face = "bold", size = 9))

p1
ggsave("../plots/species_wide_dxy.pdf", plot = p1, width = 7.5, height = 2.5, units = "in")












###### p2 - plot heat map of dxy



heatmap_species_wide_dxy<-plot_species_wide_dxy
chrom_regions<-read.table("../processed_data/Chrom_regions/ct_ch_domains.tsv",
                          header = TRUE)
chrom_regions <-chrom_regions %>%
  dplyr::rename(chrom=CHROM)
chrom_regions$start<-as.numeric(chrom_regions$start)
chrom_regions$end<-as.numeric(chrom_regions$end)



colnames(heatmap_species_wide_dxy)







# add location column in heatmap_species_wide_dxy
heatmap_species_wide_dxy <- heatmap_species_wide_dxy %>%
  mutate(location = NA)

# for every row of chrom_regions
for (i in 1:nrow(chrom_regions)) {
  # extract every row of chrom_regions 
  region <- chrom_regions[i, ]
  
  # in heatmap_species_wide_dxy, find certain raw，set location
  heatmap_species_wide_dxy <- heatmap_species_wide_dxy %>%
    mutate(location = ifelse(
      chrom == region$chrom & x >= region$start & x <= region$end,
      region$Location,
      location
    ))
}


# View(heatmap_species_wide_dxy)



plot_heatmap_dxy<-heatmap_species_wide_dxy %>%
  group_by(chrom, location) %>%
  summarize(mean_dxy = mean(dxy, na.rm = TRUE)) %>%
  na.omit()


plot_heatmap_dxy_wide<- reshape2::dcast(plot_heatmap_dxy, chrom ~ location, value.var = "mean_dxy")
row.names(plot_heatmap_dxy_wide)<-plot_heatmap_dxy_wide$chrom
plot_heatmap_dxy_wide<-plot_heatmap_dxy_wide[,-1]
plot_heatmap_dxy_wide <- plot_heatmap_dxy_wide %>%
  dplyr::mutate(Full = rowMeans(.))
plot_heatmap_dxy_wide<-plot_heatmap_dxy_wide%>%
  rename(Center=center,
         "Left arm"=left_arm,
         "Left tip"=left_tip,
         "Right arm"=right_arm,
         "Right tip"=right_tip,)


bk1 = unique(c(seq(0.03,0.15, length=50)))
p2<- pheatmap::pheatmap(plot_heatmap_dxy_wide, 
                        scale = 'none',
                        cluster_cols = F,
                        breaks = bk1,
                        cluster_row = F,
                        show_colnames     = T,
                        show_rownames     = T,
                        main = "Dxy (%)",
                        color = colorRampPalette(c("#F7FBFF", 
                                                   "#BEBADA"))(50), 
                        legend = T,
                        border_color = 'lightgray',
                        fontsize = 6,
                        fontsize_col = 6,
                        angle_col=45)

ggsave("../plots/Dxy_heatmap.pdf", plot = p2, width = 3.75, height = 2.5, units = "in")






######## p3 - plot Dxy vs. pi Correlation

# load pi stat
window_diversity<-read.csv("../processed_data/pi_theta_d/chromosome_windows_diversity.csv",
                           header = TRUE,
                           row.names = 1)
species_wide_pi<-window_diversity %>%
  dplyr::filter(stat_type == "pi") %>%
  dplyr::select(chrom,x,window_start,window_stop,stat) %>%
  dplyr::rename(pi=stat,window_start_pi=window_start,window_stop_pi=window_stop) %>%
  dplyr::mutate(window_start_kb = round(window_start_pi / 1000),
         window_stop_kb= round(window_stop_pi / 1000))



# set threshold for merge (1 kb)
merge_condition <- function(x, y, tolerance = 1) {
  abs(x - y) <= tolerance
}

# merge under this threshold
merged_dxy_pi <- species_wide_dxy %>%
  left_join(species_wide_pi, 
            by = c("chrom"),
            suffix = c("_dxy", "_pi")) %>%
  filter(merge_condition(window_start_kb_dxy, window_start_kb_pi) &
           merge_condition(window_stop_kb_dxy, window_stop_kb_pi))





p3 <-ggplot2::ggplot(merged_dxy_pi, aes(x = pi, y = dxy)) +
    geom_point(size = 1, alpha = 0.4, stroke = 0, color = "black") +  
    geom_smooth(method = "lm", se = FALSE, color="red") +
    theme_classic() +
    # labs(x = "Geographic distance (km)", y = "Average dxy (%)") +
    ggpubr::stat_cor(method = "pearson", color = "black") +
    theme(legend.position = "none") +
    labs(title = "Dxy vs. \u03C0 Correlation")+
    labs(x="\u03C0",y="Dxy (%)")+
    ylim(0, 0.6)
p3

ggsave("../plots/Dxy_pi_cor_raw.pdf", plot = p3, width = 3.75, height = 2.5, units = "in")


# 
# # filter the details about the two outliers (low pi but high Dxy value):
# 
# two_outliers<-merged_dxy_pi %>% 
#   dplyr::filter(pi<0.003 & dxy >0.35) %>% 
#   dplyr::select(chrom, window_start_dxy, window_stop_dxy,
#                 dxy, window_start_pi, window_stop_pi, pi)
# 
# write.csv(file = "two_outliers.csv",two_outliers,quote = FALSE,
#           row.names = FALSE)
# 

  




