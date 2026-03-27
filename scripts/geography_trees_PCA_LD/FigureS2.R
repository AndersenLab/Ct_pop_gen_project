rm(list = ls())

library(extrafont)
library(dplyr)
library(ggplot2)
library(cowplot)

isotype_map<-readRDS("../../processed_data/geo_info/isotype_map.rds")
geo_distance_isotype<-readRDS("../../processed_data/geo_info/p_isotype_hist.rds")
geo_distance_strain<-readRDS("../../processed_data/geo_info/p_strain_hist.rds")

p1 <- isotype_map + ggtitle(NULL) 
p2 <- geo_distance_isotype + ggtitle(NULL)
p3 <- geo_distance_strain + ggtitle(NULL)

distance_plots<-plot_grid(
  p2,p3,
  ncol = 2,
  labels = c("c", "d"),
  label_x = c(-0.015, 0)
)

combined2 <- plot_grid(
  p1, distance_plots,
  ncol = 1,
  rel_heights = c(2.153522, 2.3)
)

ggsave("../../figures/FigureS2_isotype_and_geo_distance.pdf", 
       combined2, 
       width = 7.5, 
       height = sum(c(2.153522,2.3)), 
       units = "in")
