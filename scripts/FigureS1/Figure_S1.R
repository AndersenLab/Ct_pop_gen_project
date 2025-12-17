rm(list = ls())

library(extrafont)
library(dplyr)
library(ggplot2)



isotype_map<-readRDS("../../processed_data/assemble_figure_S1/isotype_map.rds")
# geo_distance<-readRDS("../../processed_data/assemble_figure_S1/Haversine_distance.rds")

geo_distance_isotype<-readRDS("../../processed_data/assemble_figure_S1/p_isotype_hist.rds")
geo_distance_strain<-readRDS("../../processed_data/assemble_figure_S1/p_strain_hist.rds")


library(cowplot)


p1 <- isotype_map + ggtitle(NULL) 
p2 <- geo_distance_isotype       + ggtitle(NULL)
p3 <- geo_distance_strain       + ggtitle(NULL)

distance_plots<-plot_grid(
  p2,p3,
  ncol = 2,
  labels = c("C", "D"),
  label_x = c(-0.015, 0)
  # label_size = 14
)
distance_plots




combined2 <- plot_grid(
  p1, distance_plots,
  ncol = 1,
  # labels = c("", ""),
  # label_size = 14,
  rel_heights = c(2.153522, 2.3)
)


print(combined2)

ggsave("../../plots/FigureS1_raw.pdf", combined2, width = 7.5, height = sum(c(2.153522,2.3)), units = "in")





