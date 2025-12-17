
rm(list = ls())

library(extrafont)
library(dplyr)
library(ggplot2)
library(cowplot)


strain_map<-readRDS("../../processed_data/assemble_figure_1/strain_map.rds")
PCA<-readRDS("../../processed_data/assemble_figure_1/Ct_all_isotypes_PCA_plot_PC1_4.rds")
tree <- readRDS("../../processed_data/assemble_figure_1/Ct_plot_0.9_tree_equal_angle_rotated.rds")


p1 <- strain_map + ggtitle(NULL)
p2 <- PCA       + ggtitle(NULL)
p3 <- tree      + ggtitle(NULL)


blank <- ggplot() + theme_void()


row2 <- plot_grid(
  blank, p2,
  ncol = 2,
  labels = c("C", "D"),
  rel_widths = c(1, 3),
  label_x = c(-0.015,0)
)


combined2 <- plot_grid(
  p1,
  row2,
  p3,
  ncol = 1,
  label_size = 14,
  rel_heights = c(2.153522, 1.8, 3.75)
)
combined2

ggsave("../../plots/Figure1_raw.pdf", combined2, width = 7.5, height = sum(c(2.153522,1.8,3.75)), units = "in")



