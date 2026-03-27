rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(maps)
library(viridis)
library(pals)
library(sf)
library(ggspatial)

#### Install this version of ggrepel when you meet any error
# remove.packages("ggrepel")
# devtools::install_version("ggrepel", version = "0.8.0")
library(ggrepel)
source("../../utilities.R")

PC_plots <- readRDS(file = "../../../processed_data/geo_info/Hawaii/Ct_PCA_Hawaii_RG.rds")

hi <- sf::st_read("../../../data/Hawaii/cb_2018_us_state_500k") %>%
  dplyr::filter(NAME == "Hawaii")
hibox = c(xmin = -160.00, ymin = 18.7, xmax = -154.4, ymax = 22.4)
hi_crop <- sf::st_crop(hi, hibox)

cdf9_7 <- readRDS(file = "../../../processed_data/geo_info/Hawaii/cdf9_7no.rds")

group_map <- cdf9_7 %>%
  dplyr::select(isotype, lineage, 
                PC1:PC6,latitude,longitude) %>% 
  rename(lat=latitude,long=longitude) %>% 
  as.data.frame()

group_map2 <- findboxes(group_map,
                        maxiter = 100,
                        xcol = "long", ycol = "lat",
                        box_padding_x = 0.05,
                        box_padding_y = 0.05,
                        point_padding_x = 0.075,
                        point_padding_y = 0.075,
                        force = 1e-6,
                        xlim = c(-159.87,-154.76),
                        ylim =c(18.03, 23.0) 
)

group_map2_join <- group_map %>%
  dplyr::left_join(group_map2, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 < 23 & y2 > 18)

hi_pca_map1 <- ggplot() + geom_sf(data = hi_crop, size = 0.25, fill = "light grey") +
  geom_segment(data = group_map2_join, aes(x = long, y = lat, xend = x2, yend = y2), colour = 'black', linewidth = 0.15) +
  geom_point(data = group_map2_join,
             aes(x=x2, y=y2, fill = as.character(lineage)), shape =21, size = 3, stroke = 0.15) +
  theme_map()+
  scale_fill_manual("", values = lineage_colors) +
  theme(text = element_text(size=12), legend.position = "none") 

hi_div_fig <- cowplot::plot_grid(PC_plots, hi_pca_map1, ncol = 1, labels = c("", "c"), rel_heights = c(2, 3))

ggsave(hi_div_fig, filename = "../../../figures/FigureS16_Hw_PCA_and_map.pdf",
       width = 7.5, height = 8.258929, units = "in")
