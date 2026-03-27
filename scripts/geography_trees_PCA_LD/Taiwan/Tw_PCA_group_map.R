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

group_map<-readRDS("../../../processed_data/geo_info/Taiwan/cdf9_7no.rds")
PC_plots <- readRDS(file = "../../../processed_data/geo_info/Taiwan/Ct_PCA_Taiwan_Lineage.rds")

lineage_colors <- c(
  Tw1 = "#E41A1C",
  Tw2 = "#FF7F00",
  Tw3 = "#FFFF33",
  Tw4 = "#4DAF4A",
  Tw5 = "#40E0D0",
  Tw6 = "#377EB8",
  Tw7 ="#984EA3"
)

world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

plot_tw <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="light grey", size=0.25)+
  theme_void()+
  lims(x=c(120, 122.5),y=c(21.8, 25.4)) 

twbox <- c(xmin = 120, ymin = 21.8, xmax = 122.5, ymax = 25.4)

tw_poly <- map_data("world") %>%
  dplyr::filter(region == "Taiwan")

group_map_tw <- group_map %>%
  dplyr::filter(long >= twbox["xmin"],
                long <= twbox["xmax"],
                lat  >= twbox["ymin"],
                lat  <= twbox["ymax"])

group_map2_tw <- findboxes(group_map_tw,
                           maxiter = 100,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 0.05,
                           box_padding_y = 0.05,
                           point_padding_x = 0.075,
                           point_padding_y = 0.075,
                           force = 1e-6,
                           xlim = c(twbox["xmin"], twbox["xmax"]),
                           ylim = c(twbox["ymin"], twbox["ymax"])
)

group_map2_tw_join <- group_map_tw %>%
  dplyr::left_join(group_map2_tw, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= twbox["ymax"] & y2 >= twbox["ymin"])

tw_pca_map1 <- ggplot() +
  geom_polygon(data = tw_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", size = 0.25) +
  geom_segment(data = group_map2_tw_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_tw_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual("", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  coord_quickmap(xlim = c(twbox["xmin"], twbox["xmax"]),
                 ylim = c(twbox["ymin"], twbox["ymax"]),
                 expand = FALSE) 

hi_div_fig <- cowplot::plot_grid(PC_plots, tw_pca_map1, ncol = 1, labels = c("", "c"), rel_heights = c(2, 3))
ggsave(hi_div_fig, 
       filename = "../../../figures/FigureS14_Taiwan_PCA_and_map.pdf", 
       width = 7.5, height = 8, units = "in")

