
rm(list=ls())


library(maps)
library(ggplot2)
library(dplyr)
library(ggthemes)

#### Install this version of ggrepel when you meet any error
# remove.packages("ggrepel")
# devtools::install_version("ggrepel", version = "0.8.0")
library(ggrepel)


source("../utilities.R")

world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica



geo_and_lineage_raw<-read.csv("../../processed_data/geo_info/geo_and_lineage.csv")


geo_and_lineage_Af<-geo_and_lineage_raw %>%
  filter(geo == "Africa")%>% 
  select(isotype,long,lat,geo,lineage)






############## panel 1 ##############

Afbox <- c(xmin = -31.20, ymin = -6, xmax = 24.62, ymax = 25)

library(maps)         # map_data
Af_poly <- map_data("world") %>%
  dplyr::filter(
    long >= -31.20, long <= 24.62,
    lat  >=  -6, lat  <=  25
  )



group_map_Af <- geo_and_lineage_Af %>%
  dplyr::filter(long >= Afbox["xmin"],
                long <= Afbox["xmax"],
                lat  >= Afbox["ymin"],
                lat  <= Afbox["ymax"])

set.seed(123)
group_map2_Af <- findboxes(group_map_Af,
                           maxiter = 1000,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 2.3,
                           box_padding_y = 2.3,
                           point_padding_x = 1.3,
                           point_padding_y = 1.3,
                           force = 1.0e-6,
                           xlim = c(Afbox["xmin"], Afbox["xmax"]),
                           ylim = c(Afbox["ymin"], Afbox["ymax"])
)

group_map2_Af_join <- group_map_Af %>%
  dplyr::left_join(group_map2_Af, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Afbox["ymax"] & y2 >= Afbox["ymin"])

Af_pca_map1 <- ggplot() +
  geom_polygon(data = Af_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  geom_segment(data = group_map2_Af_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_Af_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.05, 0.09))+
  coord_quickmap(xlim = c(Afbox["xmin"], Afbox["xmax"]),
                 ylim = c(Afbox["ymin"], Afbox["ymax"]),
                 expand = FALSE)+
  geom_text(aes(x = -24, y = 12, label = paste0("Cape Verde", "\n", "Islands")),
            fontface = "bold", size = 5,inherit.aes = FALSE)+
  geom_text(aes(x = -1, y = 0, label = "Sao Tome Island"),
            fontface = "bold", size = 5,inherit.aes = FALSE)
  

Af_pca_map1







############## panel 2 ##############

Afbox <- c(xmin = 40.7, ymin = -23, xmax = 59, ymax = -15)

library(maps)         # map_data
Af_poly <- map_data("world") %>%
  dplyr::filter(
    long >= 40.7, long <= 59,
    lat  >=  -23, lat  <= -15
  )



group_map_Af <- geo_and_lineage_Af %>%
  dplyr::filter(long >= Afbox["xmin"],
                long <= Afbox["xmax"],
                lat  >= Afbox["ymin"],
                lat  <= Afbox["ymax"])

set.seed(123)
group_map2_Af <- findboxes(group_map_Af,
                           maxiter = 1000,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 0.8,
                           box_padding_y = 0.8,
                           point_padding_x = 0.5,
                           point_padding_y = 0.5,
                           force = 1.0e-6,
                           xlim = c(Afbox["xmin"], Afbox["xmax"]),
                           ylim = c(Afbox["ymin"], Afbox["ymax"])
)

group_map2_Af_join <- group_map_Af %>%
  dplyr::left_join(group_map2_Af, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Afbox["ymax"] & y2 >= Afbox["ymin"])

Af_pca_map2 <- ggplot() +
  geom_polygon(data = Af_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  geom_segment(data = group_map2_Af_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_Af_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.06, 0.09))+
  coord_quickmap(xlim = c(Afbox["xmin"], Afbox["xmax"]),
                 ylim = c(Afbox["ymin"], Afbox["ymax"]),
                 expand = FALSE)+
  geom_text(aes(x = 52.5, y = -17, label = "Nosy Boraha Island"),
            fontface = "bold", size = 5,inherit.aes = FALSE)+
  geom_text(aes(x = 53.2, y = -21, label = "Reunion Island"),
            fontface = "bold", size = 5,inherit.aes = FALSE)

Af_pca_map2

# ggsave(Af_pca_map2, 
#        filename = "Ct_Af_map_2.pdf", 
#        width = 7.5, 
#        height = 7.5/3219*2674, 
#        units = "in")














library(cowplot)

Af_map_combined <- plot_grid(
  Af_pca_map1,
  Af_pca_map2,
  ncol = 1,
  labels = c("A", "B"),
  align = "v",
  axis = "tb" 
)

Af_map_combined

ggsave(
  filename = "../../plots/FigureS17_Ct_Af_maps_combined.pdf",
  plot = Af_map_combined,
  width = 7.5,
  height = 8,
  units = "in"
)






