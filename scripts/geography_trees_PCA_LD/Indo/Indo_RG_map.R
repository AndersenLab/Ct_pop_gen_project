rm(list=ls())

library(maps)
library(ggplot2)
library(dplyr)
library(ggthemes)

#### Install this version of ggrepel when you meet any error
# remove.packages("ggrepel")
# devtools::install_version("ggrepel", version = "0.8.0")
library(ggrepel)

source("../../utilities.R")

world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

geo_and_lineage_raw<-read.csv("../../../processed_data/geo_info/geo_and_lineage.csv")
geo_and_lineage_Malay<-geo_and_lineage_raw %>%
  dplyr::filter(geo == "Indonesia")%>% 
  dplyr::select(isotype,long,lat,geo,lineage)
Malaybox <- c(xmin = 110.3762, ymin = -10.1462, xmax = 123.6027, ymax = -0.0699)

Malay_poly <- map_data("world") %>%
  dplyr::filter(
    long >= 110.3762, long <= 123.6027,
    lat  >=  -10.1462, lat  <=  -0.0699
  )

group_map_Malay <- geo_and_lineage_Malay %>%
  dplyr::filter(long >= Malaybox["xmin"],
                long <= Malaybox["xmax"],
                lat  >= Malaybox["ymin"],
                lat  <= Malaybox["ymax"])

set.seed(123)
group_map2_Malay <- findboxes(group_map_Malay,
                           maxiter = 100,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 0.3,
                           box_padding_y = 0.3,
                           point_padding_x = 0.3,
                           point_padding_y = 0.3,
                           force = 1.9e-6,
                           xlim = c(Malaybox["xmin"], Malaybox["xmax"]),
                           ylim = c(Malaybox["ymin"], Malaybox["ymax"])
)

group_map2_Malay_join <- group_map_Malay %>%
  dplyr::left_join(group_map2_Malay, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Malaybox["ymax"] & y2 >= Malaybox["ymin"])

Malay_pca_map1 <- ggplot() +
  geom_polygon(data = Malay_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  geom_segment(data = group_map2_Malay_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_Malay_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.03, -0.02))+
  coord_quickmap(xlim = c(Malaybox["xmin"], Malaybox["xmax"]),
                 ylim = c(Malaybox["ymin"], Malaybox["ymax"]),
                 expand = FALSE)

Malay_pca_map1

ggsave(Malay_pca_map1, 
       filename = "../../../figures/FigureS25_Indo_map.pdf", 
       width = 7.5, 
       height = 7.5/3219*2674, 
       units = "in")






