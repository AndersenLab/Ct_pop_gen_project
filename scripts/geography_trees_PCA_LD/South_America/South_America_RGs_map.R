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

geo_and_lineage_South<-geo_and_lineage_raw %>%
  dplyr::filter(geo == "South America")%>% 
  dplyr::select(isotype,long,lat,geo,lineage)

Southbox <- c(xmin = -85, ymin = -21, xmax = -37, ymax = 12)

South_poly <- map_data("world") %>%
  dplyr::filter(
    long >= -85, long <= -37,
    lat  >=  -21, lat  <=  12
  )

group_map_South <- geo_and_lineage_South %>%
  dplyr::filter(long >= Southbox["xmin"],
                long <= Southbox["xmax"],
                lat  >= Southbox["ymin"],
                lat  <= Southbox["ymax"])

set.seed(123)
group_map2_South <- findboxes(group_map_South,
                           maxiter = 30,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 2.9,
                           box_padding_y = 2.9,
                           point_padding_x = 1.3,
                           point_padding_y = 1.3,
                           force = 1.9e-5,
                           xlim = c(Southbox["xmin"], Southbox["xmax"]),
                           ylim = c(Southbox["ymin"], Southbox["ymax"])
)

group_map2_South_join <- group_map_South %>%
  dplyr::left_join(group_map2_South, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Southbox["ymax"] & y2 >= Southbox["ymin"])

South_pca_map1 <- ggplot() +
  geom_polygon(data = South_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  geom_segment(data = group_map2_South_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_South_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.05, 0.09))+
  coord_quickmap(xlim = c(Southbox["xmin"], Southbox["xmax"]),
                 ylim = c(Southbox["ymin"], Southbox["ymax"]),
                 expand = FALSE)

ggsave(South_pca_map1, 
       filename = "../../../figures/FigureS21_SA_map.pdf", 
       width = 7.5, 
       height = 7.5/48*33, 
       units = "in")






