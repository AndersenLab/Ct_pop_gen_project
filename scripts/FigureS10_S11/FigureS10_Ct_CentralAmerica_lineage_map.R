
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

geo_and_lineage_Cen<-geo_and_lineage_raw %>%
  filter(geo == "Central America")%>% 
  select(isotype,long,lat,geo,lineage)





Cenbox <- c(xmin = -85, ymin = 6, xmax = -79, ymax = 12)


library(maps) 
Cen_poly <- map_data("world") %>%
  dplyr::filter(
    long >= -85, long <= -79,
    lat  >=  6, lat  <=  12
  )



group_map_Cen <- geo_and_lineage_Cen %>%
  dplyr::filter(long >= Cenbox["xmin"],
                long <= Cenbox["xmax"],
                lat  >= Cenbox["ymin"],
                lat  <= Cenbox["ymax"])


group_map2_Cen <- findboxes(group_map_Cen,
                           maxiter = 100,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 0.5,
                           box_padding_y = 0.5,
                           point_padding_x = 0.9,
                           point_padding_y = 0.9,
                           force = 1e-6,
                           xlim = c(Cenbox["xmin"], Cenbox["xmax"]),
                           ylim = c(Cenbox["ymin"], Cenbox["ymax"])
)


group_map2_Cen_join <- group_map_Cen %>%
  dplyr::left_join(group_map2_Cen, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Cenbox["ymax"] & y2 >= Cenbox["ymin"])

Cen_pca_map1 <- ggplot() +
  geom_polygon(data = Cen_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  geom_segment(data = group_map2_Cen_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_Cen_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.08, 0.12))+
  coord_quickmap(xlim = c(Cenbox["xmin"], Cenbox["xmax"]),
                 ylim = c(Cenbox["ymin"], Cenbox["ymax"]),
                 expand = FALSE)

Cen_pca_map1

ggsave(Cen_pca_map1, 
       filename = "../../plots/FigureS10_Ct_Cen_map.pdf", 
       width = 7.5, 
       height = 7.5, 
       units = "in")






