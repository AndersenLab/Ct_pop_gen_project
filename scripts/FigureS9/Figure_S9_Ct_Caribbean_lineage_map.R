
rm(list=ls())


library(maps)
library(ggplot2)
library(dplyr)
library(ggthemes)


source("../utilities.R")

world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica


geo_and_lineage_raw<-read.csv("../../processed_data/geo_info/geo_and_lineage.csv")


geo_and_lineage_car<-geo_and_lineage_raw %>%
  filter(geo == "Caribbean")%>% 
  select(isotype,long,lat,geo,lineage)





Carbox <- c(xmin = -67, ymin = 12, xmax = -59, ymax = 19)


library(maps)
Car_poly <- map_data("world") %>%
  dplyr::filter(
    long >= -67, long <= -59,
    lat  >=  12, lat  <=  19
  )


group_map_Car <- geo_and_lineage_car %>%
  dplyr::filter(long >= Carbox["xmin"],
                long <= Carbox["xmax"],
                lat  >= Carbox["ymin"],
                lat  <= Carbox["ymax"])


group_map2_Car <- findboxes(group_map_Car,
                           maxiter = 100,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 0.5,
                           box_padding_y = 0.5,
                           point_padding_x = 0.9,
                           point_padding_y = 0.9,
                           force = 1e-6,
                           xlim = c(Carbox["xmin"], Carbox["xmax"]),
                           ylim = c(Carbox["ymin"], Carbox["ymax"])
)


group_map2_Car_join <- group_map_Car %>%
  dplyr::left_join(group_map2_Car, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Carbox["ymax"] & y2 >= Carbox["ymin"])

Car_pca_map1 <- ggplot() +
  geom_polygon(data = Car_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  geom_segment(data = group_map2_Car_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_Car_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.08, 0.12))+
  coord_quickmap(xlim = c(Carbox["xmin"], Carbox["xmax"]),
                 ylim = c(Carbox["ymin"], Carbox["ymax"]),
                 expand = FALSE)

Car_pca_map1

ggsave(Car_pca_map1, filename = "../../plots/FigureS9_Ct_Car_map.pdf", width = 7.5, height = 7.5/8*7, units = "in")






