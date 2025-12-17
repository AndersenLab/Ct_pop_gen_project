
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


geo_and_lineage_Micro<-geo_and_lineage_raw %>%
  filter(geo == "Micronesia")%>% 
  select(isotype,long,lat,geo,lineage)




Microbox <- c(xmin = 158.0656, ymin = 6.7557, xmax = 158.3875, ymax = 7.0231)

library(maps)
Micro_poly <- map_data("world") %>%
  dplyr::filter(
    long >= 158.0656, long <= 158.3875,
    lat  >=  6.7557, lat  <=  7.0231
  )



group_map_Micro <- geo_and_lineage_Micro %>%
  dplyr::filter(long >= Microbox["xmin"],
                long <= Microbox["xmax"],
                lat  >= Microbox["ymin"],
                lat  <= Microbox["ymax"])

set.seed(123)
group_map2_Micro <- findboxes(group_map_Micro,
                           maxiter = 100,
                           xcol = "long", ycol = "lat",
                           box_padding_x = 0.003,
                           box_padding_y = 0.003,
                           point_padding_x = 0.003,
                           point_padding_y = 0.003,
                           force = 1.9e-6,
                           xlim = c(Microbox["xmin"], Microbox["xmax"]),
                           ylim = c(Microbox["ymin"], Microbox["ymax"])
)

group_map2_Micro_join <- group_map_Micro %>%
  dplyr::left_join(group_map2_Micro, by = c("isotype" = "isotype")) %>%
  dplyr::filter(y2 <= Microbox["ymax"] & y2 >= Microbox["ymin"])

Micro_pca_map1 <- ggplot() +
  geom_polygon(data = Micro_poly, aes(x = long, y = lat, group = group),
               fill = "light grey", color = "black", linewidth = 0.25) +
  # connecting segments from real coord -> boxed coord
  geom_segment(data = group_map2_Micro_join,
               aes(x = long, y = lat, xend = x2, yend = y2),
               colour = "black", linewidth = 0.15) +
  geom_point(data = group_map2_Micro_join,
             aes(x = x2, y = y2, fill = as.character(lineage)),
             shape = 21, size = 3, stroke = 0.15) +
  theme_map() +
  scale_fill_manual(name = "Group", values = lineage_colors) +
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.05, 0.09))+
  coord_quickmap(xlim = c(Microbox["xmin"], Microbox["xmax"]),
                 ylim = c(Microbox["ymin"], Microbox["ymax"]),
                 expand = FALSE)

Micro_pca_map1

ggsave(Micro_pca_map1, 
       filename = "../../plots/FigureS14_Ct_Micro_map.pdf", 
       width = 7.5, 
       height = 7.5/3219*2674, 
       units = "in")






