rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(maps)

source("../utilities.R")

Ct_raw_data<-read.csv("../../data/20250627_c_tropicalis_strain_data.csv")
Cb_raw_data<-read.csv("../../data/meta_data_Ce_Cb/Cb/20250626_c_briggsae_strain_data.csv")
Ce_raw_data<-read.csv("../../data/meta_data_Ce_Cb/Ce/20250625_c_elegans_strain_data.csv")

Ct_data<-Ct_raw_data %>% 
  dplyr::select(species,strain,latitude,longitude)
Cb_data<-Cb_raw_data %>% 
  dplyr::select(species,strain,latitude,longitude)
Ce_data<-Ce_raw_data %>% 
  dplyr::select(species,strain,latitude,longitude)

all_data<-rbind(Ct_data,Cb_data,Ce_data) %>% 
  na.omit()

world <- map_data("world")
world <- world[world$region != "Antarctica", ]

binned <- all_data %>%
  dplyr::mutate(
      lon_bin = round(longitude, 0),
      lat_bin = round(latitude, 0)
    ) %>%
  dplyr::count(species, lon_bin, lat_bin, name = "n")

binned <- binned %>%
  dplyr::mutate(n_cap = pmin(n, 40))

binned$species <- factor(
  binned$species,
  levels = c(
    "Caenorhabditis elegans",
    "Caenorhabditis briggsae",
    "Caenorhabditis tropicalis"
  )
)

plot_world <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(x = long, y = lat, map_id = region),
    color = "gray95", fill = "gray85", linewidth = 0.05
  ) +
  geom_point(
    data = binned,
    aes(x = lon_bin, y = lat_bin,
        size = n_cap, color = species),
    alpha = 0.5,shape =16
  ) +
  scale_size_continuous(
    limits = c(1, 40),
    range  = c(0.25, 5),
    breaks = c(1, 10, 20, 30, 40),
    labels = c("1", "10", "20", "30", ">=40"), 
    name   = "Strains per \n1° x 1° grid"
  ) +
  scale_color_manual(values = c(
    "Caenorhabditis elegans"    = "#DB6333",
    "Caenorhabditis briggsae"   = "#53886C",
    "Caenorhabditis tropicalis" = "#0719BC"
  )) +
  facet_wrap(~ species, ncol = 1) +
  coord_quickmap(xlim = c(-180, 180),
                 ylim = c(-60, 85),
                 expand = FALSE) +
  theme_void() +
  theme(legend.position = "right") +
  guides(color = "none")+
  theme(
    legend.position = "right",
    strip.text = element_text(
      size = 9,
      face = "bold.italic",
      margin = margin(t = 4, b = 4)
    )
  )

ggsave("../../figures/FigureS1_all_strains_map_3_species.pdf", plot = plot_world, width = 7, height = 8, units = "in")





