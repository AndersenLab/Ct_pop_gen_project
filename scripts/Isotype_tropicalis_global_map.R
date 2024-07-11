## Modified from Lee et.al. 2021, Nature E&E
rm(list=ls())

# library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(maps)
library(ggpubr)
# library(gsheet)


# Load input data. File was downloaded on April 2nd 2024.
# raw_xlsx<-readxl::read_excel("../data/C. tropicalis WI strain info.xlsx")
# # Or load the data from google drive 
# raw_xlsx <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/edit")

raw_csv<-read.csv("../data/20231201_c_tropicalis_strain_data.csv")

indep_isotype_info<-raw_csv %>% 
dplyr::filter(strain == isotype)  %>% 
  dplyr::select(isotype,latitude,longitude)%>%
  dplyr::filter(!is.na(latitude) & !is.na(longitude)) %>%
  dplyr::rename("lat"="latitude", "long"="longitude") 

indep_isotype_info$lat<-as.numeric(indep_isotype_info$lat)
indep_isotype_info$long<-as.numeric(indep_isotype_info$long)


# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# all_tro_isotype<- readxl::read_excel("strains_sequenced_and_done.xlsx",1)


# we don't use this form anymore since I filter the sequenced and done isotype (n=707)
# all_tro_isotype<-read.csv(file = "C. tropicalis WI isotype info - C. tropicalis WI isotype info.csv",
#                           header = TRUE,
#                           sep = ',')


#add the two missing data to the latest isotype list “ECA2954” and “ECA2955"
# 
# indep_isotype_info<-all_tro_isotype %>%
#   dplyr::mutate(isotype = ifelse(strain %in% c("ECA2954", "ECA2955"), strain, isotype)) %>%
#   dplyr::filter(strain == isotype)  %>% 
#   dplyr::select(isotype,latitude,longitude)%>%
#   dplyr::filter(!is.na(latitude) & !is.na(longitude)) %>%
#   dplyr::rename("lat"="latitude", "long"="longitude") 

#1. Central America
df_ca_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -97.8  & long < -76 & lat > 9 & lat < 10.5) 
ca_isotype <- as.character(df_ca_isotype$isotype)

#2. South America
df_sa_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -72  & long < -52 & lat > -13 & lat < 6) 
sa_isotype <- as.character(df_sa_isotype$isotype)

#3. Hawaii
df_hw_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -170  & long < -150)
hw_isotype <- as.character(df_hw_isotype$isotype)

#4. New zealand
df_nz_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 165  & long < 180) %>%
  dplyr::filter(lat > -50 & lat < -25)
nz_isotype <- as.character(df_nz_isotype$isotype)


#5. Australia
df_au_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 110  & long < 160) %>%
  dplyr::filter(lat > -50 & lat < -10) 
au_isotype <- as.character(df_au_isotype$isotype)

#6. Caribbean
df_car_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -80  & long < -50) %>%
  dplyr::filter(lat > 10 & lat < 20) 
car_isotype <- as.character(df_car_isotype$isotype)


#7. Africa
df_af_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -30  & long < 70) %>%
  dplyr::filter(lat > -40 & lat < 30) 
af_isotype <- as.character(df_af_isotype$isotype)

#8. Europe
df_eu_isotype <- indep_isotype_info %>%
  dplyr::filter(long > -10  & long < 25) %>%
  dplyr::filter(lat > 35 & lat < 57) 
eu_isotype <- as.character(df_eu_isotype$isotype)


#9. Taiwan
df_tw_isotype <- indep_isotype_info %>%
  dplyr::filter(long > 50  & long < 150) %>%
  dplyr::filter(lat > 12 & lat < 75) 
tw_isotype <- as.character(df_tw_isotype$isotype)


#10. others
df_etc_isotype <- indep_isotype_info %>%
  dplyr::filter(!isotype %in% c(af_isotype, au_isotype, car_isotype,
                               eu_isotype,hw_isotype, ca_isotype, 
                               nz_isotype,sa_isotype, tw_isotype))
etc_isotype <- as.character(df_etc_isotype$isotype)



### Merge                  

indep_isotype_info_geo <- indep_isotype_info %>%
  dplyr::mutate(geo = ifelse(isotype %in% hw_isotype, "Hawaii", 
                             ifelse(isotype %in% ca_isotype, "Central America", 
                                    ifelse(isotype %in% au_isotype, "Australia", 
                                           ifelse(isotype %in% af_isotype, "Africa", 
                                                  ifelse(isotype %in% car_isotype, "Caribbean",
                                                         ifelse(isotype %in% eu_isotype, "Europe", 
                                                                ifelse(isotype %in% sa_isotype, "South America", 
                                                                       ifelse(isotype %in% nz_isotype, "New Zealand", 
                                                                              ifelse(isotype %in% tw_isotype, "Taiwan", "Unknown")))))))))) 


unique(indep_isotype_info_geo$geo)

# export indep_isotype_info_geo
write.csv(indep_isotype_info_geo,
          # file = "test.csv",
          file= "../processed_data/Geo_info/indep_isotype_info_geo.csv",
          quote = FALSE,
          row.names = FALSE)


library(RColorBrewer)
display.brewer.all()
display.brewer.pal(7,"Set2")
brewer.pal(7, "Set2")
display.brewer.pal(7, "RdBu")
brewer.pal(7, "RdBu")

geo.colours <- c("Hawaii"="#66C2A5", "Australia"="#FC8D62", "Central America"="#8DA0CB",
                 "South America"="#E78AC3", "Africa"="#A6D854", "Caribbean"="#FFD92F", 
                 "Taiwan" = "#E5C494")



### Geographic distribution of all isotype 
world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica


plot_world_raw <-ggplot2::ggplot()+ geom_map(data=world, map=world,
                                             aes(x=long, y=lat, map_id=region),
                                             color="black", fill="#F7F7F7", size=0.2)+
  geom_point(data = indep_isotype_info_geo, aes(x=long, y=lat, color = geo), shape =16, 
             # alpha = 0.3,
             size = 2, ) +
  scale_color_manual(values = geo.colours) +
  lims(x=c(-180,191),y=c(-58,84)) +
  theme_void()+
  #  theme(text = element_text(size=12)) +
  labs(color = NULL)+
  theme(legend.position = c(0.075, 0.2))+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank(),   # Conceal Tick Marks
        legend.text = element_text(size = 5))+
  theme(legend.key.size = unit(0.4, "cm"))

# ggthemes::theme_map()


plot_world_raw


# the proper wide:height is 2.612
# save image 3*7.83 inches



### Pie charts of the isotype
geo_freq <- indep_isotype_info_geo %>%
  dplyr::group_by(geo) %>%
  dplyr::summarize(frequency = n()) %>%
  arrange(desc(frequency))

pie_plot<-ggplot(geo_freq, aes(x = "", y = frequency, fill = geo)) +
  geom_bar(stat = "identity", width = 1, size = 0.1, color = "white") +
  coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  # geom_text(aes(label = frequency), position = position_stack(vjust = 0.5)) +
  geom_text(aes(x = 1.8, label = frequency, color = geo),
            position = position_stack(vjust = 0.5),
            size = 2,
            hjust = 0.5) +
  scale_fill_manual(values = geo.colours)+
  scale_color_manual(values = geo.colours)
# save image 4*4 inches
pie_plot

# ggsave("isotype_pie_plot.pdf", plot = pie_plot, width = 4, height = 4, units = "in")



plot_world<-plot_world_raw + 
  ggplot2::annotation_custom(ggplotGrob(pie_plot), 
                             xmin = -160, xmax = -80, 
                             ymin = -60, ymax = 20)
plot_world
# ggsave("isotype_world_plot.pdf", plot = plot_world, width = 7.5, height = 2.87, units = "in")






### Plot small maps

# Get frequencies
n_sa<-geo_freq %>% 
  dplyr::filter(geo=="South America") %>% 
  dplyr::select(frequency) %>% 
  pull()

n_hw<-geo_freq %>% 
  dplyr::filter(geo=="Hawaii") %>% 
  dplyr::select(frequency) %>% 
  pull()

n_ca<-geo_freq %>% 
  dplyr::filter(geo=="Central America") %>% 
  dplyr::select(frequency) %>% 
  pull()

n_car<-geo_freq %>% 
  dplyr::filter(geo=="Caribbean") %>% 
  dplyr::select(frequency) %>% 
  pull()

n_tw<-geo_freq %>% 
  dplyr::filter(geo=="Taiwan") %>% 
  dplyr::select(frequency) %>% 
  pull()

# small maps _ Hawaii
plot_hw <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="#F7F7F7", size=0.2)+
  geom_point(data = indep_isotype_info_geo, aes(x=long, y=lat,color = geo), shape =16, 
             # alpha = 0.5,
             size =1, ) +
  scale_color_manual(values = geo.colours) +
  theme_void()+
  
  theme(text = element_text(size=12)) +
  lims(x=c(-160, -155),y=c(18,23)) +
  theme(legend.position = "none")+
  geom_text(aes(x = -160, y = 18, label = paste("Hawaii \nn =",n_hw)), hjust = 0, vjust = 0, size = 2)+
  labs(color = NULL)+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank())   # Conceal Tick Marks

plot_hw


# small maps _ South America
plot_sa <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="#F7F7F7", size=0.2)+
  geom_point(data = indep_isotype_info_geo, aes(x=long, y=lat,color = geo), shape =16, 
             # alpha = 0.5,
             size =1) +
  scale_color_manual(values = geo.colours) +
  theme_void()+
  
  theme(text = element_text(size=12)) +
  lims(x=c(-85, -50),y=c(-27,8)) +
  theme(legend.position = "none")+
  geom_text(aes(x = -85, y = -27, label = paste("South America \nn =",n_sa)), hjust = 0, vjust = 0, size = 2)+
  labs(color = NULL)+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank())   # Conceal Tick Marks

plot_sa






# small maps _ Central America
plot_ca <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="#F7F7F7", size=0.2)+
  geom_point(data = indep_isotype_info_geo, aes(x=long, y=lat,color = geo), shape =16, 
             # alpha = 0.5, 
             size =1 ) +
  scale_color_manual(values = geo.colours) +
  theme_void()+
  labs(color = NULL)+
  
  theme(text = element_text(size=12)) +
  lims(x=c(-85, -79),y=c(6, 12)) +
  theme(legend.position = "none")+
  geom_text(aes(x = -85, y = 6, label = paste("Central America \nn =",n_ca)), hjust = 0, vjust = 0, size = 2)+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank())   # Conceal Tick Marks

plot_ca






# small maps _ Caribbean
plot_car <- ggplot()+ geom_map(data=world, map=world,
                               aes(x=long, y=lat, map_id=region),
                               color="black", fill="#F7F7F7", size=0.2)+
  geom_point(data = indep_isotype_info_geo, aes(x=long, y=lat,color = geo), shape =16, 
             # alpha = 0.5, 
             size =1, ) +
  scale_color_manual(values = geo.colours) +
  theme_void()+
  labs(color = NULL)+
  
  theme(text = element_text(size=12)) +
  lims(x=c(-67, -60),y=c(12, 19)) +
  theme(legend.position = "none")+
  geom_text(aes(x = -67, y = 12, label = paste("Caribbean \nn =",n_car)), hjust = 0, vjust = 0, size = 2)+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank())   # Conceal Tick Marks

plot_car





# small maps _ Taiwan
plot_tw <- ggplot()+ geom_map(data=world, map=world,
                              aes(x=long, y=lat, map_id=region),
                              color="black", fill="#F7F7F7", size=0.2)+
  geom_point(data = indep_isotype_info_geo, aes(x=long, y=lat,color = geo), shape =16, 
             # alpha = 0.5,
             size =1, ) +
  scale_color_manual(values = geo.colours) +
  theme_void()+
  labs(color = NULL)+
  
  theme(text = element_text(size=12)) +
  lims(x=c(119.5, 122),y=c(21.65, 25.15)) +
  theme(legend.position = "none")+
  geom_text(aes(x = 119.5, y = 21.65, label = paste("Taiwan \nn =",n_tw)), hjust = 0, vjust = 0, size = 2)+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank())   # Conceal Tick Marks

plot_tw






# Assemble the 6 plots (except pie chart)
p <- ggpubr::ggarrange(ggpubr::ggarrange(plot_world,ncol = 1, labels = c("a"),
                                         widths = c(1)),
                       ggpubr::ggarrange(plot_hw,
                                         plot_sa,
                                         plot_car,
                                         plot_ca,
                                         plot_tw,
                                         ncol = 5, labels = c("b","c","d","e","f"),
                                         widths = c(rep(0.2,5))),
                       nrow = 2, heights = c(0.657,0.343))
p

#save image 4.37 * 7.5 inches

ggsave("../plots/isotype_only_map.pdf", plot = p, width = 7.5, height = 4.37, units = "in")















