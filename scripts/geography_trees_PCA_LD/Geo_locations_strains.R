rm(list = ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(maps)
library(ggpubr)
library(RColorBrewer)

source("../utilities.R")

raw_data<-read.csv("../../data/20250627_c_tropicalis_strain_data.csv")
indep_strain_info<- raw_data

indep_strain_info<-raw_data %>% 
  dplyr::rename(lat=latitude) %>% 
  dplyr::rename(long=longitude) %>% 
  dplyr::select(strain,isotype,lat,long)

indep_strain_info$lat<-as.numeric(indep_strain_info$lat)
indep_strain_info$long<-as.numeric(indep_strain_info$long)

#1. Central America
df_ca_strain <- indep_strain_info %>%
  dplyr::filter(long > -97.8  & long < -76 & lat > 9 & lat < 10.5) 
ca_strain <- as.character(df_ca_strain$isotype)

#2. South America
df_sa_strain <- indep_strain_info %>%
  dplyr::filter(long > -72  & long < -52 & lat > -13 & lat < 6) 
sa_strain <- as.character(df_sa_strain$isotype)

#3. Hawaii
df_hw_strain <- indep_strain_info %>%
  dplyr::filter(long > -176  & long < -129 & lat > 4 & lat < 46)  
hw_strain <- as.character(df_hw_strain$isotype)

#5. Australia
df_au_strain <- indep_strain_info %>%
  dplyr::filter(long > 110  & long < 160) %>%
  dplyr::filter(lat > -50 & lat < -10) 
au_strain <- as.character(df_au_strain$isotype)

#6. Caribbean
df_car_strain <- indep_strain_info %>%
  dplyr::filter(long > -80  & long < -50) %>%
  dplyr::filter(lat > 10 & lat < 20) 
car_strain <- as.character(df_car_strain$isotype)

#7. Africa
df_af_strain <- indep_strain_info %>%
  dplyr::filter(long > -30  & long < 70) %>%
  dplyr::filter(lat > -40 & lat < 30) 
af_strain <- as.character(df_af_strain$isotype)

#9. Taiwan
df_tw_strain <- indep_strain_info %>%
  dplyr::filter(long > 50  & long < 150) %>%
  dplyr::filter(lat > 12 & lat < 75) 
tw_strain <- as.character(df_tw_strain$isotype)

#11. Malay Archipelago
df_mal_strain <- indep_strain_info %>%
  dplyr::filter(long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605)  # Malay Archipelago
mal_strain <- as.character(df_mal_strain$isotype)

#12. Micronesia
df_mic_strain <- indep_strain_info %>%
  dplyr::filter(long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335) # Micronesia
mic_strain <- as.character(df_mic_strain$isotype)

#0. others
df_etc_strain <- indep_strain_info %>%
  dplyr::filter(!isotype %in% c(af_strain, au_strain, car_strain,
                                hw_strain, ca_strain, 
                                sa_strain, tw_strain,
                                mal_strain, mic_strain
  ))
etc_strain <- as.character(df_etc_strain$isotype)



### Merge                  
indep_strain_info_geo <- indep_strain_info %>%
  dplyr::mutate(geo = ifelse(isotype %in% hw_strain, "Hawaii", 
                             ifelse(isotype %in% ca_strain, "Central America", 
                                    ifelse(isotype %in% au_strain, "Australia", 
                                           ifelse(isotype %in% af_strain, "Africa", 
                                                  ifelse(isotype %in% car_strain, "Caribbean",
                                                                ifelse(isotype %in% sa_strain, "South America", 
                                                                              ifelse(isotype %in% tw_strain, "Taiwan", 
                                                                                     ifelse(isotype %in% mal_strain,"Indonesia",
                                                                                            ifelse(isotype %in% mic_strain,"Micronesia","Unknown"))))))))))
unique(indep_strain_info_geo$geo)


######## Look into the differences #######
etc_strain_location <- indep_strain_info %>% 
  dplyr::filter(strain %in% etc_strain)
etc_strain_location

## Geographic distribution of all strains 
world <- map_data('world')
world <- world[world$region != "Antarctica",]

plot_world <-ggplot2::ggplot()+ geom_map(data=world, map=world,
                                         aes(x=long, y=lat, map_id=region),
                                         color="gray95", fill="gray80", linewidth=0.05)+
  geom_point(data = indep_strain_info_geo, aes(x= long, y=lat, color = geo), shape =16, size = 1.2, alpha = 1.0) +
  scale_color_manual(values = geo.colours) +
  lims(x=c(-180,191),y=c(-58,84)) +
  theme_void()+
  labs(color = NULL)+
  theme(legend.position = "none")+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 5))+
  theme(legend.key.size = unit(0.4, "cm"))

### Pie charts of the strains
geo_freq <- indep_strain_info_geo %>%
  dplyr::group_by(geo) %>%
  dplyr::summarize(frequency = n()) %>%
  arrange(desc(frequency))

plot_freq <- ggplot(geo_freq, aes(x = "", y = frequency, fill = geo)) +
 geom_bar(stat = "identity", width = 1,) +
 coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
      axis.title = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none") +
 geom_text(aes (x = 1.3, label = frequency),
              color = "black",
           position = position_stack(vjust = 0.5),
           size = 3,
         hjust = 0.5) +
scale_fill_manual(values = geo.colours)+
 scale_color_manual(values = geo.colours)

plot_freq<- plot_freq +
 geom_text(aes(x = 1.65, label = geo),
           color = "black",
                      position = position_stack(vjust = 0.5),
                       hjust = 0.5, vjust = 0.5,
                      size = 3)

#### plot freq ###### 
tmp_vec <- seq(1,nrow(geo_freq),1)
geo_freq2 <- geo_freq %>% 
  dplyr::arrange(frequency)
geo_freq2$group <-tmp_vec  

plot_freq2 <- ggplot2::ggplot(geo_freq2) +
  geom_rect(aes(xmin=0,xmax=frequency,ymin=group-0.2,ymax=group+0.2,fill=geo)) + 
  geom_text(aes(x=1,y=group+0.42,label=geo),
            fontface = "bold", hjust=0, size = 2.5) +
  geom_text(aes(x=frequency+20,y=group,label=as.character(frequency)),size = 2.5) + 
  scale_fill_manual(values=geo.colours) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none')

# Assemble the 2 plots
p <- ggpubr::ggarrange(ggpubr::ggarrange(plot_freq2,
                                         plot_world,
                                         ncol = 2, 
                                         labels = c("a","b"),
                                         widths = c(0.25, 0.75),
                                         label.x = c(-0.02,0)))

saveRDS(p, file = "../../processed_data/geo_info/strain_map.rds")


