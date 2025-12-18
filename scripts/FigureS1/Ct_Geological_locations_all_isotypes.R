rm(list = ls())
getwd()

list.files()

library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(maps)
library(ggpubr)
library(maps)

source("../utilities.R")




raw_data<-read.csv("../../data/20250627_c_tropicalis_strain_data.csv")

indep_isotype_info<- raw_data

sample_list<- read.table("../../processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt")
nrow(sample_list)
#622

indep_isotype_info<-raw_data %>% 
  filter(strain %in% sample_list$V1) %>% 
  rename(lat=latitude) %>% 
  rename(long=longitude) %>% 
  mutate(isotype=strain) %>% 
  select(strain,isotype,lat,long)
  

indep_isotype_info$lat<-as.numeric(indep_isotype_info$lat)
indep_isotype_info$long<-as.numeric(indep_isotype_info$long)
 


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
   dplyr::filter(long > -176  & long < -129 & lat > 4 & lat < 46)  
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
 
 
 
 
 #11. Malay Archipelago
 df_mal_isotype <- indep_isotype_info %>%
   dplyr::filter(long > 96.607511  & long < 130.884854 & lat > -11.145781 & lat < 5.561605)  # Malay Archipelago
 mal_isotype <- as.character(df_mal_isotype$isotype)
 
 
 
 #12. Micronesia
 df_mic_isotype <- indep_isotype_info %>%
   dplyr::filter(long > 157.681622  & long < 158.489117 & lat > 6.683663 & lat < 7.185335) # Micronesia
 mic_isotype <- as.character(df_mic_isotype$isotype)
 
 
 
 
 
 
 
 
 #0. others
 df_etc_isotype <- indep_isotype_info %>%
   dplyr::filter(!isotype %in% c(af_isotype, au_isotype, car_isotype,
                                eu_isotype,hw_isotype, ca_isotype, 
                                nz_isotype,sa_isotype, tw_isotype,
                                #### pac_isotype,
                                mal_isotype, mic_isotype
                                ))
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
                                                                               ifelse(isotype %in% tw_isotype, "Taiwan", 
                                                                                      # ifelse(isotype %in% pac_isotype, "Pacific", 
                                                                                             ifelse(isotype %in% mal_isotype,"Indonesia",
                                                                                                    ifelse(isotype %in% mic_isotype,"Micronesia","Unknown"))))))))))))
 unique(indep_isotype_info_geo$geo)
 

######## Look into the differences #######
etc_isotype_location <- indep_isotype_info %>% 
  dplyr::filter(isotype %in% etc_isotype)
etc_isotype_location
# write.csv(etc_isotype_location,"undefined_region_isotypes.csv", quote = FALSE,row.names = FALSE)
  
  


### output
write.csv(indep_isotype_info_geo,
          "../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv",
          quote = FALSE, row.names = FALSE)


library(RColorBrewer) #palette
display.brewer.all()


display.brewer.pal(8,"Set2")
brewer.pal(8, "Set2")
display.brewer.pal(10, "Set3")
brewer.pal(10, "Set3")
display.brewer.pal(9, "Set1")
brewer.pal(9,"Set1")


world <- map_data('world')
world <- world[world$region != "Antarctica",] # intercourse antarctica

plot_world <-ggplot2::ggplot()+ geom_map(data=world, map=world,
                                         aes(x=long, y=lat, map_id=region),
                                         color="gray95", fill="gray80", linewidth=0.05)+
  geom_point(data = indep_isotype_info_geo, aes(x= long, y=lat, color = geo), shape =16, size = 1.2, alpha = 1.0) +
  scale_color_manual(values = geo.colours) +
  lims(x=c(-180,191),y=c(-58,84)) +
  theme_void()+
  #  theme(text = element_text(size=12)) +
  labs(color = NULL)+
  theme(legend.position = "none")+
  theme(axis.text = element_blank(),    # Conceal Tick Marks
        axis.title = element_blank(),   # Conceal Tick Marks
        legend.text = element_text(size = 5))+
  theme(legend.key.size = unit(0.4, "cm"))

# ggthemes::theme_map()

plot_world


#ggsave("Cb_isotypes_map.pdf",plot = plot_world,width = 7.5, height = 4.37, units = "in")




### Pie charts of the isotypes
geo_freq <- indep_isotype_info_geo %>%
  dplyr::group_by(geo) %>%
  dplyr::summarize(frequency = n()) %>%
  arrange(desc(frequency))


##CHANGE THE SAVING DIRECTORIES TO RF FOR BELOW
write.csv(file = "../../processed_data/geo_info/Ct_isotype_geo_freq.csv", geo_freq,quote = FALSE,
          row.names = FALSE)

plot_freq <- ggplot(geo_freq, aes(x = "", y = frequency, fill = geo)) +
 geom_bar(stat = "identity", width = 1,) +
 coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
      axis.title = element_blank(),
       panel.grid = element_blank(),
       legend.position = "none") +
 # geom_text(aes(label = frequency), position = position_stack(vjust = 0.5)) +
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
                       #angle = 45,
                       hjust = 0.5, vjust = 0.5,
                      size = 3)
                     
plot_freq

####ALT plot freq ###### 

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

plot_freq2


p <- ggpubr::ggarrange(ggpubr::ggarrange(plot_freq2,
                                         plot_world,
                                         ncol = 2, 
                                         labels = c("a","b"),
                                         widths = c(0.25, 0.75),
                                         label.x = c(-0.02, 0)
))
p

saveRDS(p, file = "../../processed_data/assemble_figure_S1/isotype_map.rds")







