
rm(list=ls())


library(ggplot2)
library(dplyr)
library(geosphere)
# library(readxl)
# library(gsheet)


# library(RColorBrewer)
# display.brewer.all()
# display.brewer.pal(7,"Set2")
# brewer.pal(7, "Set2")
# display.brewer.pal(7,"Set1")
# brewer.pal(7, "Set1")



# Load input data. File was downloaded on April 2nd 2024.
raw_xlsx<-read.csv("../data/20231201_c_tropicalis_strain_data.csv")

# raw_xlsx<-readxl::read_excel("../data/C. tropicalis WI strain info.xlsx")
# # Or load the data from google drive 
# raw_xlsx <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/edit")




#### 1. The Haversine distance between every strain pairwise ####

# filter the data
tmp_task_1_data<-raw_xlsx 
tmp_task_1_data<-raw_xlsx %>%
  dplyr::select(strain,longitude,latitude) 
tmp_task_1_data<- as.data.frame(tmp_task_1_data)
row.names(tmp_task_1_data)<-raw_xlsx$strain
tmp_task_1_data<-tmp_task_1_data %>%
  dplyr::filter(longitude != "NA" | latitude != "NA")
tmp_task_1_data$longitude<-as.numeric(tmp_task_1_data$longitude)
tmp_task_1_data$latitude<-as.numeric(tmp_task_1_data$latitude)
task_1_data<-tmp_task_1_data[,-1]
dim(task_1_data)[1]
# 727 strains


# calculate pairwise Haversine distance
task_1_matrix <- geosphere::distm(task_1_data, fun = distVincentyEllipsoid)
task_1_matrix[upper.tri(task_1_matrix,diag=TRUE)]<-NA
row.names(task_1_matrix)<-row.names(task_1_data)
colnames(task_1_matrix)<-row.names(task_1_data)


# Melt the matrix 
task_1_distance<-task_1_matrix %>% 
  reshape2::melt() %>% 
  dplyr::filter(value!="NA") %>% 
  # dplyr::mutate(paired = stringr::str_c(Var1,Var2,sep = "_") ) %>%
  # dplyr::select(paired,value) %>%
  dplyr::mutate(value=value/1000) %>% # change the unite into kilometer (km)
  dplyr::rename(geo_distance=value)
colnames(task_1_distance)<-c("strain1","strain2","geo_distance")






###### 3. The Haversine distance between strains within an isotype ######

# Filter the strains which have a isotype
task_3_data<-raw_xlsx %>%
  dplyr::select(strain,isotype) %>%
  na.omit() %>% 
  dplyr::group_by(isotype) %>%
  dplyr::filter(n() > 1)
dim(task_3_data)[1]
# 279 strains

colnames(task_3_data)


grouped_data <- task_3_data %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(strain_combinations = list(
    apply(combn(strain, 2), 2, function(x) paste(x, collapse = "_"))
  )) %>%
  tidyr::unnest(strain_combinations) %>% 
  as.data.frame()


# grouped_data <- task_3_data %>%
#   group_by(isotype) %>%
#   summarise(strain_combinations = list(combn(strain, 2, simplify = FALSE))) %>%
#   unnest(strain_combinations)



all_distance<-task_1_distance
all_distance<-all_distance %>% 
  dplyr::mutate(strain_combinations_1= paste(strain1, strain2, sep = "_")) %>% 
  dplyr::mutate(strain_combinations_2= paste(strain2, strain1, sep = "_"))

task_3_distance<-all_distance %>% 
  dplyr::filter(strain_combinations_1 %in% grouped_data$strain_combinations)


# task_3_distance_tmp_1<- 
#   dplyr::left_join(all_distance, grouped_data,
#                    by = c("strain_combinations_1" = "strain_combinations")) %>% 
#   na.omit()

task_3_distance_tmp_2<- 
  dplyr::left_join(all_distance, grouped_data,
                   by = c("strain_combinations_2" = "strain_combinations")) %>% 
  na.omit()

#check if the new file's row number sum is the same as grouped_data's
nrow(grouped_data)
# nrow(task_3_distance_tmp_1)
nrow(task_3_distance_tmp_2)


task_3_distance<-task_3_distance_tmp_2 %>% 
  dplyr::select(geo_distance,strain_combinations_2) %>% 
  rename(strain_combinations=strain_combinations_2)






# plot the histogram
p3_1<-ggplot2::ggplot(task_3_distance, aes(x = geo_distance)) +
  geom_histogram(bins = 20, fill = "#BEBADA", color = "black") +
  stat_bin(bins = 20,geom = "text", aes(label = ..count..), vjust = -0.5, size = 3, color = "black") +
  labs(title = expression("The Vincenty distance between" ~italic("C. tropicalis")~ "strains within an isotype"), x = "Distance (km)", y = "Frequency") +
  theme_classic()

p3_1

# ggsave("task_3_1.pdf", plot = p3_1, width = 7.5, height = 4, units = "in")






# There is only 1 pair of strains seems like outliers: 
# NIC58 and QG3042, NIC58 is the isotype reference strain

# generate a list to have a look at the outlier
outliner<- raw_xlsx %>% 
  dplyr::filter(strain %in% c("NIC58","QG3042"))

write.table(file = "../tables/C_tropicalis_outliners_1.tsv",
            outliner,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')


# If we remove the strains, then re-plot the histogram


task_3_distance_2<-task_3_distance %>%
  dplyr::filter(strain_combinations != "NIC58_QG3042")

# plot the histogram
p3_2<-ggplot2::ggplot(task_3_distance_2, aes(x = geo_distance)) +
  geom_histogram(bins = 20, fill = "#BEBADA", color = "black") +
  stat_bin(bins = 20,geom = "text", aes(label = ..count..), vjust = -0.5, size = 3, color = "black") +
  labs(title = "Remove NIC58 vs. QG3042", x = "Distance (km)", y = "Frequency") +
  theme_classic()

p3_2

# ggsave("task_3_2.pdf", plot = p3_2, width = 7.5, height = 4, units = "in")



# Now have a look at the longest distance strain comparison: 
longest<- raw_xlsx %>% 
  dplyr::filter(strain %in% c("NIC547","NIC580"))






###### task  end ######



p_all <- ggpubr::ggarrange(p3_1,p3_2,
                           # vjust = 20,
                           ncol = 1, nrow = 2,
                           labels = c("a","b"))
p_all

# ggsave("Geographical_distance.pdf", plot = p_all, width = 7.5, height = 7.5, units = "in")





#### merge sunset figure
p3_2_revised<-p3_2 +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))+
  ylim(0,260)+
  theme(plot.title = element_text(size = 11,vjust=-1),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11))
p3_2_revised

p_annotation<-p3_1 + 
  ggplot2::annotation_custom(ggplotGrob(p3_2_revised), 
                             xmin = 1000, xmax = 3000, 
                             ymin = 70, ymax = 280)
p_annotation
ggsave("../plots/Ct_Geographical_distance_annotation.pdf", plot = p_annotation, width = 7.5, height = 4, units = "in")





