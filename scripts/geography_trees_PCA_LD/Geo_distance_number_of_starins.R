rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(geosphere)
library(gsheet)

raw_xlsx<-read.csv("../../data/20250627_c_tropicalis_strain_data.csv")

n_single_isotypes <- raw_xlsx %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(n_strains = n()) %>%
  dplyr::filter(n_strains == 1) %>% 
  nrow()

condordance_groups_raw<-read.table("../../data/isotype_groups.tsv",
                        header = TRUE)

raw_xlsx<-raw_xlsx %>% 
  dplyr::filter(strain %in% condordance_groups_raw$strain)
setdiff(condordance_groups_raw$strain,raw_xlsx$strain)

raw_xlsx<-raw_xlsx %>% 
  dplyr::select(-isotype) %>% 
  dplyr::left_join(condordance_groups_raw, by = c("strain"))

###### The Haversine distance between strains within an isotype ######
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

# calculate pairwise Haversine distance
task_1_matrix <- geosphere::distm(task_1_data, fun = distHaversine)
task_1_matrix[upper.tri(task_1_matrix,diag=TRUE)]<-NA
row.names(task_1_matrix)<-row.names(task_1_data)
colnames(task_1_matrix)<-row.names(task_1_data)

# Melt the matrix 
task_1_distance<-task_1_matrix %>% 
  reshape2::melt() %>% 
  dplyr::filter(value!="NA") %>% 
  dplyr::mutate(value=value/1000) %>%
  dplyr::rename(geo_distance=value)
colnames(task_1_distance)<-c("strain1","strain2","geo_distance")

# Filter the strains which have a isotype
task_3_data<-raw_xlsx %>%
  dplyr::select(strain,isotype) %>%
  na.omit() %>% 
  dplyr::group_by(isotype) %>%
  dplyr::filter(n() > 1)

grouped_data <- task_3_data %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(strain_combinations = list(
  apply(combn(strain, 2), 2, function(x) paste(x, collapse = "_"))
    )) %>%
  tidyr::unnest(strain_combinations) %>% 
  as.data.frame()

all_distance<-task_1_distance
all_distance<-all_distance %>% 
  dplyr::mutate(strain_combinations_1= paste(strain1, strain2, sep = "_")) %>% 
  dplyr::mutate(strain_combinations_2= paste(strain2, strain1, sep = "_"))

task_3_distance<-all_distance %>% 
  dplyr::filter(strain_combinations_1 %in% grouped_data$strain_combinations)

task_3_distance_tmp_2<- 
  dplyr::left_join(all_distance, grouped_data,
                   by = c("strain_combinations_2" = "strain_combinations")) %>% 
  na.omit()

#check if the new file's row number sum is the same as grouped_data's
nrow(grouped_data)
# nrow(task_3_distance_tmp_1)
nrow(task_3_distance_tmp_2)

none_single_isotype<-task_3_distance_tmp_2$isotype %>% 
  unique() %>% 
  length() 
none_single_isotype
# [1] 101
## number of single isotype 622-101 = 521

results_data<-task_3_distance_tmp_2 %>% 
  dplyr::group_by(isotype) %>%
  dplyr::mutate(max_geo_distance = max(geo_distance, na.rm = TRUE)) %>%
  dplyr::ungroup()

isotype_max <- task_3_distance_tmp_2 %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(max_geo_distance = max(geo_distance, na.rm = TRUE))

isotype_max <- isotype_max %>%
  dplyr::mutate(distance_bin = cut(
    max_geo_distance,
    breaks = c(0, 1, 5, 10, 15, 20), 
    labels = c("< 1", "1-5", "5-10", "10-15", "15-20"), 
    right = FALSE
  ))

freq_table <- isotype_max %>%
  dplyr::count(distance_bin)

p_isotype_hist<-ggplot(freq_table, aes(x = distance_bin, y = n)) +
  geom_bar(stat = "identity", fill = "gray40", alpha = 0.7) +
  theme_classic() +
  geom_text(aes(label = n), 
            vjust = -0.3,
            size = 4) + 
  labs(
    x = "Max distance (km)",
    y = "Number of isotypes"
  )+
  ylim(0,92)

saveRDS(p_isotype_hist, file = "../../processed_data/geo_info/p_isotype_hist.rds")

isotype_strain_count <- raw_xlsx %>%
  dplyr::group_by(isotype) %>%
  dplyr::summarise(n_strains = n())

freq_table_strain <- isotype_max %>%
  dplyr::left_join(isotype_strain_count, by = "isotype") %>% 
  dplyr::select(distance_bin,n_strains) %>% 
  dplyr::group_by(distance_bin) %>%
  dplyr::summarise(total_strains = sum(n_strains, na.rm = TRUE))

p_strain_hist<-ggplot(freq_table_strain, aes(x = distance_bin, y = total_strains)) +
  geom_bar(stat = "identity", fill = "gray40", alpha = 0.7) +
  theme_classic() +
  geom_text(aes(label = total_strains), 
            vjust = -0.3,
            size = 4) + 
  labs(
    x = "Max distance (km)",
    y = "Number of strains"
  )+
  ylim(0,225)

saveRDS(p_strain_hist, file = "../../processed_data/geo_info/p_strain_hist.rds")



