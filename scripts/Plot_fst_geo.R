# 

rm(list=ls())


library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)




# load raw data
fst_geo_all_raw<-read.table(file = "../processed_data/fst_dxy/fst/geo/C_tro_fst_geo_fst.txt",
                           sep = '\t',
                           header = TRUE)

##### 1. Fst geo #######


result<-fst_geo_all_raw %>% 
  na.omit() %>% 
  dplyr::select(pop1,pop2,avg_wc_fst) %>% 
  dplyr::group_by(pop1, pop2) %>% 
  dplyr::summarise(Mean_fst = mean(avg_wc_fst)) %>% 
  dplyr::ungroup()

result_tmp1<-result %>% 
  dplyr::mutate(pop1_tmp=pop2) %>% 
  dplyr::mutate(pop2_tmp=pop1) %>% 
  dplyr::select(pop1_tmp,pop2_tmp,Mean_fst) %>% 
  dplyr::rename(pop1=pop1_tmp) %>% 
  dplyr::rename(pop2=pop2_tmp)

result_tmp2<-rbind(result,result_tmp1)

result_tmp3<-result_tmp2 %>% 
  select(pop1) %>% 
  unique() %>% 
  mutate(pop2=pop1) %>% 
  mutate(Mean_fst=NA)
  
result_tmp4<-rbind(result_tmp2,result_tmp3)

# make the matrix
result_wide <-tidyr::spread(result_tmp4, pop2, Mean_fst)
result_wide<-as.data.frame(result_wide)
raw_names<-result_wide$pop1
result_wide<-result_wide %>%
  dplyr::select(!(pop1))
row.names(result_wide)<-raw_names
result_wide[lower.tri(result_wide)]<-NA
result_wide$pop<-row.names(result_wide)


# convert to long format
result_long <- tidyr::pivot_longer(result_wide, cols = -pop, names_to = "region", values_to = "Mean_Fst")

result_long <- na.omit(result_long)

# plot heat map
p_fst_geo<-ggplot(result_long, aes(x = region, y = pop, fill = Mean_Fst)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient(low = "lightyellow", high = "orange") +  
  ggplot2::theme_minimal() +  
  ggplot2::labs(x = "Region", y = "Population") +  
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
  theme(legend.position = c(0.2, 0.7))+
  labs(x = NULL, y = NULL)+
  geom_text(aes(label = round(Mean_Fst, 3)), color = "black", size = 3)  # add text labels

p_fst_geo

# ggsave("../plots/Fst_geo_raw.pdf", plot = p_fst_geo, width = 7.5, height = 5, units = "in")

# ggsave("Fst_geo_raw_ppt.pdf", plot = p_fst_geo, width = 5.1, height = 3.4, units = "in")




####  remove under sampled regions Autralia and Africa
result_long_remove_AA<-result_long %>% 
  dplyr::filter(pop!="Australia" & pop!="Africa") %>%
  dplyr::filter(region!="Australia" & region!="Africa")


# re-plot heat map
p_fst_geo_remove_AA<-ggplot(result_long_remove_AA, aes(x = region, y = pop, fill = Mean_Fst)) +
  ggplot2::geom_tile(color = "white") +  
  ggplot2::scale_fill_gradient(low = "lightyellow", high = "orange") +  
  ggplot2::theme_minimal() +  
  ggplot2::labs(x = "Region", y = "Population") +  
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
  theme(legend.position = c(0.2, 0.7))+
  labs(x = NULL, y = NULL)+
  geom_text(aes(label = round(Mean_Fst, 3)), color = "black", size = 3)  # add text labels

p_fst_geo_remove_AA

ggsave("../plots/Fst_geo_remove_AA_raw.pdf", plot = p_fst_geo_remove_AA, width = 7.5, height = 5, units = "in")
# ggsave("Fst_geo_remove_AA_raw_ppt.pdf", plot = p_fst_geo_remove_AA, width = 5.1, height = 3.4, units = "in")









