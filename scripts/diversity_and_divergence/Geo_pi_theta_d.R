rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

calculate_diversity<-function(custom_region,diversity_type){

  file_path <- c('../../processed_data/pi_theta_d_geo/region/chromosome_windows_diversity.csv')
  new_file_path <- gsub("region", custom_region, file_path)
  Geo_stat_raw <- read.csv(new_file_path)
  
diversity <- Geo_stat_raw %>%
  dplyr::filter(stat_type == diversity_type) %>%
  dplyr::select(stat) %>% 
  na.omit()
  
output<-mean(diversity$stat)
return(output)
    
}

Africa_theta<-calculate_diversity("Africa","theta")
Australia_theta<-calculate_diversity("Australia","theta")
Caribbean_theta<-calculate_diversity("Caribbean","theta")
Central_America_theta<-calculate_diversity("Central_America","theta")
Hawaii_theta<-calculate_diversity("Hawaii","theta")
South_America_theta<-calculate_diversity("South_America","theta")
Taiwan_theta<-calculate_diversity("Taiwan","theta")
Indonesia_theta<-calculate_diversity("Indonesia","theta")
Micronesia_theta<-calculate_diversity("Micronesia","theta")

Africa_pi<-calculate_diversity("Africa","pi")
Australia_pi<-calculate_diversity("Australia","pi")
Caribbean_pi<-calculate_diversity("Caribbean","pi")
Central_America_pi<-calculate_diversity("Central_America","pi")
Hawaii_pi<-calculate_diversity("Hawaii","pi")
South_America_pi<-calculate_diversity("South_America","pi")
Taiwan_pi<-calculate_diversity("Taiwan","pi")
Indonesia_pi<-calculate_diversity("Indonesia","pi")
Micronesia_pi<-calculate_diversity("Micronesia","pi")

Africa_d<-calculate_diversity("Africa","d")
Australia_d<-calculate_diversity("Australia","d")
Caribbean_d<-calculate_diversity("Caribbean","d")
Central_America_d<-calculate_diversity("Central_America","d")
Hawaii_d<-calculate_diversity("Hawaii","d")
South_America_d<-calculate_diversity("South_America","d")
Taiwan_d<-calculate_diversity("Taiwan","d")
Indonesia_d<-calculate_diversity("Indonesia","d")
Micronesia_d<-calculate_diversity("Micronesia","d")

species_wide_stats<-read.csv(file = "../../processed_data/pi_theta_d/chromosome_windows_diversity.csv")

species_pi <- species_wide_stats %>%
  dplyr::filter(stat_type == "pi") %>%
  dplyr::select(stat) %>% 
  na.omit()
species_pi<-mean(species_pi$stat)

species_theta <- species_wide_stats %>%
  dplyr::filter(stat_type == "theta") %>%
  dplyr::select(stat) %>% 
  na.omit()
species_theta<-mean(species_theta$stat)

species_d <- species_wide_stats %>%
  dplyr::filter(stat_type == "d") %>%
  dplyr::select(stat) %>% 
  na.omit()
species_d<-mean(species_d$stat)

table_geo_p_theta_d <- data.frame(
  region = rep(c("All","Africa", "Australia", "Caribbean", "Central_America", "Hawaii", "Indonesia", "Micronesia","South_America", "Taiwan"), each = 3),
  stat = rep(c("theta", "pi","d"), times = 10),
  values = c(species_theta, species_pi, species_d, 
             Africa_theta, Africa_pi, Africa_d,
             Australia_theta, Australia_pi, Australia_d, 
             Caribbean_theta, Caribbean_pi, Caribbean_d,
             Central_America_theta, Central_America_pi, Central_America_d,
             Hawaii_theta, Hawaii_pi, Hawaii_d,
             Indonesia_theta, Indonesia_pi, Indonesia_d,
             Micronesia_theta, Micronesia_pi, Micronesia_d,
             South_America_theta, South_America_pi, South_America_d,
             Taiwan_theta, Taiwan_pi, Taiwan_d)
)

table_geo_p_theta_d$values<-round(table_geo_p_theta_d$values,6)

table_geo_p_theta_d<-table_geo_p_theta_d %>% 
  dplyr::arrange(stat)

geo_diversity_result <- table_geo_p_theta_d %>%
  tidyr::pivot_wider(names_from = stat, values_from = values) %>%
  dplyr::filter(!is.na(pi) & !is.na(theta)& !is.na(d)) %>%
  dplyr::select(region, pi, theta, d) %>% 
  dplyr::rename(Region=region,"\u03C0 value"=pi,"\u03B8 value"=theta, "Tajima's D"=d) %>% 
  dplyr::mutate(Region = ifelse(Region == "Indonesia",
                         "Indonesia",
                         Region))

geo_freq <- read.csv(file = "../../processed_data/geo_info/Ct_isotype_geo_freq.csv")
geo_freq_col<-geo_freq %>%
  dplyr::mutate(geo = ifelse(geo == "South America", "South_America", geo)) %>% 
  dplyr::mutate(geo = ifelse(geo == "Central America", "Central_America", geo))

geo_freq_col<-rbind(geo_freq_col, c("All", sum(geo_freq_col$freq))) 
geo_freq_col$geo <- factor(geo_freq_col$geo, 
                           levels = c("All", "Africa", "Australia", "Caribbean",
                                      "Central_America", "Hawaii", "South_America", 
                                      "Indonesia","Micronesia", "Taiwan"))
geo_freq_col<-geo_freq_col %>% 
  dplyr::arrange(geo) %>% 
  dplyr::rename(Region=geo,"Number of strains"=frequency)

all_result<-dplyr::left_join(geo_diversity_result,
                             geo_freq_col,
                             by = c("Region"))

write.csv(file = "../../tables/TableS4_geo_pi_theta_d.csv",
          all_result,
          quote = FALSE,
          row.names = FALSE)


