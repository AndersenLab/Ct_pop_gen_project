rm(list = ls())

library(terra)
library(sf)
library(dplyr)



### geo info 
Ct_df_raw<-read.csv("../../processed_data/geo_info/Ct_indep_isotype_info_geo.csv")

Ct_isotype_list_Taiwan<-Ct_df_raw %>% 
  filter(geo == "Taiwan")

write.csv(Ct_isotype_list_Taiwan,
          "../../processed_data/geo_info/Ct_isotype_list_Taiwan.csv",
          row.names = FALSE)

geo_Taiwan <- "../../processed_data/geo_info/Ct_isotype_list_Taiwan.csv"



library(terra)
library(sf)
library(dplyr)




extract_gis_data <- function(raster_file, geo_list_csv) {
  r <- rast(raster_file) 
  print(r)
  crs(r)
  
  samples <- read.csv(geo_list_csv, stringsAsFactors = FALSE)
  samples$long <- as.numeric(as.character(samples$long))
  samples$lat  <- as.numeric(as.character(samples$lat))
  
  pts <- st_as_sf(samples, coords = c("long", "lat"), 
                  crs = 4326,
                  remove = FALSE)
  
  pts_proj <- st_transform(pts, crs(r))
  
  v <- terra::vect(pts_proj)
  vals <- terra::extract(r, v)
  
  out <- cbind(st_drop_geometry(pts), vals)
  if ("ID" %in% names(out)) out <- out %>% select(-ID)
  return(out)
}







### elevation 
elevation_tif <- "../../data/Taiwan/Taiwan_environmental_dataset-master/GeoTIFF/G1km_TWD97-121_DTM_/G1km_TWD97-121_DTM_ELE.tif"
elevation_results <- extract_gis_data(elevation_tif, geo_Taiwan) %>% 
  rename(elevation=`G1km_TWD97-121_DTM_ELE`)
head(elevation_results)


# Area Solar Radiation, ASR (Watt-hour(WH)/m2)
ASR_tif <- "../../data/Taiwan/Taiwan_environmental_dataset-master/GeoTIFF/G1km_TWD97-121_DTM_/G1km_TWD97-121_DTM_ASR.tif"
ASR_results <- extract_gis_data(ASR_tif, geo_Taiwan)%>% 
  rename(ASR=`G1km_TWD97-121_DTM_ASR`)
head(ASR_results)


# annual_rainfall Bio12
annual_rainfall_tif<-"../../data/Taiwan/Taiwan_environmental_dataset-master/GeoTIFF/G1km_TWD97-121_Climate_2010s_/G1km_TWD97-121_Climate_2010s_Bio12.tif"
annual_rainfall_results <- extract_gis_data(annual_rainfall_tif, geo_Taiwan) %>% 
  rename(annual_rainfall=`G1km_TWD97-121_Climate_2010s_Bio12`)
head(annual_rainfall_results)

# annual_temperature Bio01
annual_temperature_tif<-"../../data/Taiwan/Taiwan_environmental_dataset-master/GeoTIFF/G1km_TWD97-121_Climate_2010s_/G1km_TWD97-121_Climate_2010s_Bio01.tif"
annual_temperature_results <- extract_gis_data(annual_temperature_tif, geo_Taiwan) %>% 
  rename(annual_temperature=`G1km_TWD97-121_Climate_2010s_Bio01`)
head(annual_temperature_results)





merged_data <- elevation_results %>%
  full_join(ASR_results, by = c("strain", "isotype", "geo","lat","long")) %>%
  full_join(annual_rainfall_results, by = c("strain", "isotype", "geo","lat","long")) %>% 
  full_join(annual_temperature_results, by = c("strain", "isotype", "geo","lat","long")) 
nrow(merged_data)
#54


write.table(merged_data,
            "../../processed_data/Taiwan/Ct_TW_isotype_GIS.txt",
            quote = FALSE,
            sep = '\t',
            col.names = TRUE,
            row.names = FALSE)








