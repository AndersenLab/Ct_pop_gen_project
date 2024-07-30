# install.packages("raster")
# install.packages("sf")

rm(list = ls())

library(glue)
library(rstudioapi)
library(dplyr)

library(raster)
library(sf)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))
getwd()

hfi_dir <- "/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/data/tif"
hfi_files <- list.files(hfi_dir, pattern = "hfp\\d{4}\\.tif$", full.names = TRUE)
hfi_rasters <- lapply(hfi_files, raster)
hfi_stack <- stack(hfi_rasters)
hfi_mean <- calc(hfi_stack, fun = mean, na.rm = TRUE)
output_file <- "/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/tif_results/hfi_mean_2000_2020.tif"
writeRaster(hfi_mean, filename = output_file, format = "GTiff", overwrite = TRUE)






################### Start mapping HFI of each sampling point ############# 




#### C. tropicalis ####

# Read the raster file of calculated mean value
hfi_mean <- raster("/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/tif_results/hfi_mean_2000_2020.tif")
plot(hfi_mean)
# Read sample data
tro_geo_data_raw <- read.csv("../processed_data/Geo_info/indep_isotype_info_geo.csv")

# Select and rename columns
tro_geo_data <- tro_geo_data_raw %>%
  dplyr::select(isotype, long, lat) %>%
  rename(Longitude = long, Latitude = lat, SampleID = isotype)

# Convert sample data to sf object
tro_samples_sf <- st_as_sf(tro_geo_data, coords = c("Longitude", "Latitude"), crs = 4326)

# Check the CRS of HFI raster data
hfi_crs <- crs(hfi_mean)
print(hfi_crs)

# If the CRS of sample data and HFI raster data are inconsistent, then transform the CRS of sample data
if (!st_crs(tro_samples_sf) == hfi_crs) {
  tro_samples_sf <- st_transform(tro_samples_sf, crs = hfi_crs)
}

# Extract the mean HFI values for sample locations
tro_geo_data$HFI <- extract(hfi_mean, as(tro_samples_sf, "Spatial"))

# have a look at HFI mean
print(tro_geo_data)

# save csv
write.csv(geo_data, "/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/tif_results/C_tro_samples_with_hfi_mean.csv", row.names = FALSE)








#### C. elegans ####


# Read sample data
ele_geo_data_raw <- read.csv("../processed_data/Geo_info/20231213_c_elegans_strain_data.csv")

# Select and rename columns
ele_geo_data <- ele_geo_data_raw %>%
  dplyr::filter(strain==isotype) %>% 
  dplyr::select(isotype, longitude, latitude) %>%  
  na.omit() %>%
  rename(Longitude = longitude, Latitude = latitude, SampleID = isotype)


# Convert sample data to sf object
ele_samples_sf <- st_as_sf(ele_geo_data, coords = c("Longitude", "Latitude"), crs = 4326)


# If the CRS of sample data and HFI raster data are inconsistent, then transform the CRS of sample data
if (!st_crs(ele_samples_sf) == hfi_crs) {
  ele_samples_sf <- st_transform(ele_samples_sf, crs = hfi_crs)
}

# Extract the mean HFI values for sample locations
ele_geo_data$HFI <- extract(hfi_mean, as(ele_samples_sf, "Spatial"))

# have a look at HFI mean
print(ele_geo_data)

# save csv
write.csv(ele_geo_data, "/home/bwang97/vast-eande106/projects/Bowen/PopGen_Tro_Project/processed_data/tif_results/C_ele_samples_with_hfi_mean.csv", row.names = FALSE)





