# Calculate Ne and outcrossing rate 



#### calculate Ne ####
rm(list = ls())


# 1. write a function to calculate Ne
calculate_Ne<-function(input_file){
# calculate mean_theta
# inputs were generated from 1. vcf_to_zarr* files; then 2. zarr_to_pi_theta_d* files
raw_input<-read.csv(input_file)

library(dplyr)
diversity <- raw_input %>%
  dplyr::filter(stat_type == "theta") %>%
  dplyr::select(stat) 

mean_theta<-mean(diversity$stat)

# population size
# calculate Ne

u <- 2.3e-9  # we use the number of C.elegans # 2022CrombieMolEcol and Teterina et al. 2023 # https://github.com/phillips-lab/CR_CE_popgen/blob/9884624ffb585c0132d29f63cab79f55ec3b720e/ancestral/9_run_Relate.srun#L78
Ne <- mean_theta/(2*4*u) # same method, mutation rate should be rescaled by 0.5. (See 2023 Teterina PLOS Genetics paper)
Ne

return(Ne)
}


# 2. calculate Ne for each type

# calculate Ne for C. tropicalis
Ne_all<- calculate_Ne (input_file = "../processed_data/pi_theta_d/chromosome_windows_diversity.csv")
Ne_af<- calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/Africa/chromosome_windows_diversity.csv")
Ne_au<- calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/Australia/chromosome_windows_diversity.csv")
Ne_car<-calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/Caribbean/chromosome_windows_diversity.csv")
Ne_ca<- calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/Central_America/chromosome_windows_diversity.csv")
Ne_hw<- calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/Hawaii/chromosome_windows_diversity.csv")
Ne_sa<- calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/South_America/chromosome_windows_diversity.csv")
Ne_tw<- calculate_Ne (input_file = "../processed_data/pi_theta_d/geo_pi_theta_d/Taiwan/chromosome_windows_diversity.csv")


# calculate population size of C. elegans all 611 isotype reference strains 
# from CaeNDR database (release 20231213)
Ne_all_Ce<- calculate_Ne (input_file = "../processed_data/Ce_diversity/chromosome_windows_diversity.csv")

# calculate population size of C. elegans 111 isotype reference strains 
# collected on Big Island, Hawaii, from CaeNDR database (release 20231213)
Ne_big_island_Ce<- calculate_Ne (input_file ="../processed_data/Ce_big_island/Ce_pi_theta_d_big_island/chromosome_windows_diversity.csv")

# calculate population size of C. elegans 203 isotype reference strains 
# Hawaii, from CaeNDR database (release 20231213)
Ne_hw_Ce<-calculate_Ne(input_file = "../processed_data/Ce_Hawaii/Ce_pi_theta_d_Hawaii/chromosome_windows_diversity.csv")



### export results
table_Ne_resultes <- data.frame(
  Samples = c("All", "Africa", "Australia","Caribbean","Central America",
             "Hawaii","South America","Taiwan","All_Ce","Big_island_Ce","Hawaii_Ce"),
  Ne_value = c(Ne_all, Ne_af, Ne_au,Ne_car,Ne_ca,
            Ne_hw,Ne_sa,Ne_tw,Ne_all_Ce,Ne_big_island_Ce,Ne_hw_Ce)
)


write.csv(file = "../tables/table_Ne_resultes.csv",
          table_Ne_resultes,
          quote = FALSE,
          row.names = FALSE)


#### end of calculate Ne ####






##################################################################








# 
# #### calculate outcrossing rate ####
# # 
# # 
# tmp<-mean(sapply(LD_all_summary$Mean, function(r) { return(1-(1-r)/(r*(1+(8*Ne*0.5))-1))}))*100
# tmp
# #[1] 99.9991
# outcrossing_rate <- 100- tmp
# outcrossing_rate
# # 0.0008955194
# 
# sd(sapply(LD_all_summary$Mean, function(r) { return(1-(1-r)/(r*(1+(8*Ne*0.5))-1))}))*100
# #[1] 0.0009459868




LD_all_summary_raw<-read.table("../processed_data/LD_all_isotype/C_tro_all.summary",
                               sep= ' ')
LD_all_summary_raw_1<-read.table("../processed_data/LD_geo/C_tro_Central America_LD.ld_LD.summary",sep= ' ')
                               

# write a function 
calculate_ocRate <- function(input_file,Ne_value){
  
  summary_data<-read.table(input_file,
               sep= ' ')
  Ne<-Ne_value
  # 1. rearrange the data set
  # remove the wrong header
  summary_data<-summary_data%>% 
    dplyr::filter(V1!="CHR_A")
  # reset the header names
  colnames(summary_data)<-c("chr_A","pos_A","chr_B","pos_B","Mean")
  #change the "23" into "X"
  summary_data$chr_A<-gsub("23","X",summary_data$chr_A)
  summary_data$chr_B<-gsub("23","X",summary_data$chr_B)
  
  
  summary_data$Distance<-(summary_data$pos_B-summary_data$pos_A)/10
  # summary_data[summary_data$chr_A != summary_data$chr_B, ]$Distance<-"NA" #remove intrachromosomal LD
  
  
  # 2. Calculate outcrossing rate
  tmp<-mean(sapply(summary_data$Mean, function(r) { return(1-(1-r)/(r*(1+(8*Ne*0.5))-1))}))*100
  outcrossing_rate <- 100- tmp
  
  sd<-sd(sapply(summary_data$Mean, function(r) { return(1-(1-r)/(r*(1+(8*Ne*0.5))-1))}))*100
  
  # export outcrossing_rate and sd
  results<-c(outcrossing_rate,sd)
  return(results)
}


ocRateSD_all<-calculate_ocRate(input_file = "../processed_data/LD_all_isotype/C_tro_all.summary",
                             Ne_value = Ne_all)
ocRateSD_af<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_Africa_LD.ld_LD.summary",
                               Ne_value = Ne_af)
ocRateSD_au<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_Australia_LD.ld_LD.summary",
                            Ne_value = Ne_au)
ocRateSD_car<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_Caribbean_LD.ld_LD.summary",
                            Ne_value = Ne_car)
ocRateSD_ca<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_Central America_LD.ld_LD.summary",
                            Ne_value = Ne_ca)
ocRateSD_hw<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_Hawaii_LD.ld_LD.summary",
                            Ne_value = Ne_hw)
ocRateSD_sa<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_South America_LD.ld_LD.summary",
                            Ne_value = Ne_sa)
ocRateSD_tw<-calculate_ocRate(input_file = "../processed_data/LD_geo/C_tro_Taiwan_LD.ld_LD.summary",
                              Ne_value = Ne_tw)



# Need LD data for all C.elegans
ocRateSD_Ce_all<-calculate_ocRate(input_file = "../processed_data/Ce_LD_all_isotype/Ce_all.summary",
                              Ne_value = Ne_all_Ce)
# Need LD data and Ne for Hawaii C.elegans 
ocRateSD_Ce_hw<-calculate_ocRate(input_file = "../processed_data/Ce_Hawaii/Ce_LD_Hawaii/Ce_Hawaii.summary",
                                  Ne_value = Ne_hw_Ce)







# generate the data.frame
Region <- c("All", "Africa", "Australia", 
            "Caribbean", "Central_America", "Hawaii", 
            "South_America", "Taiwan","C_elegans_All",
            "C_elegans_Hawaii")
Outcrossing_Rate <- c(ocRateSD_all[1], ocRateSD_af[1], 
                      ocRateSD_au[1], ocRateSD_car[1], 
                      ocRateSD_ca[1], ocRateSD_hw[1], 
                      ocRateSD_sa[1], ocRateSD_tw[1],
                      ocRateSD_Ce_all[1],ocRateSD_Ce_hw[1])
SD <- c(ocRateSD_all[2], ocRateSD_af[2], 
        ocRateSD_au[2], ocRateSD_car[2], 
        ocRateSD_ca[2], ocRateSD_hw[2], 
        ocRateSD_sa[2], ocRateSD_tw[2],
        ocRateSD_Ce_all[2],ocRateSD_Ce_hw[2])

outcrossing_rate_table <- data.frame(Region = Region, Outcrossing_Rate = Outcrossing_Rate, SD = SD)

write.csv(file = "../tables/outcrossing_rate_table.csv",
          outcrossing_rate_table,quote = FALSE,row.names = FALSE)









