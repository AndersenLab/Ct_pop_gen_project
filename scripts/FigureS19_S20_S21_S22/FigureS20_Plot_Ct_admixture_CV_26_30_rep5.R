
rm(list=ls())



#Load necessary packages
# remotes::install_local("pophelper-master.zip")
library(pophelper)
library(tidyverse)
library(ggthemes)

source("../utilities.R")


# My colors 
ancestry.colours <- setNames(
  c("#4b0200","#da000f", "#ff6d93", "#d45700",
    "#563900", "#ffe5ca", "#ffb914", "#ffda90", "#77c000",
    "#01e51b", "#00491e", "#01dea2", "#82fffa", "#00a4b1",
    "#5cb9ff", "#000f2d", "#0141b9", "#9c87ff", "#cbb1ff",
    "#f479ff", "#5f0058", "#ad0041","#ffb4a8"
    
  ),
  c(LETTERS[1:23])
)

extra_cols <- setNames(
  c(
    "#FF1493",
    "#00FF7F",
    "#8B00FF",
    "#FFD700",
    "#1E90FF",
    "black",
    "grey"
  ),
  c("X", "Y", "Z", "AA", "AB","AC","AD")
)

ancestry.colours <- c(ancestry.colours, extra_cols)


isotype_geo_info<-read.csv("../../processed_data/Geo_info/Ct_indep_isotype_info_geo.csv") 

sample_list<-read.table("../../processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt")
library(dplyr)

desired <- as.character(sample_list$V1)

isotype_geo_info_ordered <- isotype_geo_info %>%
  filter(isotype %in% desired) %>%
  mutate(isotype = factor(isotype, levels = desired)) %>%
  arrange(isotype)


sample_names <- isotype_geo_info_ordered$isotype


best_k <- data.frame(K.28 = 26:30)
best_k






# generate K summary plot - Supplemental figure XX
k_values <- best_k$K.28
best_k_value <- as.numeric(sub("K.", "", colnames(best_k)[1]))


############ which replicate (1:10) #####

which_replicate<-c(5)

for (which_replicate in which_replicate) {
  
  admix_plots <- list()
  for(kpops in 1:length(grep(".Q", list.files("../../processed_data/Ct_admixture/"), value = T))){
    K <- as.numeric(sub(".Q", "", strsplit(grep(".Q", list.files("../../processed_data/Ct_admixture/"), value = T)[kpops], split = "\\_")[[1]][3]))
    
    if (!(K %in% k_values)) {
      next  
    }
    
    # load Q files
    qfile_name <- grep(pattern = glue::glue("_{K}_\\d+\\.Q$"), value = T, x = list.files("../../processed_data/Ct_admixture/"))
    qfile <- pophelper::readQ(files = paste0("../../processed_data/Ct_admixture/",qfile_name))[[which_replicate]]
    # add pop names
    names_pool <- c(LETTERS, paste0("A", LETTERS))
    colnames(qfile) <- names_pool[1:K]
    

    qfile <- qfile %>%
      dplyr::mutate(samples = sample_names)
    
    # make long and determin order of plotting
    long_admix_pops <- qfile %>%
      dplyr::mutate(samples = sample_names) %>%
      tidyr::gather(cluster, frac_cluster, -samples) %>%
      dplyr::group_by(samples) %>%
      dplyr::mutate(max_frac = max(frac_cluster)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(cluster, max_frac) %>%
      dplyr::mutate(samples = factor(samples, levels = unique(samples)))
    
    # establish plot order of strains based on anc pop and max fraction
    plot_order <- long_admix_pops %>%
      dplyr::filter(frac_cluster == max_frac) %>%
      dplyr::arrange(cluster, -max_frac) %>%
      dplyr::mutate(samples = factor(samples, levels = unique(samples)))
    
    
    admix_plots[[kpops]] <-long_admix_pops %>%
      dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
      ggplot() +
      geom_bar(stat = "identity", 
               aes(x = ordered_samples, 
                   y = frac_cluster, 
                   fill = cluster)) +
      scale_fill_manual(values = ancestry.colours) +
      labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
      theme_bw() +
      theme(axis.text.x=element_blank(),    
            axis.text.y=element_blank(),
            axis.title.y = element_text(size = 8, angle = 90, vjust = .5),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank(),
            plot.margin = unit(c(0,0,0,0), units = "cm"),
            legend.position = "bottom"
            )+
      guides(fill = guide_legend(nrow = 2, ncol = 15))

    if(!exists("representative_K_strains")){
      representative_K_strains <- dplyr::filter(plot_order, max_frac  > 0.999
      ) %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(sample_n = 1:n()) %>%
        dplyr::top_n(3, sample_n) %>%
        dplyr::mutate(K_size = K)
    } else {
      representative_K_strains <- dplyr::filter(plot_order, frac_cluster > 0.999
      ) %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(sample_n = 1:n()) %>%
        dplyr::top_n(3, sample_n) %>%
        dplyr::mutate(K_size = K) %>%
        dplyr::bind_rows(representative_K_strains, .)
    }
    
    admix_plots <- Filter(Negate(is.null), admix_plots)
    
  }


samples <- isotype_geo_info_ordered[, "isotype"]  # sample names
groups <- isotype_geo_info_ordered[, "geo"]    # group names

sample_colors <- geo.colours[groups]



K_admixture_plot_1 <- admix_plots[[1]] + theme(legend.position = "none")
K_admixture_plot_2 <- admix_plots[[11]] + theme(legend.position = "none")
K_admixture_plot_3 <- admix_plots[[21]] + theme(legend.position = "none")
K_admixture_plot_4 <- admix_plots[[31]] + theme(legend.position = "none")

legend <- cowplot::get_legend(admix_plots[[41]])
K_admixture_plot_5 <- admix_plots[[41]] + theme(legend.position = "none")  


admixture_plots <- cowplot::plot_grid(K_admixture_plot_1,
                                      K_admixture_plot_2,
                                      K_admixture_plot_3,
                                      K_admixture_plot_4,
                                      K_admixture_plot_5,
                                      ncol = 1)
admixture_plots



final_admixture_plots <- cowplot::plot_grid(admixture_plots, 
                                            legend, 
                                            ncol = 1,
                                            rel_heights = c(1, 0.1))
final_admixture_plots

ggsave(final_admixture_plots, 
       filename = paste0("../../plots/FigureS20_raw_Ct_5CVs_admixture_replicate_",which_replicate,".pdf"), height = 7.5, width = 7.5, useDingbats=FALSE)


}









