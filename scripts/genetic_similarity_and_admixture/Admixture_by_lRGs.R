rm(list=ls())

library(pophelper) # remotes::install_local("pophelper-master.zip")
library(tidyverse)
library(ggthemes)
library(dplyr)
library(ggpubr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("../utilities.R")

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
desired <- as.character(sample_list$V1)

isotype_geo_info_ordered <- isotype_geo_info %>%
  dplyr::filter(isotype %in% desired) %>%
  dplyr::mutate(isotype = factor(isotype, levels = desired)) %>%
  dplyr::arrange(isotype)
sample_names <- isotype_geo_info_ordered$isotype

isotype_by_lineage_raw<-read.csv("../../processed_data/geo_info/Ct_lineage_all.csv") %>% 
  dplyr::rename(isotype=sample) %>% 
  dplyr::rename(Lineage=lineage)

isotype_by_lineage<-isotype_by_lineage_raw %>% 
  dplyr::select(isotype,Lineage) %>% 
  dplyr::filter(isotype %in% sample_names)

###### ###### ###### ###### ###### ###### ###### 
########      Plot figure by RGs     ###### 
#########  plot CV = best_K (=28) figure ###### 
###### ###### ###### ###### ###### ###### ###### 
best_k_value<-28
which_replicate<-c(1:10)
final_plots <- list()
for (which_replicate in which_replicate) {
  
  best_k_qfile<-read.table(paste0("../../processed_data/Ct_admixture_k28/K28_Processed_Ancestry_replicate_",which_replicate,".tsv"),
                           header = TRUE)
  
  best_k_long_admix_pops <- best_k_qfile %>%
    dplyr::mutate(samples = sample_names) %>%
    tidyr::gather(cluster, frac_cluster, -samples) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(max_frac = max(frac_cluster)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>%
    dplyr::left_join(isotype_geo_info_ordered, by = c("samples" = "isotype")) %>% 
    dplyr::left_join(isotype_by_lineage, by = c("samples" = "isotype"))
  
  write.csv(best_k_long_admix_pops,
            paste0("../../processed_data/Ct_admixture_k28/K28_best_k_long_admix_pops_replicate_",which_replicate,".csv")
            )
  
  best_k_plot_order <- best_k_long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  best_k_admix_plots<-best_k_long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = best_k_plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {best_k_value}")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1,
                                     size = 1,
                                     margin = margin(t = 0)),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(size = 8,angle = 90, vjust = .5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_line(linewidth = 0.1),
          axis.ticks.length = unit(0.02, "cm"),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))+
    coord_cartesian(expand = FALSE) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.3, "cm"))+
    guides(fill = guide_legend(nrow = 2)) +
    facet_grid(.~ Lineage, scales = "free_x", space = "free_x")+
    theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                      margin = margin(b = 1, t = 3)),
          strip.background = element_blank(),
          panel.spacing = unit(0.5, "lines")
    )

  pie_best_k_long_admix_pops <- best_k_long_admix_pops %>%
    dplyr::filter(frac_cluster != "0.000010") %>% 
    dplyr::group_by(Lineage, cluster) %>%
    dplyr::summarise(total_frac = sum(frac_cluster), .groups = 'drop') %>%
    dplyr::group_by(Lineage) %>%
    dplyr::mutate(percent = total_frac / sum(total_frac) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster = factor(cluster, levels = sort(unique(cluster))))
  
  stacked_bar <- ggplot(pie_best_k_long_admix_pops, 
                        aes(x = Lineage, y = percent, fill = cluster, order = cluster)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = ancestry.colours) +
    labs(
      x = NULL,
      y = "Percentage (%)",
      fill = "Cluster"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )+
    labs(title = "All isotypes") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 10))
  
  combined_plot <- ggpubr::ggarrange(best_k_admix_plots, stacked_bar,
                                     ncol = 1, nrow = 2, 
                                     heights = c(1, 1))

  rep_strains_best_k_long_admix_pops<-best_k_long_admix_pops %>% 
    dplyr::filter(max_frac > 0.999)
  
  rep_strains_best_k_plot_order <- rep_strains_best_k_long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))

  rep_strains_best_k_admix_plots<-rep_strains_best_k_long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = best_k_plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {best_k_value}")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1,
                                     size = 1,
                                     margin = margin(t = 0)),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(size = 8, angle = 90, vjust = .5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_line(linewidth = 0.1),
          axis.ticks.length = unit(0.02, "cm"),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))+
    coord_cartesian(expand = FALSE) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.3, "cm"))+
    guides(fill = guide_legend(nrow = 2)) +
    facet_grid(.~ Lineage, scales = "free_x", space = "free_x")+
    theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                      margin = margin(b = 1, t = 3)),
          strip.background = element_blank(),
          panel.spacing = unit(0.1, "lines")
    )

  non_admixed_isotype_list<-rep_strains_best_k_plot_order %>% 
    dplyr::select(-frac_cluster,-max_frac)
  
  write.table(non_admixed_isotype_list,
              paste0("../../processed_data/Geo_info/non_admixed_isotype_replicate_",which_replicate,".txt"),
              sep = '\t',
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
 
  rep_strains_prep <- rep_strains_best_k_long_admix_pops %>%
    dplyr::group_by(Lineage, cluster) %>%
    dplyr::summarise(total_frac = sum(frac_cluster), .groups = 'drop') %>%
    dplyr::group_by(Lineage) %>%
    dplyr::mutate(percent = total_frac / sum(total_frac) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster = factor(cluster, levels = levels(pie_best_k_long_admix_pops$cluster)))
  
  stacked_bar_rep_strains <- ggplot(rep_strains_prep, 
                                    aes(x = Lineage, y = percent, fill = cluster, order = cluster)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = ancestry.colours) +
    labs(
      x = NULL,
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none"
    )+ 
    labs(title = "Non-admixed isotypes") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 10))
  
  rep_strains_combined_plot <- ggpubr::ggarrange(rep_strains_best_k_admix_plots, 
                                                 stacked_bar_rep_strains,
                                                 ncol = 1, nrow = 2, 
                                                 heights = c(1, 1)) # 设置上下图的高度比例

  count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
    dplyr::count(Lineage, cluster, name = "n_samples")
  
  p_heatmap_rep_strains<-ggplot(count_df_heatmap_rep_strains, aes(x = cluster, y = Lineage, fill = n_samples)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(
      x = "Cluster",
      y = "Lineage",
      fill = "Number of Samples",
      title = "Samples per Cluster by Lineage"
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank()
    )
  
  count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
    dplyr::count(Lineage, cluster, name = "n_samples") %>%
    dplyr::mutate(
      geo = factor(Lineage, levels = sort(unique(Lineage), decreasing = TRUE))
    )
  
  p_heatmap_rep_strains <- ggplot(count_df_heatmap_rep_strains,
                                  aes(x = cluster, y = Lineage, fill = n_samples)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(n_samples > 0, as.character(n_samples), "")),
              size = 3, color = "black") +
    scale_fill_gradient(low = "lightblue", high = "#1cb6ff") +
    labs(
      x     = NULL,
      y     = NULL,
      fill  = "Number of Samples",
      title = "Non-admixed isotypes per subpoulation by lineage"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank(),
      plot.title.position   = "plot",      
      legend.position = "bottom",
      plot.title = element_text(
        hjust = 0.65,
        face = "bold",
        size = 10 
      )
    )

  second_row <- plot_grid(
    stacked_bar, 
    stacked_bar_rep_strains, 
    p_heatmap_rep_strains,
    ncol         = 3,
    rel_widths   = c(1, 1, 2),  
    labels      = c("b","c","d")
  )
  
  final_plot <- plot_grid(
    best_k_admix_plots,
    second_row,
    ncol        = 1,
    rel_heights = c(1, 1),
    labels      = c("a","") 
  )
  
  if (which_replicate == 5){
  ggsave(paste0("../../figures/raw_FigureS10_admixture_by_RGs_replicate",which_replicate,".pdf"), final_plot, width = 10, height = 8, useDingbats = FALSE)
  }
  
  best_k_qfile<-read.table(paste0("../../processed_data/Ct_admixture_k28/K28_Processed_Ancestry_replicate_",which_replicate,".tsv"),
                           header = TRUE)
  
  best_k_long_admix_pops <- best_k_qfile %>%
    dplyr::mutate(samples = sample_names) %>%
    tidyr::gather(cluster, frac_cluster, -samples) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(max_frac = max(frac_cluster)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples))) %>%
    dplyr::left_join(isotype_geo_info_ordered, by = c("samples" = "isotype")) %>% 
    dplyr::left_join(isotype_by_lineage, by = c("samples" = "isotype"))
  
  best_k_long_admix_pops<-best_k_long_admix_pops %>% 
    dplyr::filter(Lineage %in% c("Tw2", "Tw3", "Tw4", "Tw5", "Tw6", "Mic1","Mic2", "Indo2"
                          ))
  
  best_k_plot_order <- best_k_long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  best_k_admix_plots<-best_k_long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = best_k_plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {best_k_value}")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1, vjust = 1,
                                     size = 1,
                                     margin = margin(t = 0)),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(size = 8,angle = 90, vjust = .5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_line(linewidth = 0.1),
          axis.ticks.length = unit(0.02, "cm"),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))+
    coord_cartesian(expand = FALSE) +
    theme(legend.position = "right",
          legend.key.size = unit(0.1, "cm"),
          legend.text = element_text(size = 4, margin = margin(l = 0, r = 1)), 
          legend.title = element_text(size = 6),
          legend.key.height = unit(0.15, "cm"), 
          legend.key.width  = unit(0.15, "cm"), 
          legend.spacing.x = unit(0.003, "cm"),  
          legend.spacing.y = unit(0, "cm"), 
          legend.margin = margin(0,0,0,0)   )+   
    
    guides(fill = guide_legend(ncol = 2)) +
    facet_grid(.~ Lineage, scales = "free_x", space = "free_x")+
    theme(strip.text.x = element_text(angle = 90, size = 6, face = "bold",
                                      margin = margin(b = 1, t = 3)),
          strip.background = element_blank(),
          panel.spacing = unit(0.1, "lines")
    )

  rep_strains_best_k_long_admix_pops<-best_k_long_admix_pops %>%
    dplyr::filter(max_frac > 0.999)
  
  rep_strains_best_k_plot_order <- rep_strains_best_k_long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  non_admixed_isotype_list<-rep_strains_best_k_plot_order %>% 
    dplyr::select(-frac_cluster,-max_frac)

  count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
    dplyr::count(Lineage, cluster, name = "n_samples")

  count_df_heatmap_rep_strains <- rep_strains_best_k_plot_order %>%
    dplyr::count(Lineage, cluster, name = "n_samples") %>%
    dplyr::mutate(
      geo = factor(Lineage, levels = sort(unique(Lineage), decreasing = TRUE))
    )
  
  p_heatmap_rep_strains <- ggplot(count_df_heatmap_rep_strains,
                                  aes(x = cluster, y = Lineage, fill = n_samples)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(n_samples > 0, as.character(n_samples), "")),
              size = 3, color = "black") +
    scale_fill_gradient(low = "lightblue", high = "#1cb6ff") +
    labs(
      x = NULL,
      y = NULL,
      fill = "n_Samples",
      title = "Non-admixed isotypes per subpop by RGs"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank(),
      plot.title.position   = "plot",      
      legend.position = "right",
      plot.title = element_text(
        hjust = 0.65,
        face = "bold",
        size = 7 
      ),
      legend.direction = "vertical",
      legend.title = element_text(size = 5),
      legend.text  = element_text(size = 4), 
      legend.key.width  = unit(0.2, "cm"), 
      legend.key.height = unit(0.3, "cm"),  
      legend.spacing.x  = unit(0.003, "cm") 
    )

  final_plot <- plot_grid(
    best_k_admix_plots,
    NULL, 
    p_heatmap_rep_strains,
    ncol = 3,
    rel_widths = c(1.4, 0.05, 1) 
  )
  
  final_plots[[ paste0("replicate_", which_replicate) ]] <- final_plot
  
}

combined_only_7_lineage <- plot_grid(plotlist = final_plots, ncol = 2, align = "v")
ggsave("../../figures/raw_FigureS11_admixture_7_RGs.pdf", 
       combined_only_7_lineage, 
       width = 10, height = 7)

