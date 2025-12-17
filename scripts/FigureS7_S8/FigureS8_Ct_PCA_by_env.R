rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(ggbeeswarm)

# remove.packages("ggrepel")
# devtools::install_version("ggrepel", version = "0.9.1", repos = "http://cran.us.r-project.org")
library(ggrepel)




source("../utilities.R")


lineage<-read.csv("../../processed_data/geo_info/Ct_lineage_all.csv")

cdf9_7 <- readRDS(file = "../../processed_data/Hawaii/cdf9_7no.rds")


cdf9_7_revised<-cdf9_7 %>% 
  rename(`mean annual rainfall`=staterf_mmann,
         `mean annual soil mositure`=sl_mst_ann,
         `mean annual ambient temperature °C`=tair_ann,
         `mean annual surface temperature °C`=tsurf_ann,
         `mean annual LAI`=lai_ann,
         `in situ substrate temperature °C`=substrate_temp,
         `in situ ambient temperature °C`=ambient_temp,
         `in situ humidity`=ambient_humidity
         ) 

Ct_pca_con_df<- cdf9_7_revised %>% 
  dplyr::select(-isotype,
                -geo, 
                -collection_island,
                -lineage) %>% 
  as.data.frame() %>% 
  mutate_all(as.numeric)


# correlation of eigen values with env parameters
PC_by_env_cors <- round(cor(Ct_pca_con_df, use = "complete.obs", method = "pearson"), 3)
PC_by_env_cors2 <- Hmisc::rcorr(as.matrix(Ct_pca_con_df), type = "pearson")
PC_by_env_cors2_flattened <- flattenCorrMatrix(PC_by_env_cors2$r, PC_by_env_cors2$P) %>%
  dplyr::filter(row %in% c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6") & !(column %in% c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))) %>%
  dplyr::mutate(bf_threshold = 0.05/nrow(.),
                sig = ifelse(p < 0.05, T, F),
                bf_sig = ifelse(p < bf_threshold, T, F)) %>%
  dplyr::select(PC = row, env_par = column, everything()) %>%
  dplyr::arrange(bf_sig)

# plots showing lineages by env_parameters
reduced <- PC_by_env_cors[7:18, 1:6]






PC_by_env_cors2_flattened_plot <- PC_by_env_cors2_flattened %>%
  mutate(label = ifelse(sig, paste0(round(cor, 2), "*"), round(cor, 2)))


pc_heatmap <- ggplot(PC_by_env_cors2_flattened_plot, aes(x = PC, y = env_par, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), 
            fontface = ifelse(PC_by_env_cors2_flattened$sig, "bold", "plain"), 
            size = 3) +
  scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#FF7F00", midpoint = 0,
                       name = "Correlation") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"), 
    axis.text.y = element_text(size = 10,color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
pc_heatmap


# cowplot::ggsave2(pc_heatmap, filename = "Ct_PCA_heatmap.pdf",width = 7.5, heigh = 7.5, units = "in")





####=====================================================================####
#                           PCA groups by env                               #
####======================================================================####


stat_df <- cdf9_7_revised %>%
  dplyr::select(-geo,-collection_island) %>%
  dplyr::filter(!is.na(lineage)) %>%
  dplyr::mutate(lineage = factor(lineage, levels = c("Af", "HC", "Hw1", "Hw2", "Hw3"))) %>%
  tidyr::gather(env_par, value, `mean annual rainfall`:elevation,
                `in situ substrate temperature °C`:`in situ humidity`) %>%
  dplyr::group_by(env_par) %>%
  dplyr::mutate(KM_pvalue = kruskal.test(value ~ lineage)[[3]]) %>%
  dplyr::ungroup()

# Kruskal–Wallis test then multiple comparisions test using Dunn's test with pvalues adjusted with Bonferroni method
options(scipen = 999)
Dunn_list <- list()

for (e in 1:length(unique(stat_df$env_par))){
  KM_df <- stat_df %>%
    dplyr::filter(env_par == (unique(stat_df$env_par)[e]))
  
  D_test <- FSA::dunnTest(KM_df$value ~ KM_df$lineage, method = "bonferroni") 
  Dunn_list[[unique(stat_df$env_par)[e]]] <- D_test
}

Dunn_list

#==============================================#
# Plot Dunn_list comparisions w/ function
#==============================================#
plot_dunn <- function(Dunn_list){
  # make plot list
  plot_list <- NULL
  
  # loop through all elements of the dunns list
  for(i in 1:length(names(Dunn_list))){
    
    # get sig value data frame
    df <-  Dunn_list[[i]][[2]] %>%
      dplyr::mutate(var = names(Dunn_list[i])) %>%
      dplyr::mutate(comp = as.character(Comparison)) %>%
      tidyr::separate(comp, into = c("row", "col"), sep = " - ") %>%
      dplyr::arrange(row, col) %>%
      dplyr::mutate(color = case_when(P.adj < 0.05 ~ "yes",
                                      P.adj >= 0.05 ~ "no"))
    # set colors
    sig_color <- c("no" = "#BF0000", "yes" = "#91BAD6")
    
    # plot comparisions
    p <- ggplot(df, aes(row, col)) +
      geom_tile(aes(fill = color)) + 
      geom_text(aes(label = round(P.adj, 3))) +
      scale_fill_manual(values = sig_color) +
      theme_bw() +
      labs(x = "", y = "", fill = "BF.sig", title = glue::glue("{unique(df$var)}"))
    
    plot_list[[i]] <- p
    
    names(plot_list)[i] <- glue::glue("{unique(df$var)}_sigplot")
  }
  return(plot_list)
}

# plot all of these comparisions
sig_plots <- plot_dunn(Dunn_list)
sig_plots_all <- cowplot::plot_grid(plotlist = sig_plots, ncol = 3)
sig_plots_all
# cowplot::ggsave2(sig_plots_all, filename = "Ct_PCA_by_cont_env_sig_plots.pdf",width = 10, heigh = 10, units = "in")
# cowplot::ggsave2(sig_plots_all, filename = "PCA_by_cont_env_sig_plots.png",width = 338, heigh = 338, units = "mm")



























#########################################
# Continuous plots
#########################################
# Plot data as box plots for each environmental parameter
plot_atemp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "in situ ambient temperature °C")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "group", x = "", y = "amb temp (°C)")
plot_atemp

plot_stemp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "in situ substrate temperature °C")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "sub temp (°C)")
plot_stemp

plot_ahum <- ggplot(data = stat_df %>% dplyr::filter(env_par == "in situ humidity")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "amb humidity (%)")  #+
plot_ahum

plot_elev <- ggplot(data = stat_df %>% dplyr::filter(env_par == "elevation")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "group", x = "", y = "elevation (m)") #+
plot_elev

plot_rain <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual rainfall")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "mean ann precip (m)") #+
plot_rain

plot_soil_moisture <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual soil mositure")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "mean ann available\nsoil moisture") 
plot_soil_moisture


plot_tair_ann <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual ambient temperature °C")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "mean ann air temp (°C)") 
plot_tair_ann

plot_tsurf_ann <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual surface temperature °C")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "mean ann surf temp (°C)")
plot_tsurf_ann


plot_leaf_area_index_ann <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual LAI")) +
  scale_fill_manual(values=lineage_colors) +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.2,linewidth=0.1) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "mean ann LAI")
plot_leaf_area_index_ann

all <- cowplot::plot_grid(plot_elev, plot_atemp,  plot_stemp, plot_ahum, plot_tair_ann, plot_tsurf_ann, plot_rain, plot_soil_moisture, plot_leaf_area_index_ann,
                          labels = c("B", "C", "D", "E", "F", "G", "H", "I", "J"),
                          label_size = 12, 
                          label_x = c(-0.025, 0, 0,-0.025, 0, 0,-0.025, 0, 0 ),
                          vjust = 1, ncol = 3,
                          align = "vh", axis = "tl")
all

# cowplot::ggsave2(all, filename = "Ct_PCA_by_cont_env.pdf",width = 7.5, heigh = 7.5, units = "in")
#cowplot::ggsave2(all, filename = "plots/Figure_5_PCA_by_cont_env.png",width = 169, heigh = 169, units = "mm")


library(cowplot)

merge_supp <- ggdraw() +
  draw_plot(
    cowplot::plot_grid(
      pc_heatmap, all,
      labels = c("A", ""),
      label_x = c(-0.015, 0),
      label_y = c(1.01, 1),
      ncol = 1,
      rel_heights = c(1, 2)
    ),
    x = 0.01,
    y = 0,
    width = 0.99,
    height = 1
  )

merge_supp

cowplot::ggsave2(merge_supp, filename = "../../plots/FigureS8_Hawaii_env_all.pdf",
                 width = 7.5, heigh = 7.5, units = "in")

