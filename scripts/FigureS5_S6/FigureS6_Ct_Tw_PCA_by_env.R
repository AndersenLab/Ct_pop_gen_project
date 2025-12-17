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




# load utilities. Findboxes
source("../utilities.R")



### redefine lineage_colors for Taiwan
lineage_colors
library(RColorBrewer)

display.brewer.all()
display.brewer.pal(9,"Set1")

lineage_colors <- c(
  Tw1 = "#E41A1C",
  Tw2 = "#FF7F00",
  Tw3 = "#FFFF33",
  Tw4 = "#4DAF4A",
  Tw5 = "#40E0D0",
  Tw6 = "#377EB8",
  Tw7 ="#984EA3"
)



lineage<-read.csv("../../processed_data/geo_info/Ct_lineage_all.csv")

cdf9_7 <- readRDS(file = "../../processed_data/Taiwan/cdf9_7no.rds")


cdf9_7_revised<-cdf9_7 %>%
  rename(`mean annual rainfall`=annual_rainfall,
         `mean annual temperature °C`=annual_temperature,
         `Area Solar Radiation WH/m2` = ASR
         )

library(tibble)
Ct_pca_con_df<- cdf9_7_revised %>% 
  column_to_rownames("strain") %>% 
  dplyr::select(-isotype,
                -geo, 
                -lineage) %>% 
  as.data.frame() %>% 
  mutate_all(as.numeric)



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
reduced <- PC_by_env_cors[7:ncol(PC_by_env_cors), 1:6]






PC_by_env_cors2_flattened_plot <- PC_by_env_cors2_flattened %>%
  mutate(label = ifelse(sig, paste0(round(cor, 2), "*"), round(cor, 2)))


pc_heatmap <- ggplot(PC_by_env_cors2_flattened_plot, aes(x = PC, y = env_par, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), fontface = ifelse(PC_by_env_cors2_flattened$sig, "bold", "plain"), size = 4) +
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

# cowplot::ggsave2(pc_heatmap, filename = "Ct_PCA_heatmap.pdf",width = 7.5, heigh = 4, units = "in")











stat_df <- cdf9_7_revised %>%
  dplyr::select(-geo,
                ) %>%
  dplyr::filter(!is.na(lineage)) %>%
  dplyr::mutate(lineage = factor(lineage, levels = c("Tw1", "Tw2", "Tw3", "Tw4", "Tw5","Tw6","Tw7"))) %>%
  tidyr::gather(env_par, value, elevation:`mean annual temperature °C`) %>%
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
    
    # add plot to plot list
    plot_list[[i]] <- p
    
    # rename plot
    names(plot_list)[i] <- glue::glue("{unique(df$var)}_sigplot")
  }
  return(plot_list)
}

# plot all of these comparisions
sig_plots <- plot_dunn(Dunn_list)
sig_plots_all <- cowplot::plot_grid(plotlist = sig_plots, ncol = 3)
sig_plots_all
# cowplot::ggsave2(sig_plots_all, filename = "Ct_PCA_by_cont_env_sig_plots.pdf",width = 10, heigh = 7, units = "in")



























#########################################
# Continuous plots
#########################################
# Plot data as box plots for each environmental parameter
plot_temp <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual temperature °C")) +
  scale_fill_manual(values=lineage_colors) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.8,linewidth=0.1) +
  theme_bw() +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "group", x = "", y = "Mean temp (°C)") #+
plot_temp



plot_rain <- ggplot(data = stat_df %>% dplyr::filter(env_par == "mean annual rainfall")) +
  scale_fill_manual(values=lineage_colors) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.8,linewidth=0.1) +
  theme_bw() +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "group", x = "", y = "Mean rainfall")
plot_rain







plot_asr <- ggplot(data = stat_df %>% dplyr::filter(env_par == "Area Solar Radiation WH/m2")) +
  scale_fill_manual(values=lineage_colors) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.8,linewidth=0.1) +
  theme_bw() +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "group", x = "", y = "ASR") 
plot_asr





plot_elevation <- ggplot(data = stat_df %>% dplyr::filter(env_par == "elevation")) +
  scale_fill_manual(values=lineage_colors) +
  stat_summary(aes(x = lineage, y = value, fill=lineage),color = "black",
               fun = mean, geom = "bar", alpha = 0.8,linewidth=0.1) +
  theme_bw() +
  geom_beeswarm(aes(x = lineage, y = value, fill = lineage), size = 1.5, stroke = 0.25, shape = 21) +
  theme(text = element_text(size = 10), panel.grid = element_blank(), strip.background = element_blank(), axis.text = element_text(size = 9, color = "black"), plot.margin = unit(c(0.25,0.25,0.25,0.25), "lines"), legend.position="none") + 
  labs(fill = "Species", x = "", y = "Elevation")
plot_elevation




all <- cowplot::plot_grid(plot_temp,plot_rain, plot_asr,  plot_elevation,
                          labels = c("B", "C", "D", "E"), 
                          label_x = c(-0.015, 0, -0.015, 0),
                          vjust = 1.5, 
                          ncol = 2, 
                          align = "vh", 
                          axis = "tl")

# cowplot::ggsave2(all, filename = "Ct_PCA_by_cont_env.pdf",width = 7.5, heigh = 5, units = "in")







merge_supp<- cowplot::plot_grid(pc_heatmap,all,
                                labels = c("A", ""),
                                label_x = c(-0.0075, 0),
                                ncol = 1, 
                                rel_heights = c(1.5, 2)
                                )
merge_supp

cowplot::ggsave2(merge_supp, filename = "../../plots/FigureS6_Taiwan_env_all.pdf",
                 width = 7.5, heigh = 7.5, units = "in")


  

