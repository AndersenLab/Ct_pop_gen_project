rm(list = ls())

library(DECIPHER)  ## BiocManager::install("DECIPHER")
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize) 
library(readr)
library(grid)


source("../utilities.R")

gtcheck_raw <- readr::read_tsv("../../processed_data/gt_check_no_HDRs/gtcheck_noHDR.tsv",
                           comment = "#", show_col_types = FALSE) 

gt_pairs <- gtcheck_raw %>%
  dplyr::filter(INFO == "DC") %>%
  tidyr::separate(
    `Time required to process one record .. 0.000615 seconds`,
    into = c("i","j","discordance","prob","sites"),
    sep = "\t"
  ) %>%
  dplyr::mutate(
    discordance = as.numeric(discordance),
    sites = as.numeric(sites),
    concordance = 1 - discordance / sites
  ) %>% 
  dplyr::select(-INFO)

gtcheck_output<-gt_pairs %>% 
  dplyr::select(discordance,sites,i,j) %>% 
  dplyr::rename(`different_alleles` = discordance) %>% 
  dplyr::rename(`total genome-wide SNVs` = sites) %>% 
  dplyr::rename(`strain1` = i) %>% 
  dplyr::rename(`strain2` = j)

gtcheck_strain<-gt_pairs %>% 
  dplyr::mutate(concordance = 1-(discordance/sites)) %>% 
  dplyr::select(i,j,concordance)

isotyep_list<-read.csv("../../processed_data/Ct_pruned_VCF_and_PCA/sample_list.txt",
                       header = FALSE)
gtcheck <- gtcheck_strain %>%
  dplyr::filter(i %in% isotyep_list$V1,
         j %in% isotyep_list$V1)

p_histogram<-ggplot()+ geom_histogram(data=gtcheck,aes(x=concordance),
                         bins = 1000)+
  scale_x_continuous(breaks = seq(0.65, 1, by = 0.01))+
  theme_bw()+
  xlab("Genetic similarity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gt_matrix <- gtcheck %>%
  tidyr::pivot_wider(names_from = j, values_from = concordance) %>%
  tibble::column_to_rownames("i") %>%
  as.matrix()

new_row <- matrix(NA, nrow = 1, ncol = ncol(gt_matrix))
rownames(new_row) <- colnames(gt_matrix)[1]
gt_matrix <- rbind(new_row, gt_matrix)

new_col <- matrix(NA, ncol = 1, nrow = nrow(gt_matrix))
colnames(new_col) <- rownames(gt_matrix)[nrow(gt_matrix)]
gt_matrix <- cbind(gt_matrix,new_col)

gt_matrix[upper.tri(gt_matrix,diag = FALSE)] <- t(gt_matrix)[upper.tri(gt_matrix,diag = FALSE)]

sorted_rows <- sort(rownames(gt_matrix))
sorted_cols <- sort(colnames(gt_matrix))

gt_matrix_plot <- gt_matrix[sorted_rows, sorted_cols]
diag(gt_matrix_plot) <- 1

geo_info_raw<-read.csv("../../processed_data/Geo_info/Ct_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw %>%
  tibble::column_to_rownames(var = "isotype") %>% 
  dplyr::select(geo) %>% 
  dplyr::rename(Geo=geo)
geo_info$Geo<-as.factor(geo_info$Geo)
filter_geo.colours <- geo.colours[names(geo.colours) %in% unique(geo_info$Geo)]
col_ordered_geo <- geo_info[sorted_cols, "Geo", drop = FALSE]

col_ha <- HeatmapAnnotation(
  Geo = col_ordered_geo$Geo,
  col = list(Geo = filter_geo.colours),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    Geo = list(title = "Geo",
               title_gp    = gpar(fontsize = 10, fontface = "bold"),
               labels_gp   = gpar(fontsize = 8),
               legend_width  = unit(0.5, "cm"),
               legend_height = unit(2, "cm")
    )
  )
)

min_concordance <- min(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)
max_concordance <- max(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)

p_heatmap <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = colorRamp2(
    c(min_concordance,
      0.924,
      max_concordance),
    c( "#3182bd", "white", "orange")
  ),
  top_annotation       = col_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),   #### internal function of dist() 
  clustering_distance_columns = function(m) as.dist(1 - m),  #### internal function of dist()
  clustering_method_rows = "average",       #### internal function of hclust()
  clustering_method_columns = "average",    #### internal function of hclust()
  show_row_names       = FALSE,
  show_column_names    = FALSE,
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3)
)

p_heatmap_2 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = colorRamp2(
    c(min_concordance,
      0.924,
      max_concordance),
    c( "#3182bd", "white", "orange")
  ),
  top_annotation       = col_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),   #### internal function of dist() 
  clustering_distance_columns = function(m) as.dist(1 - m),  #### internal function of dist()
  clustering_method_rows = "average",       #### internal function of hclust()
  clustering_method_columns = "average",    #### internal function of hclust()
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 0.5),
  column_names_side    = "top",
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3)
)

phylo_vals <- as.vector(gt_matrix_plot)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]  # remove NAs
breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)

phylo_col_fun <- circlize::colorRamp2(
  breaks,
  c("#4575B4", "#87C6C2", "#FFFFE0","#F4D166","#D73027")
)

p_heatmap_3 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = phylo_col_fun,
  top_annotation       = col_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),   #### internal function of dist() 
  clustering_distance_columns = function(m) as.dist(1 - m),  #### internal function of dist()
  clustering_method_rows = "average",       #### internal function of hclust()
  clustering_method_columns = "average",    #### internal function of hclust()
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 0.5),
  column_names_side    = "top",
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3)
)

ht_list <- draw(p_heatmap_3)
col_idx <- column_order(ht_list)
ordered_colnames <- colnames(gt_matrix_plot)[col_idx]

df_lineage_all <- read_csv("../../processed_data/geo_info/Ct_lineage_all.csv")
lineage_vector <- setNames(df_lineage_all$lineage, 
                           df_lineage_all$sample)
lineage_vector <- lineage_vector[colnames(gt_matrix_plot)]
col_ha_lineage <- HeatmapAnnotation(
  Group = lineage_vector,
  col = list(Group = lineage_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(title = "Relatedness\ngroup"),
  annotation_height = unit(6, "mm")
)

col_ha_with_name <- HeatmapAnnotation(
  Geo = col_ordered_geo$Geo,
  col = list(Geo = filter_geo.colours),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_height = unit(6, "mm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    Geo = list(title = "Geo",
               title_gp    = gpar(fontsize = 10, fontface = "bold"),
               labels_gp   = gpar(fontsize = 8),
               legend_width  = unit(0.5, "cm"),
               legend_height = unit(2, "cm")
    )
  )
)

combined_top_ha <- c(col_ha_lineage, col_ha_with_name)

p_heatmap_4 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = phylo_col_fun,
  top_annotation       = combined_top_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),   #### internal function of dist() 
  clustering_distance_columns = function(m) as.dist(1 - m),  #### internal function of dist()
  clustering_method_rows = "average",       #### internal function of hclust()
  clustering_method_columns = "average",    #### internal function of hclust()
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 0.5),
  column_names_side    = "top",
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3)
)

col_ha_lineage <- HeatmapAnnotation(
  `Relatedness group` = lineage_vector,
  col = list(`Relatedness group` = lineage_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  show_legend = FALSE,
  annotation_height = unit(6, "mm")
)

col_ha_with_name <- HeatmapAnnotation(
  Geo = col_ordered_geo$Geo,
  col = list(Geo = filter_geo.colours),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  show_legend = FALSE,
  annotation_height = unit(6, "mm")
)

combined_top_ha <- c(col_ha_lineage, col_ha_with_name)

p_heatmap_5 <- Heatmap(
  matrix = gt_matrix_plot,
  name = "Genetic\nsimilarity",
  col = phylo_col_fun,
  top_annotation = combined_top_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_gp = gpar(lwd = 0.3),
  column_dend_gp = gpar(lwd = 0.3)
)

make_box_label_graphics <- function(cols, labels, fontsize = 8) {
  lapply(seq_along(cols), function(i) {
    force(i)
    function(x, y, w, h) {
      grid.rect(
        x = x,
        y = y + unit(1.6, "mm"),
        width = unit(3.5, "mm"),
        height = unit(3.5, "mm"),
        gp = gpar(fill = cols[i], col = "black", lwd = 0.4)
      )
      ## 文字
      grid.text(
        labels[i],
        x = x,
        y = y - unit(2.8, "mm"),
        gp = gpar(fontsize = fontsize)
      )
    }
  })
}

geo_levels <- c(
  "Caribbean",
  "Hawaii",
  "Australia",
  "Taiwan",
  "South America",
  "Central America",
  "Indonesia",
  "Micronesia",
  "Africa"
)
geo_cols <- unname(filter_geo.colours[geo_levels])

lgd_geo <- Legend(
  title = "Geo              ",
  at = geo_levels,
  labels = rep("", length(geo_levels)),
  graphics = make_box_label_graphics(geo_cols, geo_levels, fontsize = 8),
  nrow = 1,
  by_row = TRUE,
  title_position = "leftcenter",
  gap = unit(21, "mm")
)

grp_levels <- c(
  "HC","Au","Hw1","Tw5", "LAC", 
  "Tw3","Tw2","Indo1","Tw1", "Indo2",
  "Mic2","Mic1", 
  "Tw4", "Tw6", "Hw2","Af",
  "Indo3","Tw7","Hw3"
)
grp_cols <- unname(lineage_colors[grp_levels])

lgd_group <- Legend(
  title = "Relatedness\ngroup",
  at = grp_levels,
  labels = rep("", length(grp_levels)),
  graphics = make_box_label_graphics(grp_cols, grp_levels, fontsize = 8),
  nrow = 1,
  by_row = TRUE,
  title_position = "leftcenter",
  gap = unit(8, "mm")
)

pdf("../../figures/raw_FigureS34_heatmap_nonHDRs.pdf", 
    width = 12, height = 10)
top_annotation_legends <- packLegend(
  lgd_group,
  lgd_geo,
  direction = "vertical",
  gap = unit(5, "mm")
)
draw(
  p_heatmap_5,
  heatmap_legend_side = "left",
  annotation_legend_side = "top",
  annotation_legend_list = list(top_annotation_legends),
  merge_legend = FALSE
)
dev.off()

