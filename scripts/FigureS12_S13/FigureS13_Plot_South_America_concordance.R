rm(list = ls())


library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(dplyr)
source("../utilities.R")




gtcheck_raw <- read.delim("../../processed_data/heatmap_hard_filtered/gtcheck.tsv") 

gtcheck_isotype<-gtcheck_raw %>%
  dplyr::mutate(concordance = 1-(discordance/sites)) %>%
  select(i,j,concordance)

geo_and_lineage_raw<-read.csv("../../processed_data/geo_info/geo_and_lineage.csv")

geo_and_lineage_South<-geo_and_lineage_raw %>%
  filter(geo == "South America")%>% 
  select(isotype,long,lat,geo,lineage)



target_lineage<-c("LAC","Af")
target_geo<-"South America"

target_isotype<-geo_and_lineage_South %>% 
  filter(geo %in% target_geo) %>% 
  filter(lineage %in% target_lineage) %>% 
  select(isotype) %>% 
  pull()






gt_matrix <- gtcheck_isotype %>%
  tidyr::pivot_wider(names_from = j, values_from = concordance) %>%
  tibble::column_to_rownames("i") %>%
  as.matrix()

# add new row
new_row <- matrix(NA, nrow = 1, ncol = ncol(gt_matrix))
rownames(new_row) <- colnames(gt_matrix)[1]
gt_matrix <- rbind(new_row,gt_matrix)

# add new col
new_col <- matrix(NA, ncol = 1, nrow = nrow(gt_matrix))
colnames(new_col) <- rownames(gt_matrix)[nrow(gt_matrix)]
gt_matrix <- cbind(gt_matrix,new_col)

# add upper.tri data
gt_matrix[upper.tri(gt_matrix,diag = FALSE)] <- t(gt_matrix)[upper.tri(gt_matrix,diag = FALSE)]

# sort
sorted_rows <- sort(rownames(gt_matrix))
sorted_cols <- sort(colnames(gt_matrix))

gt_matrix_plot_tmp <- gt_matrix[sorted_rows, sorted_cols]
diag(gt_matrix_plot_tmp ) <- 1

# View(gt_matrix_plot_tmp)



gt_matrix_plot<-gt_matrix_plot_tmp %>% 
  as.data.frame() %>% 
  mutate(isotype = rownames(gt_matrix_plot_tmp)) %>% 
  filter(isotype %in% target_isotype) %>% 
  select(target_isotype)%>%
  select(sort(names(.))) %>% 
  as.matrix()

# sort
sorted_rows <- sort(rownames(gt_matrix_plot))
sorted_cols <- sort(colnames(gt_matrix_plot))







#### threshold
library(tidyr)

gt_matrix_plot_long_tmp <- as.data.frame(gt_matrix_plot) %>%
  mutate(Row = rownames(.)) %>%
  pivot_longer(-Row, names_to = "Col", values_to = "Value") %>% 
  filter(Value!=1)

threshold_histogram<-ggplot()+ 
  geom_histogram(data=gt_matrix_plot_long_tmp,
                 aes(x=Value),
                 bins = 100)+
  scale_x_continuous(breaks = seq(min(gt_matrix_plot_long_tmp$Value, 
                                      na.rm = TRUE),
                                  max(gt_matrix_plot_long_tmp$Value, 
                                      na.rm = TRUE),
                                  by = 0.005))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))+
  theme(panel.grid = element_blank())

threshold_histogram

# ggsave("threshold_histogram.pdf",threshold_histogram,
#        width = 7,height = 3, units = "in")

#### threshold end






########### Annotations #######
### Geo
geo_info_raw<-read.csv("../../processed_data/Geo_info/Ct_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw %>%
  filter(strain %in% rownames(gt_matrix_plot)) %>% 
  tibble::column_to_rownames(var = "isotype") %>% 
  dplyr::select(geo) %>% 
  dplyr::rename(Geo=geo)
geo_info$Geo<-as.factor(geo_info$Geo)


filter_geo.colours <- geo.colours[names(geo.colours) %in% unique(geo_info$Geo)]
col_ordered_geo <- geo_info[sorted_cols, "Geo", drop = FALSE]






###### three breaks #####

#get cc values
phylo_vals <- as.vector(gt_matrix_plot)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]  # remove NAs

# Set the breakpoints by quantile
breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)

#draw color scale
phylo_col_fun <- circlize::colorRamp2(
  breaks,
  #seq(min(phylo_vals, na.rm = TRUE), max(phylo_vals, na.rm = TRUE), length.out = 5),
  c("#4575B4", "#87C6C2", "#FFFFE0","#F4D166","#D73027")
)




# define column annotation: Geo
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




p_heatmap_1 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = phylo_col_fun,
  # col                  = colorRamp2(
  #   c(min_concordance,
  #     0.924,
  #     # mean(c(min_concordance, max_concordance)),
  #     max_concordance),
  #   c( "#3182bd", "white", "orange")
  # ),
  top_annotation       = col_ha,
  # left_annotation      = row_ha,
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



# dev.off()
# pdf(paste0("Ct_isotype_", target_isotype, "_concordance_heatmap.pdf"), width = 7, height = 6)
# 
# draw(p_heatmap_1, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right")
# dev.off()
###### three breaks end #####







##### plot using threshold ### 

# calculate min and max as before
min_concordance <- min(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)
max_concordance <- max(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)

p_heatmap_2 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = colorRamp2(
    c(min_concordance,
      0.996,
      # 0.96,
      max_concordance),
    c("#3182bd", "white", "orange")
  ),
  top_annotation       = col_ha,
  # left_annotation      = row_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),   #### internal function of dist() 
  clustering_distance_columns = function(m) as.dist(1 - m),  #### internal function of dist()
  clustering_method_rows = "average",       #### internal function of hclust()
  clustering_method_columns = "average",    #### internal function of hclust()
  show_row_names       = FALSE,
  show_column_names    = FALSE,
  # show_column_names    = TRUE,
  row_dend_gp          = gpar(lwd = 0.3),
  column_dend_gp       = gpar(lwd = 0.3),
  
  # column_title         = paste0(target_isotype," isotypes"),
  column_title_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_title_rot     = 0

)




# dev.off()
# pdf(paste0("Ct_isotype_", target_isotype, "_concordance_heatmap_2.pdf"), width = 7, height = 6)
# 
# draw(p_heatmap_2, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right")
# dev.off()










###### with lat and long #####

## lat, long, geo annotaion
isotype_anno <- geo_info_raw %>%
  filter(isotype %in% target_isotype) %>%
  left_join((geo_and_lineage_raw %>% select(isotype,lineage)), by = "isotype") %>% 
  tibble::column_to_rownames(var = "isotype") %>%
  dplyr::select(geo, lineage, lat, long,)




# factor region
isotype_anno$geo <- as.factor(isotype_anno$geo)

col_ordered_geo_target  <- isotype_anno[sorted_cols, "geo", drop = FALSE]
col_ordered_lineage_target  <- isotype_anno[sorted_cols, "lineage", drop = FALSE]
col_ordered_lat_target  <- as.numeric(isotype_anno[sorted_cols, "lat"])
col_ordered_long_target <- as.numeric(isotype_anno[sorted_cols, "long"])

# # colours
filter_geo.colours_target <- geo.colours[names(geo.colours) %in% unique(col_ordered_geo_target$geo)]
filter_lineage.colours_target <- lineage_colors[names(lineage_colors) %in% unique(col_ordered_lineage_target$lineage)]






### minimize the range 
lat_col_fun <- colorRamp2(
  c(min(col_ordered_lat_target), 
    mean(c(min(col_ordered_lat_target), max(col_ordered_lat_target)),na.rm = TRUE), 
    max(col_ordered_lat_target)),
  c("#DFF5E1", "#4CAF50", "#1B5E20")
  # c("whitesmoke", "#6baed6", "#08306b")
)

long_col_fun <- colorRamp2(
  c(min(col_ordered_long_target), 
    mean(c(min(col_ordered_long_target),max(col_ordered_long_target)), na.rm = TRUE), 
    max(col_ordered_long_target)),
  c("whitesmoke", "#fa9752", "#67000d")
)





# define column annotation: Geo + Lat + Long
col_ha_target <- HeatmapAnnotation(
  Geo  = col_ordered_geo_target$geo,
  Group = col_ordered_lineage_target$lineage,
  Lat  = col_ordered_lat_target,
  Long = col_ordered_long_target,
  col = list(
    Geo  = filter_geo.colours_target,
    Group = filter_lineage.colours_target,
    Lat  = lat_col_fun,
    Long = long_col_fun
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_legend_param = list(
    Geo = list(
      title = "Geo",
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_width  = unit(0.6, "cm"),
      legend_height = unit(2, "cm")
    ),
    Group = list(
      title = "Group",
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_width  = unit(0.6, "cm"),
      legend_height = unit(2, "cm")
    ),
    Lat = list(
      title = "Latitude",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
      # ,at = pretty(lat_rng, n = 3)
    ),
    Long = list(
      title = "Longitude",
      title_gp  = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
      # ,
      # at = pretty(long_rng, n = 3)
    )
  )
)


table(col_ordered_geo_target$geo)



##### use the same values defined by quantile in the main heatmap 
breaks <-c(0.6865332,0.8853456,0.9068701,
           0.9537771,1.0000000 )
#draw color scale
phylo_col_fun <- circlize::colorRamp2(
  breaks,
  #seq(min(phylo_vals, na.rm = TRUE), max(phylo_vals, na.rm = TRUE), length.out = 5),
  c("#4575B4", "#87C6C2", "#FFFFE0","#F4D166","#D73027")
)



library(grid)

p_heatmap_3 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Genetic\nsimilarity",
  col                  =phylo_col_fun,
  # col                  = colorRamp2(
  #   c(min_concordance,
  #     0.997,
  #     # mean(c(min_concordance, max_concordance)),
  #     max_concordance),
  #   c( "#3182bd", "white", "orange")
  # ),
  top_annotation       = col_ha_target,
  # left_annotation      = row_ha,
  cluster_rows         = TRUE,
  cluster_columns      = TRUE,
  clustering_distance_rows = function(m) as.dist(1 - m),
  clustering_distance_columns = function(m) as.dist(1 - m),
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  show_row_names       = FALSE,
  show_column_names    = TRUE,
  row_dend_gp      = gpar(lwd = 0.3),
  column_dend_gp   = gpar(lwd = 0.3),
  
  
  # show_column_names    = TRUE,
  column_names_rot     = 90,
  column_names_gp      = gpar(fontsize = 1.5),
  column_names_side    = "top",
  column_names_max_height = unit(0.15, "cm"),
  
  
  # column_title         = paste0(target_isotype," isotypes"),
  column_title_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_title_rot     = 0
  
)


dev.off()
pdf(paste0("../../plots/FigureS13_Ct_",target_isotype,"_concordance_heatmap_3.pdf"), width = 7, height = 6)


draw(p_heatmap_3, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()






