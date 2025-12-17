rm(list = ls())

library(DECIPHER)  ## BiocManager::install("DECIPHER")
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize) 


source("../utilities.R")


gtcheck_raw <- read.delim("../../processed_data/heatmap_hard_filtered/gtcheck.tsv") 

#### output ###
gtcheck_output<-gtcheck_raw %>% 
  select(discordance,sites,i,j) %>% 
  rename(`identical alleles` = discordance) %>% 
  rename(`total genome-wide SNVs` = sites) %>% 
  rename(`strain1` = i) %>% 
  rename(`strain2` = j)
write.csv(gtcheck_output,"../../tables/TableS2.csv",row.names = FALSE)
#### output end ###



gtcheck<-gtcheck_raw %>% 
  dplyr::mutate(concordance = 1-(discordance/sites)) %>% 
  select(i,j,concordance)

gtcheck$concordance %>% max()
# [1] 0.9999146
gtcheck$concordance %>% min()
# [1] 0.6865332


#### threshold
p_histogram<-ggplot()+ geom_histogram(data=gtcheck,aes(x=concordance),
                         bins = 1000)+
  scale_x_continuous(breaks = seq(0.65, 1, by = 0.01))+
  theme_bw()+
  xlab("Genetic similarity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_histogram

# ggsave("Ct_gt_check_histogram.pdf",p_histogram,
#        width = 7,height = 4, units = "in")



# convert into matrix
gt_matrix <- gtcheck %>%
  tidyr::pivot_wider(names_from = j, values_from = concordance) %>%
  tibble::column_to_rownames("i") %>%
  as.matrix()




# add new row
new_row <- matrix(NA, nrow = 1, ncol = ncol(gt_matrix))
rownames(new_row) <- colnames(gt_matrix)[1]
gt_matrix <- rbind(new_row, gt_matrix)

# add new col
new_col <- matrix(NA, ncol = 1, nrow = nrow(gt_matrix))
colnames(new_col) <- rownames(gt_matrix)[nrow(gt_matrix)]
gt_matrix <- cbind(gt_matrix,new_col)





# add upper.tri data
gt_matrix[upper.tri(gt_matrix,diag = FALSE)] <- t(gt_matrix)[upper.tri(gt_matrix,diag = FALSE)]




# sort
sorted_rows <- sort(rownames(gt_matrix))
sorted_cols <- sort(colnames(gt_matrix))

gt_matrix_plot <- gt_matrix[sorted_rows, sorted_cols]
diag(gt_matrix_plot) <- 1






########### Annotations #######
### Geo
geo_info_raw<-read.csv("../../processed_data/Geo_info/Ct_indep_isotype_info_geo.csv")
geo_info<-geo_info_raw %>%
  tibble::column_to_rownames(var = "isotype") %>% 
  dplyr::select(geo) %>% 
  dplyr::rename(Geo=geo)
geo_info$Geo<-as.factor(geo_info$Geo)
# geo.colours<-as.list(geo.colours)
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


# calculate min and max
min_concordance <- min(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)
max_concordance <- max(gt_matrix_plot[gt_matrix_plot != 0], na.rm = TRUE)










################################################
############ three breaks ######################
################################################



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

phylo_col_fun <- circlize::colorRamp2(
  breaks,
  c("#4575B4", "#87C6C2", "#FFFFE0","#F4D166","#D73027")
)



p_heatmap_3 <- Heatmap(
  matrix               = gt_matrix_plot,
  name                 = "Concordance",
  col                  = phylo_col_fun,
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


#### heatmap order
ht_list <- draw(p_heatmap_3)
col_idx <- column_order(ht_list)
ordered_colnames <- colnames(gt_matrix_plot)[col_idx]


extract_sample_list<-function(start_isotype,end_isotype){
start <- match(start_isotype, ordered_colnames)
end   <- match(end_isotype, ordered_colnames)
if(is.na(start) || is.na(end)) stop("no such sample")

if(start <= end) subset_ordered <- ordered_colnames[start:end] else subset_ordered <- ordered_colnames[end:start]
subset_ordered
return(subset_ordered)
}


LAC<-extract_sample_list(start_isotype = "NIC1257",
                           end_isotype = "NIC1939")
Tw1<-extract_sample_list(start_isotype = "NIC1658",
                         end_isotype = "NIC1603")
Tw2<-extract_sample_list(start_isotype = "NIC1560",
                         end_isotype = "TWN2103")
Ma1<-c("HPT56")
Tw3<-extract_sample_list(start_isotype = "NIC1651",
                         end_isotype = "NIC1656")
Ma2<-c("HPT25","HPT39")
Tw4<-extract_sample_list(start_isotype = "NIC895",
                         end_isotype = "NIC898")

Mic1<-extract_sample_list(start_isotype = "QG4867",
                          end_isotype = "QG4739")
Mic2<-extract_sample_list(start_isotype = "QG5551",
                          end_isotype = "QG5319")
HC<-extract_sample_list(start_isotype = "NIC1950",
                         end_isotype = "ECA794")
Hw1<-extract_sample_list(start_isotype = "ECA803",
                        end_isotype = "ECA1518")
Tw5<-c("BRC20391","NIC1563")
Au<-extract_sample_list(start_isotype = "QG2899",
                        end_isotype = "QG2897")
Tw6<-c("NIC1592","NIC1594")
Hw2<-extract_sample_list(start_isotype = "ECA1449",
                         end_isotype = "ECA1453")
Af<-extract_sample_list(start_isotype = "NIC1534",
                         end_isotype = "JU3170")
Ma3<-c("HPT33")
Tw7<-c("NIC1679")
Hw3<-extract_sample_list(start_isotype = "ECA1307",
                         end_isotype = "ECA1299")


lineage_samples_list <- list(
  LAC = LAC,
  Tw1 = Tw1,
  Tw2 = Tw2,
  Ma1 = Ma1,
  Tw3 = Tw3,
  Ma2 = Ma2,
  Tw4 = Tw4,
  # Mic = Mic,
  Mic1 = Mic1,
  Mic2 = Mic2,
  HC = HC,
  Hw1 = Hw1,
  Tw5 = Tw5,
  Au = Au,
  Tw6 = Tw6,
  Hw2 = Hw2,
  Af = Af,
  Ma3 = Ma3,
  Tw7 = Tw7,
  Hw3 = Hw3
)

df_lineage_all <- dplyr::bind_rows(
  lapply(names(lineage_samples_list), function(lineage_name) {
    tibble(sample = lineage_samples_list[[lineage_name]], lineage = lineage_name)
  })
)


df_lineage_all$lineage %>% unique() %>% length()

write.csv(df_lineage_all,
          "../../processed_data/geo_info/Ct_lineage_all.csv",
          row.names = FALSE)












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


dev.off()
pdf("../../plots/Figure_2_raw_Ct_hard_filtered_heatmap_with_lineage.pdf", width = 10, height = 9)


draw(p_heatmap_4, merge_legend = TRUE, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()









library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(circlize)
library(tibble)


upper_mat <- gt_matrix_plot
upper_mat[lower.tri(upper_mat, diag = TRUE)] <- NA 


upper_long <- as.data.frame(upper_mat) %>%
  mutate(Row = rownames(upper_mat)) %>%
  pivot_longer(
    cols = -Row,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  drop_na(Value)


upper_long_lineage <- upper_long %>%
  left_join(df_lineage_all, by = c("Row" = "sample")) %>%
  rename(lineage_row = lineage) %>%
  left_join(df_lineage_all, by = c("Column" = "sample")) %>%
  rename(lineage_col = lineage)


lineage_mean_sym <- upper_long_lineage %>%
  mutate(
    lineage_min = pmin(lineage_row, lineage_col),
    lineage_max = pmax(lineage_row, lineage_col)
  ) %>%
  group_by(lineage_min, lineage_max) %>%
  summarise(mean_value = mean(Value), .groups = "drop")


mat_lineage <- lineage_mean_sym %>%
  pivot_wider(
    names_from  = lineage_max,
    values_from = mean_value,
    values_fill = NA
  ) %>%
  column_to_rownames("lineage_min") %>%
  as.matrix()


mat_lineage <- rbind(mat_lineage, Tw7 = NA)


mat_full <- mat_lineage
mat_full[lower.tri(mat_full)] <- t(mat_full)[lower.tri(mat_full)]


df_long <- mat_full %>%
  as.data.frame() %>%
  mutate(Row = rownames(mat_full)) %>%
  pivot_longer(
    cols = -Row,
    names_to = "Column",
    values_to = "Value"
  )


lower_flag <- lower.tri(mat_full, diag = TRUE)

df_lower <- df_long %>%
  filter(
    lower_flag[
      cbind(
        match(Row, rownames(mat_full)),
        match(Column, colnames(mat_full))
      )
    ]
  ) %>%
  mutate(
    Row    = factor(Row,    levels = rev(rownames(mat_full))),
    Column = factor(Column, levels = colnames(mat_full))
  )


phylo_vals <- df_lower$Value
phylo_vals <- phylo_vals[!is.na(phylo_vals)]

breaks <- c(
  min(phylo_vals, na.rm = TRUE),
  quantile(phylo_vals, 0.25, na.rm = TRUE),
  median(phylo_vals, na.rm = TRUE),
  quantile(phylo_vals, 0.75, na.rm = TRUE),
  max(phylo_vals, na.rm = TRUE)
)

cols <- c("#4575B4", "#87C6C2", "#FFFFE0", "#F4D166", "#D73027")
vals <- scales::rescale(breaks)


p <- ggplot(df_lower, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = cols,
    values  = vals,
    limits  = range(breaks),
    name    = "Genetic\nsimilarity"
  ) +
  geom_text(
    aes(label = round(Value, 2)), 
    na.rm = TRUE,
    size = 3
  ) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    legend.position = c(0.75, 0.75)
  )

p

ggsave("../../plots/Figure_S4_mean_concordance.pdf",
       plot = p, width = 7.5, height = 7.5, 
       units = "in", device = 'pdf')


#






df_lineage_all_tmp<- df_lineage_all %>% 
  rename(isotype=sample)
  
geo_and_lineage<-geo_info_raw %>% 
  left_join(df_lineage_all_tmp,by = "isotype")

write.csv(geo_and_lineage,
          "../../processed_data/geo_info/geo_and_lineage.csv",
          row.names = FALSE,
          quote = FALSE)


LAC_geo_and_lineage<-geo_and_lineage %>% 
  filter(lineage == "LAC")

LAC_geo_and_lineage_car<-LAC_geo_and_lineage %>% 
  filter(LAC_geo_and_lineage$geo == "Caribbean")
nrow(LAC_geo_and_lineage_car)
# [1] 47





HC_geo_and_lineage<-geo_and_lineage %>% 
  filter(lineage == "HC")
nrow(HC_geo_and_lineage)
# [1] 35

Hw_geo_and_lineage_car<-HC_geo_and_lineage %>% 
  filter(geo == "Hawaii")
nrow(Hw_geo_and_lineage_car)
# [1] 18

Ca_geo_and_lineage_car<-HC_geo_and_lineage %>% 
  filter(geo == "Caribbean")
nrow(Ca_geo_and_lineage_car)
# [1] 17




Hw1_geo_and_lineage<-geo_and_lineage %>% 
  filter(lineage == "Hw1")
nrow(Hw1_geo_and_lineage)
# [1] 23



