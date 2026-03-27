rm(list = ls())

library(dplyr)
library(readr)
library(geosphere)
library(vegan)
library(ggplot2)
library(reshape2)

gtcheck_isotype_raw <- read.table("../../data/gtcheck.txt",
                                  header = TRUE)

gtcheck_isotype <- gtcheck_isotype_raw %>%
  mutate(concordance = 1 - (discordance / sites)) %>%
  select(i, j, concordance)

Ct_raw_csv <- readr::read_tsv("../../tables/Table_S1.tsv")

Ct_isotyeps_csv <- Ct_raw_csv %>%
  dplyr::filter(strain == isotype) %>%
  dplyr::select(isotype, latitude, longitude) %>%
  na.omit() %>%
  dplyr::mutate(
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  )

coords <- Ct_isotyeps_csv[, c("longitude", "latitude")]
geo_mat <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000  # km
rownames(geo_mat) <- Ct_isotyeps_csv$isotype
colnames(geo_mat) <- Ct_isotyeps_csv$isotype

gt2 <- gtcheck_isotype %>%
  dplyr::mutate(
    iso1 = pmin(as.character(i), as.character(j)),
    iso2 = pmax(as.character(i), as.character(j))
  ) %>%
  dplyr::select(iso1, iso2, concordance) %>%
  dplyr::distinct()

iso_keep <- intersect(rownames(geo_mat), unique(c(gt2$iso1, gt2$iso2)))
geo_mat2 <- geo_mat[iso_keep, iso_keep]
gen_mat2 <- matrix(NA_real_, nrow = length(iso_keep), ncol = length(iso_keep),
                   dimnames = list(iso_keep, iso_keep))

idx1 <- match(gt2$iso1, iso_keep)
idx2 <- match(gt2$iso2, iso_keep)
ok <- !is.na(idx1) & !is.na(idx2)

gen_mat2[cbind(idx1[ok], idx2[ok])] <- 1 - gt2$concordance[ok]
gen_mat2[cbind(idx2[ok], idx1[ok])] <- 1 - gt2$concordance[ok]
diag(gen_mat2) <- 0

keep_complete <- iso_keep[rowSums(is.na(gen_mat2)) == 0]
geo_mat3 <- geo_mat2[keep_complete, keep_complete]
gen_mat3 <- gen_mat2[keep_complete, keep_complete]

mant <- vegan::mantel(as.dist(gen_mat3), as.dist(geo_mat3),
                      method = "spearman", permutations = 999)
mant_r <- unname(mant$statistic)
mant_p <- unname(mant$signif)

mant_label <- paste0("Mantel Spearman r = ", round(mant_r, 3),
                     "\nP = ", signif(mant_p, 3)
                     )

idx <- which(lower.tri(geo_mat3), arr.ind = TRUE)
plot_df <- data.frame(
  isoA = rownames(geo_mat3)[idx[, 1]],
  isoB = colnames(geo_mat3)[idx[, 2]],
  geo_km = geo_mat3[idx],
  gen_dist = gen_mat3[idx]
) %>%
  filter(!is.na(geo_km), !is.na(gen_dist))

p <- ggplot(plot_df, aes(x = geo_km, y = gen_dist)) +
  geom_point(alpha = 0.35, size = 0.5, shape = 16) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text",
           x = Inf, y = Inf,
           label = mant_label,
           hjust = 1.05, vjust = 1.1,
           size = 4.5) +
  labs(
    x = "Geographic distance (km)",
    y = "Genetic distance (1 - genetic similarity)",
    title = NULL
  ) +
  theme_classic()

ggsave("../../figures/FigureS12_IBD_isotype_pairs_Mantel.png",
       plot = p,
       width = 7,
       height = 7,
       dpi = 300)







rm(list = ls())

library(dplyr)
library(readr)
library(geosphere)
library(vegan)
library(ggplot2)

gtcheck_isotype_raw <- read.table("../../data/gtcheck.txt",
                                  header = TRUE)

gtcheck_isotype <- gtcheck_isotype_raw %>%
  dplyr::mutate(concordance = 1 - (discordance / sites)) %>%
  dplyr::select(i, j, concordance)

Ct_raw_csv <- readr::read_tsv("../../tables/Table_S1.tsv")

Ct_isotypes <- Ct_raw_csv %>%
  dplyr::filter(strain == isotype) %>%
  dplyr::select(isotype, latitude, longitude, Relatedness_group) %>%
  na.omit() %>%
  dplyr::mutate(
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude)
  )

gt2 <- gtcheck_isotype %>%
  dplyr::mutate(
    iso1 = pmin(as.character(i), as.character(j)),
    iso2 = pmax(as.character(i), as.character(j))
  ) %>%
  dplyr::select(iso1, iso2, concordance) %>%
  dplyr::distinct()

run_ibd_one_rg <- function(rg_name, df_iso, gt_pairs, permutations = 999) {
  
  df_rg <- df_iso %>% filter(Relatedness_group == rg_name)
  iso_rg <- df_rg$isotype
  
  if (length(iso_rg) < 3) {
    return(list(rg = rg_name, kept_complete = 0,
                mantel_r = NA_real_, mantel_p = NA_real_,
                plot_df = NULL))
  }
  
  coords <- df_rg[, c("longitude", "latitude")]
  geo_mat <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
  rownames(geo_mat) <- df_rg$isotype
  colnames(geo_mat) <- df_rg$isotype
  
  gt_rg <- gt_pairs %>% filter(iso1 %in% iso_rg, iso2 %in% iso_rg)
  
  iso_keep <- intersect(rownames(geo_mat), unique(c(gt_rg$iso1, gt_rg$iso2)))
  if (length(iso_keep) < 3) {
    return(list(rg = rg_name, kept_complete = 0,
                mantel_r = NA_real_, mantel_p = NA_real_,
                plot_df = NULL))
  }
  
  geo_mat2 <- geo_mat[iso_keep, iso_keep]
  
  gen_mat2 <- matrix(NA_real_, nrow = length(iso_keep), ncol = length(iso_keep),
                     dimnames = list(iso_keep, iso_keep))
  
  idx1 <- match(gt_rg$iso1, iso_keep)
  idx2 <- match(gt_rg$iso2, iso_keep)
  ok <- !is.na(idx1) & !is.na(idx2)
  
  gen_mat2[cbind(idx1[ok], idx2[ok])] <- 1 - gt_rg$concordance[ok]
  gen_mat2[cbind(idx2[ok], idx1[ok])] <- 1 - gt_rg$concordance[ok]
  diag(gen_mat2) <- 0
  
  keep_complete <- iso_keep[rowSums(is.na(gen_mat2)) == 0]
  if (length(keep_complete) < 3) {
    return(list(rg = rg_name, kept_complete = length(keep_complete),
                mantel_r = NA_real_, mantel_p = NA_real_,
                plot_df = NULL))
  }
  
  geo_mat3 <- geo_mat2[keep_complete, keep_complete]
  gen_mat3 <- gen_mat2[keep_complete, keep_complete]

  mant <- vegan::mantel(as.dist(gen_mat3), as.dist(geo_mat3),
                        method = "spearman", permutations = permutations)
  mantel_r <- unname(mant$statistic)
  mantel_p <- unname(mant$signif)
  
  idx <- which(lower.tri(geo_mat3), arr.ind = TRUE)
  plot_df <- data.frame(
    RG = rg_name,
    geo_km = geo_mat3[idx],
    gen_dist = gen_mat3[idx]
  ) %>%
    dplyr::filter(!is.na(geo_km), !is.na(gen_dist)) %>%
    dplyr::mutate(
      mantel_r = mantel_r,
      mantel_p = mantel_p,
      N = length(keep_complete)
    )
  
  list(rg = rg_name,
       kept_complete = length(keep_complete),
       mantel_r = mantel_r,
       mantel_p = mantel_p,
       plot_df = plot_df)
}

rg_levels <- sort(unique(Ct_isotypes$Relatedness_group))
res_list <- lapply(rg_levels, run_ibd_one_rg,
                   df_iso = Ct_isotypes,
                   gt_pairs = gt2,
                   permutations = 999)

res_df <- bind_rows(lapply(res_list, function(x) {
  data.frame(
    rg = x$rg,
    N = x$kept_complete,
    mantel_r = x$mantel_r,
    mantel_p = x$mantel_p
  )
})) %>%
  arrange(desc(N))

plot_all_df <- bind_rows(lapply(res_list, function(x) x$plot_df)) %>%
  dplyr::filter(!is.na(RG))

plot_all_df2 <- plot_all_df %>%
  dplyr::filter(N > 10)

rg_keep <- plot_all_df2 %>%
  dplyr::group_by(RG) %>%
  dplyr::summarise(N = first(N), .groups="drop") %>%
  dplyr::arrange(desc(N)) %>%
  dplyr::pull(RG)

plot_all_df2$RG <- factor(plot_all_df2$RG, levels = rg_keep)

labels_df <- plot_all_df2 %>%
  dplyr::group_by(RG) %>%
  dplyr::summarise(
    mantel_r = first(mantel_r),
    mantel_p = first(mantel_p),
    N = first(N),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0(
      "Mantel Spearman r = ", round(mantel_r, 3),
      "\n P = ", signif(mantel_p, 2),
      "\n N = ", N
    )
  )

p_big <- ggplot(plot_all_df2, aes(x = geo_km, y = gen_dist)) +
  geom_point(alpha = 0.25, size = 0.4, shape = 16) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_text(
    data = labels_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 3.4,
    color = "red",
    alpha = 0.7
  ) +
  facet_wrap(~RG, scales = "free", ncol = 3) +
  labs(
    x = "Geographic distance (km)",
    y = "Genetic distance (1 - concordance)",
    title = NULL
  ) +
  theme_classic() +
  theme(strip.text = element_text(size = 11))

ggsave("../../figures/FigureS13_IBD_RGs.png",
       plot = p_big, 
       width = 7, height = 7, 
       dpi = 300)



