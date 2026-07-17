rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(purrr)

## ============================================================
## Helper functions
## ============================================================

make_ct_region_rects <- function(domains_file) {
  
  domains_raw <- readr::read_tsv(domains_file, show_col_types = FALSE)
  
  domains_wide <- domains_raw %>%
    dplyr::filter(Location %in% c("left_arm", "right_arm")) %>%
    dplyr::mutate(arm = ifelse(Location == "left_arm", "left", "right")) %>%
    dplyr::select(CHROM, arm, start, end) %>%
    tidyr::pivot_wider(
      names_from = arm,
      values_from = c(start, end),
      names_glue = "{arm}_{.value}"
    ) %>%
    dplyr::ungroup()
  
  region_rects <- domains_wide %>%
    tidyr::pivot_longer(
      cols = c(left_start, right_end),
      names_to = "region_side",
      values_to = "x"
    ) %>%
    dplyr::mutate(
      region = "Tip",
      side = ifelse(region_side == "left_start", "Left", "Right"),
      xmin = ifelse(region_side == "left_start", 0, x),
      xmax = ifelse(region_side == "left_start", x, Inf)
    ) %>%
    dplyr::select(CHROM, region, side, xmin, xmax) %>%
    dplyr::bind_rows(
      domains_wide %>%
        dplyr::transmute(
          CHROM,
          region = "Arm",
          side = "Left",
          xmin = left_start,
          xmax = left_end
        ),
      domains_wide %>%
        dplyr::transmute(
          CHROM,
          region = "Arm",
          side = "Right",
          xmin = right_start,
          xmax = right_end
        ),
      domains_wide %>%
        dplyr::transmute(
          CHROM,
          region = "Center",
          side = "Center",
          xmin = left_end,
          xmax = right_start
        )
    ) %>%
    dplyr::mutate(
      ymin = -Inf,
      ymax = Inf,
      xmin = xmin / 1e6,
      xmax = xmax / 1e6
    )
  
  return(region_rects)
}


make_cb_region_rects <- function(domains_file) {
  
  genome_domain_raw <- readr::read_csv(domains_file, show_col_types = FALSE)
  
  region_rects <- genome_domain_raw %>%
    dplyr::mutate(
      CHROM = chrom,
      region = dplyr::case_when(
        grepl("tip$", sub_region) ~ "Tip",
        grepl("arm$", sub_region) ~ "Arm",
        sub_region == "center" ~ "Center",
        TRUE ~ NA_character_
      ),
      side = dplyr::case_when(
        grepl("^left", sub_region) ~ "Left",
        grepl("^right", sub_region) ~ "Right",
        sub_region == "center" ~ "Center",
        TRUE ~ NA_character_
      ),
      xmin = start / 1e6,
      xmax = stop / 1e6,
      ymin = -Inf,
      ymax = Inf
    ) %>%
    dplyr::filter(!is.na(region)) %>%
    dplyr::select(CHROM, region, side, xmin, xmax, ymin, ymax)
  
  return(region_rects)
}


make_hdr_dt <- function(hdr_file, source_filter = NULL, collapse_by_chrom = FALSE) {
  
  hdrs <- readr::read_tsv(hdr_file, show_col_types = FALSE)
  
  if (!is.null(source_filter)) {
    hdrs <- hdrs %>%
      dplyr::filter(source == source_filter)
  }
  
  if (collapse_by_chrom) {
    
    hdr_dt <- hdrs %>%
      dplyr::group_by(CHROM) %>%
      dplyr::arrange(minStart, maxEnd, .by_group = TRUE) %>%
      dplyr::mutate(
        new_group = cumsum(minStart > dplyr::lag(cummax(maxEnd), default = first(minStart) - 1))
      ) %>%
      dplyr::group_by(CHROM, new_group) %>%
      dplyr::summarise(
        minStart = min(minStart),
        maxEnd = max(maxEnd),
        .groups = "drop"
      )
  }
  
  hdr_dt <- hdrs %>%
    dplyr::mutate(
      tstart = minStart / 1e6,
      tend   = maxEnd / 1e6
    ) %>%
    dplyr::select(CHROM, tstart, tend) %>%
    dplyr::distinct() %>%
    as.data.table()
  
  data.table::setkey(hdr_dt, CHROM, tstart, tend)
  
  return(hdr_dt)
}


annotate_overlaps <- function(combined_df, region_rects) {
  
  combined_dt <- as.data.table(
    combined_df %>%
      dplyr::mutate(
        window_start = window_start / 1e6,
        window_stop  = window_stop / 1e6
      )
  )
  
  region_dt <- as.data.table(region_rects)
  
  data.table::setnames(combined_dt, c("window_start", "window_stop"), c("start", "end"))
  data.table::setnames(region_dt, c("xmin", "xmax"), c("start", "end"))
  
  data.table::setkey(region_dt, CHROM, start, end)
  data.table::setkey(combined_dt, chrom, start, end)
  
  annotated_dt <- data.table::foverlaps(
    region_dt,
    combined_dt,
    by.x = c("CHROM", "start", "end"),
    by.y = c("chrom", "start", "end"),
    nomatch = 0,
    type = "any"
  )
  
  return(annotated_dt)
}


process_dxy_annotation <- function(dxy_df, region_rects, hdr_dt) {
  
  annotated_dxy_lin <- annotate_overlaps(dxy_df, region_rects) %>%
    dplyr::mutate(
      region_side = ifelse(side != "Center", paste0(side, "_", region), region)
    ) %>%
    dplyr::filter(region != "Tip") %>%
    dplyr::mutate(
      region_side = factor(
        region_side,
        levels = c("Left_Arm", "Center", "Right_Arm"),
        labels = c("Left Arm", "Center", "Right Arm")
      )
    )
  
  data.table::setDT(annotated_dxy_lin)
  
  x <- data.table::copy(annotated_dxy_lin)[, `:=`(qstart = start, qend = end)]
  
  data.table::setkey(x, CHROM, qstart, qend)
  data.table::setkey(hdr_dt, CHROM, tstart, tend)
  
  ov_dxy <- data.table::foverlaps(
    x = x,
    y = hdr_dt,
    by.x = c("CHROM", "qstart", "qend"),
    by.y = c("CHROM", "tstart", "tend"),
    type = "any",
    nomatch = 0L
  )
  
  ov_dxy[, `:=`(
    ov_len = pmax(0, pmin(qend, tend) - pmax(qstart, tstart)),
    int_len = pmax(0, qend - qstart)
  )]
  
  ov_dxy[, prop := data.table::fifelse(int_len > 0, ov_len / int_len, 0)]
  
  ov_dxy_best <- as.data.frame(ov_dxy) %>%
    dplyr::select(CHROM, start, end, prop) %>%
    dplyr::group_by(CHROM, start, end) %>%
    dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  anno_dxy <- as.data.frame(annotated_dxy_lin) %>%
    dplyr::left_join(ov_dxy_best, by = c("CHROM", "start", "end")) %>%
    dplyr::mutate(
      hdstatus = ifelse(is.na(prop) | prop < 0.5, "non-HDR", "HDR"),
      hdstatus = factor(hdstatus, levels = c("non-HDR", "HDR")),
      stat_type = "Dxy"
    )
  
  return(anno_dxy)
}


make_dxy_plot <- function(anno_dxy_df) {
  
  p <- ggplot(anno_dxy_df) +
    geom_jitter(
      aes(
        y = hdstatus,
        x = stat,
        colour = hdstatus
      ),
      size = 0.5,
      alpha = 0.2
    ) +
    geom_boxplot(
      aes(
        y = hdstatus,
        x = stat
      ),
      colour = "grey40",
      outliers = FALSE,
      width = 0.2,
      alpha = 0.6
    ) +
    facet_grid(comp_lin ~ region) +
    theme_bw(base_family = "Helvetica") +
    ylab("") +
    xlab(expression(D[xy])) +
    scale_colour_manual(values = c("non-HDR" = "black", "HDR" = "#0719BC")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 5, 5, 5),
      text = element_text(family = "Helvetica"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}


summarise_nonhdr_dxy <- function(anno_dxy_df, comp_levels, focal_pop, species_name) {
  
  out <- anno_dxy_df %>%
    dplyr::filter(hdstatus == "non-HDR") %>%
    dplyr::group_by(comp_lin) %>%
    dplyr::summarise(
      n_windows   = dplyr::n(),
      mean_dxy    = mean(stat, na.rm = TRUE),
      median_dxy  = median(stat, na.rm = TRUE),
      sd_dxy      = sd(stat, na.rm = TRUE),
      min_dxy     = min(stat, na.rm = TRUE),
      max_dxy     = max(stat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      comp_lin = factor(comp_lin, levels = comp_levels)
    ) %>%
    dplyr::arrange(comp_lin) %>%
    dplyr::rename(
      mean_dxy_nonHDR = mean_dxy,
      group = comp_lin
    ) %>%
    dplyr::select(group, mean_dxy_nonHDR) %>%
    dplyr::mutate(
      group = paste0(group, " and ", focal_pop),
      species = species_name
    ) %>%
    dplyr::select(species, group, mean_dxy_nonHDR)
  
  return(out)
}


check_dxy_distribution <- function(dxy_df) {
  
  out <- dxy_df %>%
    dplyr::group_by(comp_lin) %>%
    dplyr::summarise(
      n_windows = dplyr::n(),
      n_NA = sum(is.na(stat)),
      n_zero = sum(stat == 0, na.rm = TRUE),
      prop_zero = n_zero / sum(!is.na(stat)),
      min_nonzero = min(stat[stat > 0], na.rm = TRUE),
      median_dxy = median(stat, na.rm = TRUE),
      mean_dxy = mean(stat, na.rm = TRUE),
      max_dxy = max(stat, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(out)
}


## ============================================================
## Ct domains and HDRs
## ============================================================

ct_region_rects <- make_ct_region_rects(
  "../../data/ct_ch_domains.tsv"
)

ct_hdr_dt <- make_hdr_dt(
  "../../tables/TableS7_HDR_CT_allStrain_5kbclust_20251201.tsv"
)


## ============================================================
## Ct Dxy
## ============================================================

ct_dxy_raw <- readr::read_tsv(
  "../../processed_data/Dxy_LAC/results/Ct_Dxy__dxy.txt",
  show_col_types = FALSE
)

focal_pop_ct <- "LAC"

ct_comp_lin_levels <- c(
  "Hw3", "Af", "Hw1", "HC",
  "Mic2", "Mic1", "Tw3", "Tw1"
)

ct_dxy_lac <- ct_dxy_raw %>%
  dplyr::filter(pop1 == focal_pop_ct | pop2 == focal_pop_ct) %>%
  dplyr::mutate(
    comp_lin = dplyr::case_when(
      pop1 == focal_pop_ct ~ pop2,
      pop2 == focal_pop_ct ~ pop1,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(comp_lin)) %>%
  dplyr::rename(
    chrom        = chromosome,
    window_start = window_pos_1,
    window_stop  = window_pos_2,
    stat         = avg_dxy
  ) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::filter(comp_lin %in% ct_comp_lin_levels) %>%
  dplyr::mutate(
    comp_lin = factor(comp_lin, levels = ct_comp_lin_levels)
  )

anno_list_ct <- split(ct_dxy_lac, ct_dxy_lac$comp_lin) %>%
  purrr::map(~ process_dxy_annotation(.x, ct_region_rects, ct_hdr_dt))

anno_dxy_ct_all <- dplyr::bind_rows(anno_list_ct) %>%
  dplyr::mutate(
    comp_lin = factor(comp_lin, levels = ct_comp_lin_levels)
  )

dxybox_ct <- make_dxy_plot(anno_dxy_ct_all)

ggsave(
  dxybox_ct,
  filename = "../../figures/FigureS31_Dxy.png",
  width = 7,
  height = 7,
  units = "in",
  dpi = 600
)

Ct_nonhdr_dxy_overall <- summarise_nonhdr_dxy(
  anno_dxy_df = anno_dxy_ct_all,
  comp_levels = ct_comp_lin_levels,
  focal_pop = focal_pop_ct,
  species_name = "C. tropicalis"
)

ct_dxy_distribution_check <- check_dxy_distribution(ct_dxy_lac)

print(ct_dxy_distribution_check)


## ============================================================
## Cb domains and HDRs
## ============================================================

cb_region_rects <- make_cb_region_rects(
  "../../data/meta_data_Ce_Cb/Cb/Cb_bounds_df.csv"
)

cb_hdr_dt <- make_hdr_dt(
  "../../data/meta_data_Ce_Cb/Cb/HDR_CB_allStrain_5kbclust_20250930.tsv",
  source_filter = "QX1410",
  collapse_by_chrom = TRUE
)


## ============================================================
## Cb Dxy
## ============================================================

cb_dxy_raw <- readr::read_tsv(
  "../../data/meta_data_Ce_Cb/Cb/Cb_Dxy/ALLCB_dxy.txt",
  show_col_types = FALSE
)

focal_pop_cb <- "Tropical"

cb_comp_lin_levels <- c(
  "AD",
  "KD",
  "TD1",
  "Temperate",
  "TH"
)

cb_dxy_focal <- cb_dxy_raw %>%
  dplyr::filter(pop1 == focal_pop_cb | pop2 == focal_pop_cb) %>%
  dplyr::mutate(
    comp_lin = dplyr::case_when(
      pop1 == focal_pop_cb ~ pop2,
      pop2 == focal_pop_cb ~ pop1,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(comp_lin)) %>%
  dplyr::rename(
    chrom        = chromosome,
    window_start = window_pos_1,
    window_stop  = window_pos_2,
    stat         = avg_dxy
  ) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::filter(comp_lin %in% cb_comp_lin_levels) %>%
  dplyr::mutate(
    comp_lin = factor(comp_lin, levels = cb_comp_lin_levels)
  )

anno_list_cb <- split(cb_dxy_focal, cb_dxy_focal$comp_lin) %>%
  purrr::map(~ process_dxy_annotation(.x, cb_region_rects, cb_hdr_dt))

anno_dxy_cb_all <- dplyr::bind_rows(anno_list_cb) %>%
  dplyr::mutate(
    comp_lin = factor(comp_lin, levels = cb_comp_lin_levels)
  )

dxybox_cb <- make_dxy_plot(anno_dxy_cb_all)

Cb_nonhdr_dxy_overall <- summarise_nonhdr_dxy(
  anno_dxy_df = anno_dxy_cb_all,
  comp_levels = cb_comp_lin_levels,
  focal_pop = focal_pop_cb,
  species_name = "C. briggsae"
)

cb_dxy_distribution_check <- check_dxy_distribution(cb_dxy_focal)

print(cb_dxy_distribution_check)

TableS10_Ct_and_Cb_Dxy_nonHDRs <- dplyr::bind_rows(
  Ct_nonhdr_dxy_overall,
  Cb_nonhdr_dxy_overall
) %>%
  dplyr::mutate(
    species = factor(
      species,
      levels = c("C. tropicalis", "C. briggsae"),
    ),
    mean_dxy_nonHDR = round(mean_dxy_nonHDR, 4)
  ) %>%
  dplyr::arrange(species) %>%
  dplyr::select(species, group, mean_dxy_nonHDR)

readr::write_csv(
  TableS10_Ct_and_Cb_Dxy_nonHDRs,
  "../../tables/TableS10_Ct_and_Cb_Dxy_nonHDRs.csv"
)

print(TableS10_Ct_and_Cb_Dxy_nonHDRs)
