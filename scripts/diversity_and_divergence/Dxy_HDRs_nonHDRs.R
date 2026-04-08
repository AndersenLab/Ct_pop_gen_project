rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(purrr)

domains_raw <- readr::read_tsv("../../data/ct_ch_domains.tsv") 

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
    xmin = ifelse(region_side == "left_start", 0, x),
    xmax = ifelse(region_side == "left_start", x, Inf)
  ) %>%
  dplyr::select(CHROM, region, xmin, xmax) %>%
  dplyr::bind_rows(
    domains_wide %>%
      dplyr::mutate(region = "Arm") %>%
      dplyr::transmute(CHROM, region, xmin = left_start, xmax = left_end),
    domains_wide %>%
      dplyr::mutate(region = "Arm") %>%
      dplyr::transmute(CHROM, region, xmin = right_start, xmax = right_end),
    domains_wide %>%
      dplyr::mutate(region = "Center") %>%
      dplyr::transmute(CHROM, region, xmin = left_end, xmax = right_start)
  ) %>%
  dplyr::mutate(
    ymin = -Inf,
    ymax = Inf,
    xmin = xmin / 1e6,
    xmax = xmax / 1e6
  )

region_rects2 <- region_rects %>%
  dplyr::arrange(region, CHROM, xmin) %>%
  dplyr::mutate(
    side = dplyr::case_when(
      region == "Center" ~ "Center",
      row_number() %% 2 == 0 ~ "Right",
      TRUE ~ "Left"
    )
  )

region_colors <- c("Tip" = "#5E3C99", "Center" = "#FDB863", "Arm" = "#4393C3")
hdrs_ct <- readr::read_tsv("../../tables/TableS6_HDR_CT_allStrain_5kbclust_20251201.tsv")

y <- hdrs_ct %>%
  dplyr::mutate(
    tstart = minStart / 1e6,
    tend   = maxEnd   / 1e6
  ) %>%
  dplyr::select(CHROM, tstart, tend) %>%
  dplyr::distinct() %>%
  as.data.table()

setkey(y, CHROM, tstart, tend)

annotate_overlaps <- function(combined_df, region_rects) {
  combined_dt <- as.data.table(combined_df %>%
                                 mutate(
                                   window_start = window_start / 1e6,
                                   window_stop  = window_stop  / 1e6
                                 ))
  region_dt <- as.data.table(region_rects)
  setnames(combined_dt, c("window_start", "window_stop"), c("start", "end"))
  setnames(region_dt, c("xmin", "xmax"), c("start", "end"))
  setkey(region_dt, CHROM, start, end)
  setkey(combined_dt, chrom, start, end)
  annotated_dt <- foverlaps(
    region_dt, combined_dt,
    by.x = c("CHROM", "start", "end"),
    by.y = c("chrom", "start", "end"),
    nomatch = 0, type = "any"
  )
  return(annotated_dt)
}

process_dxy_annotation <- function(dxy_tro_lin, region_rects2) {
  annotated_dxy_lin <- annotate_overlaps(dxy_tro_lin, region_rects2) %>%
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
  
  setDT(annotated_dxy_lin)
  x <- copy(annotated_dxy_lin)[, `:=`(qstart = start, qend = end)]
  setkey(x, CHROM, qstart, qend)
  setkey(y, CHROM, tstart, tend)
  
  ov_dxy <- foverlaps(
    x = x, y = y,
    by.x = c("CHROM", "qstart", "qend"),
    by.y = c("CHROM", "tstart", "tend"),
    type = "any", nomatch = 0L
  )
  ov_dxy[, `:=`(
    ov_len = pmax(0, pmin(qend, tend) - pmax(qstart, tstart)),
    int_len = pmax(0, qend - qstart)
  )]
  ov_dxy[, prop := fifelse(int_len > 0, ov_len / int_len, 0)]
  ov_dxy_best <- as.data.frame(ov_dxy) %>%
    dplyr::select(CHROM, start, end, prop) %>%
    dplyr::group_by(CHROM, start, end) %>%
    dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  anno_dxy <- as.data.frame(annotated_dxy_lin) %>%
    dplyr::left_join(ov_dxy_best, by = c("CHROM", "start", "end")) %>%
    dplyr::mutate(
      hdstatus = ifelse(is.na(prop) | prop < 0.5, "non‑HDR", "HDR"),
      hdstatus = factor(hdstatus, levels = c("non‑HDR", "HDR")),
      stat_type = "Dxy"
    )
  return(anno_dxy)
}

dxy_files <- list(
  Hw3  = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Hw3__dxy.txt"),
  Af   = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Af__dxy.txt"),
  HC   = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_HC__dxy.txt"),
  Hw1  = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Hw1__dxy.txt"),
  Mic1 = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Mic1__dxy.txt"),
  Mic2 = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Mic2__dxy.txt"),
  Tw1  = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Tw1__dxy.txt"),
  Tw3  = readr::read_delim("../../processed_data/Dxy_LAC/results/LAC_Tw3__dxy.txt")
)

anno_list <- lapply(names(dxy_files), function(name) {
  df <- dxy_files[[name]] %>%
    dplyr::rename(
      chrom        = chromosome,
      window_start = window_pos_1,
      window_stop  = window_pos_2,
      stat         = avg_dxy
    ) %>%
    dplyr::mutate(comp_lin = name)
  process_dxy_annotation(df, region_rects2)
})

anno_dxy_all <- dplyr::bind_rows(anno_list)
anno_dxy_all <- anno_dxy_all %>%
  mutate(comp_lin = factor(comp_lin, levels = c("Hw3","Af","Hw1","HC",
                                                "Mic2","Mic1","Tw3","Tw1")))

dxybox_hw3 <- ggplot(anno_dxy_all) +
  geom_jitter(aes(y = hdstatus, x = stat, 
                  colour = hdstatus), 
              size = 0.5,
              alpha = 0.2) +
  geom_boxplot(aes(y = hdstatus, x = stat),
               colour = "grey40",
               outliers = FALSE,
               width = 0.2,
               alpha = 0.6) +
  facet_grid(comp_lin ~ region) +
  theme_bw(base_family = "Helvetica") +
  ylab("") +
  xlab(expression(D[xy])) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_colour_manual(values = c("non‑HDR" = "black", "HDR" = "#0719BC")) +
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

ggsave(dxybox_hw3,
       filename = "../../figures/FigureS32_Dxy.png",
       width = 7, height = 7,
       units = "in", dpi = 600)

nonhdr_dxy_overall <- anno_dxy_all %>%
  dplyr::filter(hdstatus == "non‑HDR") %>%
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
    comp_lin = factor(comp_lin, levels = c("Hw3","Af","Hw1","HC",
                                           "Mic2","Mic1","Tw3","Tw1"))
  ) %>%
  dplyr::arrange(comp_lin)

Ct_nonhdr_dxy_overall<-nonhdr_dxy_overall %>% 
  rename(mean_dxy_nonHDR=mean_dxy,
         group=comp_lin) %>% 
  select(group,mean_dxy_nonHDR) %>% 
  mutate(group = paste0(group, " and LAC"),
         mean_dxy_nonHDR = round(mean_dxy_nonHDR, 3)) %>% 
  dplyr::mutate(species = "C. tropicalis")


## Cb ##
files <- c(
  "../../data/meta_data_Ce_Cb/Cb/Cb_Dxy/Tropical_AD__dxy.txt",
  "../../data/meta_data_Ce_Cb/Cb/Cb_Dxy/Tropical_KD__dxy.txt",
  "../../data/meta_data_Ce_Cb/Cb/Cb_Dxy/Tropical_TD1__dxy.txt",
  "../../data/meta_data_Ce_Cb/Cb/Cb_Dxy/Tropical_Temperate__dxy.txt",
  "../../data/meta_data_Ce_Cb/Cb/Cb_Dxy/Tropical_TH__dxy.txt"
)

all_dxy_raw <- map_dfr(files, ~ {
  df <- read_tsv(.x, show_col_types = FALSE)
  df$group <- gsub(".*Tropical_(.*)__dxy.*", "\\1", .x)
  df
})

genome_domain_raw<-read.csv("../../data/meta_data_Ce_Cb/Cb/Cb_bounds_df.csv",
                            header = TRUE)

genome_domain <- genome_domain_raw %>%
  mutate(
    category = case_when(
      grepl("tip$", sub_region)    ~ "Tip",
      grepl("arm$", sub_region)    ~ "Arm",
      sub_region == "center"       ~ "Center",
      TRUE                         ~ NA_character_
    ),
    xmin = start / 1e6,
    xmax = stop  / 1e6
  ) %>%
  filter(!is.na(category)) %>%
  mutate(category = factor(category, levels = c("Tip", "Arm", "Center")))

all_dxy <- all_dxy_raw %>%
  mutate(
    x = (window_pos_1 + window_pos_2) / 2,
    x_mb = x / 1e6,
    group = as.factor(group)
  ) %>% 
  filter(!(chromosome %in% "MtDNA"))

chrom_levels <- unique(all_dxy$chromosome)
all_dxy <- all_dxy %>% mutate(chrom = factor(chromosome, levels = chrom_levels))
genome_domain <- genome_domain %>% mutate(chrom = factor(chrom, levels = chrom_levels))

hdrs <- readr::read_tsv("../../data/meta_data_Ce_Cb/Cb/HDR_CB_allStrain_5kbclust_20250930.tsv",
                        show_col_types = FALSE)

getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k = 1
    j = 1
    while (k == 1) {
      checkIntersect <- temp %>%
        dplyr::arrange(CHROM, minStart) %>%
        dplyr::mutate(check = ifelse(dplyr::lead(minStart) <= maxEnd, TRUE, FALSE)) %>%
        dplyr::mutate(check = ifelse(is.na(check), FALSE, check))
      
      if (nrow(checkIntersect %>% dplyr::filter(check == TRUE)) == 0) {
        k = 0
      } else {
        temp <- checkIntersect %>%
          dplyr::mutate(gid = data.table::rleid(check)) %>%
          dplyr::mutate(gid = ifelse((check == FALSE | is.na(check)) & dplyr::lag(check) == TRUE,
                                     dplyr::lag(gid), gid))
        
        collapse <- temp %>%
          dplyr::filter(check == TRUE | (check == FALSE & dplyr::lag(check) == TRUE)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart = min(minStart)) %>%
          dplyr::mutate(newEnd = max(maxEnd)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid, .keep_all = TRUE) %>%
          dplyr::mutate(minStart = newStart, maxEnd = newEnd) %>%
          dplyr::select(-newEnd, -newStart)
        
        retain <- temp %>%
          dplyr::filter(check == FALSE & dplyr::lag(check) == FALSE)
        
        temp <- rbind(collapse, retain) %>%
          dplyr::select(-gid, -check)
        
        j = j + 1
      }
    }
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

collapsed_tropical <- plyr::ldply(
  getRegFreq(
    hdrs %>%
      dplyr::filter(source == "QX1410") %>%
      dplyr::group_split(CHROM)
  ),
  data.frame
) %>%
  dplyr::mutate(divSize = maxEnd - minStart) %>%
  dplyr::select(-STRAIN)

all_dxy_overlap <- all_dxy %>%
  dplyr::transmute(
    CHROM = chromosome,
    start = window_pos_1 / 1e6,
    end   = window_pos_2 / 1e6,
    x,
    x_mb,
    group,
    avg_dxy
  )

setDT(all_dxy_overlap)
setDT(collapsed_tropical)

x_dt <- copy(all_dxy_overlap)[, rowid := .I]
y_dt <- copy(collapsed_tropical)[
  , `:=`(tstart = minStart / 1e6,
         tend   = maxEnd / 1e6)
][, .(CHROM, tstart, tend)]

setkey(x_dt, CHROM, start, end)
setkey(y_dt, CHROM, tstart, tend)

ov <- foverlaps(
  x = x_dt,
  y = y_dt,
  by.x = c("CHROM", "start", "end"),
  by.y = c("CHROM", "tstart", "tend"),
  type = "any",
  nomatch = 0L
)

ov[, `:=`(
  ov_len = pmax(0, pmin(end, tend) - pmax(start, tstart)),
  int_len = pmax(0, end - start)
)]

ov[, prop := fifelse(int_len > 0, ov_len / int_len, 0)]

ov_best <- as.data.frame(ov) %>%
  dplyr::select(rowid, prop) %>%
  dplyr::group_by(rowid) %>%
  dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

all_dxy_annotated <- as.data.frame(x_dt) %>%
  dplyr::left_join(ov_best, by = "rowid") %>%
  dplyr::mutate(
    hdstatus = ifelse(is.na(prop) | prop < 0.5, "non-HDR", "HDR")
  )

mean_nonHDR_by_group_raw <- all_dxy_annotated %>%
  dplyr::filter(hdstatus == "non-HDR") %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    mean_dxy_nonHDR = mean(avg_dxy, na.rm = TRUE),
    n_windows_nonHDR = dplyr::n(),
    n_na = sum(is.na(avg_dxy))
  ) %>%
  dplyr::arrange(desc(mean_dxy_nonHDR))

mean_nonHDR_by_group <- mean_nonHDR_by_group_raw %>%
  dplyr::select(group, mean_dxy_nonHDR) %>% 
  mutate(group = paste0(group, " and Tropical"),
         mean_dxy_nonHDR = round(mean_dxy_nonHDR, 3)) %>% 
  dplyr::mutate(species = "C. briggsae")

output<-rbind(Ct_nonhdr_dxy_overall,mean_nonHDR_by_group)
write_csv(output,"../../tables/TableS9_Ct_and_Cb_Dxy_nonHDRs.csv")

