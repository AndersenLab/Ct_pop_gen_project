rm(list = ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(readr)

source("../utilities.R")

read_pixy_1to1_4fold <- function() {
  pi_df <- readr::read_tsv(
    "../../processed_data/HDR_stats/pi_theta_d_1to1_4fold_pixy/Ct_1to1_4fold_GLOBAL_pi.txt",
    show_col_types = FALSE
  ) %>%
    dplyr::transmute(
      chrom = chromosome,
      x = (window_pos_1 + window_pos_2) / 2,
      window_start = window_pos_1,
      window_stop = window_pos_2,
      stat_type = "pi",
      stat = avg_pi
    )
  
  theta_df <- readr::read_tsv(
    "../../processed_data/HDR_stats/pi_theta_d_1to1_4fold_pixy/Ct_1to1_4fold_GLOBAL_watterson_theta.txt",
    show_col_types = FALSE
  ) %>%
    dplyr::transmute(
      chrom = chromosome,
      x = (window_pos_1 + window_pos_2) / 2,
      window_start = window_pos_1,
      window_stop = window_pos_2,
      stat_type = "theta",
      stat = avg_watterson_theta
    )
  
  dplyr::bind_rows(pi_df, theta_df) %>%
    dplyr::filter(chrom != "MtDNA")
}

window_diversity <- read_pixy_1to1_4fold() %>%
  dplyr::filter(!(is.na(stat)) & (stat != 0))

### read genome domain
genome_domain_raw <- readr::read_tsv("../../data/ct_ch_domains.tsv", show_col_types = FALSE)

genome_domain <- genome_domain_raw %>%
  dplyr::rename(chrom = CHROM) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      grepl("tip$", Location) ~ "Tip",
      grepl("arm$", Location) ~ "Arm",
      Location == "center" ~ "Center",
      TRUE ~ NA_character_
    ),
    xmin = start / 1e6,
    xmax = end / 1e6
  ) %>%
  dplyr::filter(!is.na(category)) %>%
  dplyr::mutate(category = factor(category, levels = c("Tip", "Arm", "Center")))

#### Plot all windowed diversity stats
windowed_div_stats_no_d <- function(windows_df) {
  diversity <- windows_df %>%
    dplyr::filter(stat_type %in% c("pi", "theta")) %>%
    dplyr::mutate(
      stat_type = dplyr::case_when(
        stat_type == "pi" ~ "Nucleotide diversity (\u03C0)",
        stat_type == "theta" ~ "Watterson's \u03B8"
      ),
      x_mb = x / 1e6
    )
  
  diversity$stat_type <- factor(
    diversity$stat_type,
    levels = c(
      "Nucleotide diversity (\u03C0)",
      "Watterson's \u03B8"
    )
  )
  
  ggplot() +
    geom_rect(
      data = genome_domain,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = category),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    geom_point(
      data = diversity,
      aes(x = x_mb, y = stat),
      color = "gray30",
      size = 0.05,
      alpha = 0.8,
      shape = 16
    ) +
    geom_smooth(
      data = diversity,
      aes(x = x_mb, y = stat),
      method = "loess",
      se = FALSE,
      span = 0.3,
      color = "lightgray"
    ) +
    facet_grid(stat_type ~ chrom, scales = "free") +
    scale_fill_manual(values = genome_domain_colors) +
    xlab("Physical genome position (Mb)") +
    ylab("Diversity statistic") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7, color = "#525252"),
      axis.title = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 6)
    )
}

result_plots_no_d <- windowed_div_stats_no_d(window_diversity)
result_plots_no_d

ggsave(
  "../../figures/FigureS42_plot_pi_theta_1to1_4fold.png",
  plot = result_plots_no_d,
  width = 7,
  height = 3,
  units = "in",
  dpi = 600
)

#### Mean pi in HDR vs non-HDR windows
getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k <- 1
    j <- 1
    
    while (k == 1) {
      checkIntersect <- temp %>%
        dplyr::arrange(CHROM, minStart) %>%
        dplyr::mutate(check = ifelse(lead(minStart) <= maxEnd, TRUE, FALSE)) %>%
        dplyr::mutate(check = ifelse(is.na(check), FALSE, check))
      
      print(nrow(checkIntersect %>% dplyr::filter(check == TRUE)))
      
      if (nrow(checkIntersect %>% dplyr::filter(check == TRUE)) == 0) {
        print("NO MORE INTERSECTS")
        k <- 0
      } else {
        temp <- checkIntersect %>%
          dplyr::mutate(gid = data.table::rleid(check)) %>%
          dplyr::mutate(
            gid = ifelse(
              (check == FALSE | is.na(check)) & lag(check) == TRUE,
              lag(gid),
              gid
            )
          )
        
        collapse <- temp %>%
          dplyr::filter(check == TRUE | (check == FALSE & lag(check) == TRUE)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart = min(minStart)) %>%
          dplyr::mutate(newEnd = max(maxEnd)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid, .keep_all = TRUE) %>%
          dplyr::mutate(minStart = newStart, maxEnd = newEnd) %>%
          dplyr::select(-newEnd, -newStart)
        
        retain <- temp %>%
          dplyr::filter(check == FALSE & dplyr::coalesce(dplyr::lag(check), FALSE) == FALSE)
        
        temp <- rbind(collapse, retain) %>%
          dplyr::select(-gid, -check)
        
        j <- j + 1
      }
    }
    
    print("WRITING TO LIST")
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  
  return(all_collapsed)
}

hdr_regions <- readr::read_tsv(
  "../../tables/TableS7_HDR_CT_allStrain_5kbclust_20251201.tsv",
  show_col_types = FALSE
) %>%
  dplyr::select(CHROM, minStart, maxEnd, STRAIN) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::arrange(CHROM, minStart, maxEnd) %>%
  dplyr::mutate(divSize = maxEnd - minStart) %>%
  dplyr::filter(divSize >= 5e3) %>%
  dplyr::select(CHROM, minStart, maxEnd, STRAIN)

nonoverlap_hdrs <- plyr::ldply(
  getRegFreq(
    hdr_regions %>%
      dplyr::arrange(CHROM, minStart) %>%
      dplyr::group_split(CHROM)
  ),
  data.frame
) %>%
  dplyr::select(-STRAIN) %>%
  dplyr::rename(start = minStart, end = maxEnd) %>%
  dplyr::mutate(divSize = end - start) %>%
  dplyr::arrange(CHROM, start, end)

###### Prepare pi windows only ######
pi_windows <- window_diversity %>%
  dplyr::filter(stat_type == "pi") %>%
  dplyr::mutate(
    rowid = dplyr::row_number(),
    CHROM = chrom,
    start = window_start / 1e6,
    end = window_stop / 1e6
  ) %>%
  dplyr::select(rowid, CHROM, start, end, stat)

###### Prepare HDRs in Mb ######
hdr_overlap <- nonoverlap_hdrs %>%
  dplyr::transmute(
    CHROM = CHROM,
    tstart = start / 1e6,
    tend = end / 1e6
  )

###### Overlap windows with HDRs ######
data.table::setDT(pi_windows)
data.table::setDT(hdr_overlap)

x_dt <- data.table::copy(pi_windows)
y_dt <- data.table::copy(hdr_overlap)

data.table::setkey(x_dt, CHROM, start, end)
data.table::setkey(y_dt, CHROM, tstart, tend)

ov <- data.table::foverlaps(
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

ov[, prop := data.table::fifelse(int_len > 0, ov_len / int_len, 0)]

ov_best <- as.data.frame(ov) %>%
  dplyr::select(rowid, prop) %>%
  dplyr::group_by(rowid) %>%
  dplyr::slice_max(prop, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

###### Add HDR / non-HDR label and calculate mean pi ######
pi_HDR_labeled <- as.data.frame(x_dt) %>%
  dplyr::left_join(ov_best, by = "rowid") %>%
  dplyr::mutate(
    hdstatus = ifelse(is.na(prop) | prop < 0.5, "non-HDR", "HDR")
  )

mean_pi_HDR_nonHDR <- pi_HDR_labeled %>%
  dplyr::group_by(hdstatus) %>%
  dplyr::summarise(
    n_windows = dplyr::n(),
    mean_pi = mean(stat, na.rm = TRUE),
    .groups = "drop"
  )

mean_pi_HDR_nonHDR


