rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(data.table)

domains_raw <- readr::read_tsv("../../data/ct_ch_domains.tsv")

domains_wide <- domains_raw %>%
  dplyr::filter(Location %in% c("left_arm", "right_arm")) %>%
  dplyr::mutate(
    arm = case_when(
      Location == "left_arm"  ~ "left",
      Location == "right_arm" ~ "right",
      TRUE                     ~ NA_character_
    )
  ) %>%
  dplyr::select(CHROM, arm, start, end) %>%
  pivot_wider(
    names_from  = arm,
    values_from = c(start, end),
    names_glue  = "{arm}_{.value}"
  ) %>%
  ungroup()

region_rects <- domains_wide %>%
  pivot_longer(
    cols      = c(left_start, right_end),
    names_to  = "region_side",
    values_to = "x"
  ) %>%
  dplyr::mutate(
    region = "Tip",
    xmin   = case_when(
      region_side == "left_start" ~ 0,
      region_side == "right_end"  ~ x,
      TRUE                         ~ NA_real_
    ),
    xmax = case_when(
      region_side == "left_start" ~ x,
      region_side == "right_end"  ~ Inf,
      TRUE                         ~ NA_real_
    )
  ) %>%
  dplyr::select(CHROM, region, xmin, xmax) %>%
  bind_rows(
    domains_wide %>%
      dplyr::mutate(region = "Arm") %>%
      transmute(CHROM, region, xmin = left_start, xmax = left_end),
    domains_wide %>%
      dplyr::mutate(region = "Arm") %>%
      transmute(CHROM, region, xmin = right_start, xmax = right_end),
    domains_wide %>%
      dplyr::mutate(region = "Center") %>%
      transmute(CHROM, region, xmin = left_end, xmax = right_start)
  ) %>%
  mutate(
    ymin = -Inf,
    ymax = Inf,
    xmin = xmin / 1e6,
    xmax = xmax / 1e6
  )

region_rects2 <- region_rects %>%
  arrange(region, CHROM, xmin) %>%
  dplyr::mutate(
    side = case_when(
      region == "Center" ~ "Center",
      row_number() %% 2 == 0 ~ "Right",
      TRUE ~ "Left"
    )
  )

hdrs_ct <- readr::read_tsv("../../tables/TableS6_HDR_CT_allStrain_5kbclust_20251201.tsv")

y_hdr <- hdrs_ct %>%
  dplyr::mutate(
    tstart = minStart / 1e6,
    tend   = maxEnd   / 1e6
  ) %>%
  dplyr::select(CHROM, tstart, tend) %>%
  distinct() %>%
  as.data.table()

setkey(y_hdr, CHROM, tstart, tend)

annotate_overlaps <- function(combined_df, region_rects) {
  combined_dt <- as.data.table(combined_df %>%
                                 dplyr::mutate(
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
    by.x    = c("CHROM", "start", "end"),
    by.y    = c("chrom", "start", "end"),
    nomatch = 0, type = "any"
  )
  return(annotated_dt)
}

compute_hdr_status <- function(annotated_dt) {
  x_dt <- copy(annotated_dt)[, `:=`(qstart = start, qend = end)]
  setkey(x_dt, CHROM, qstart, qend)
  setkey(y_hdr, CHROM, tstart, tend)
  ov_dt <- foverlaps(
    x = x_dt,
    y = y_hdr,
    by.x = c("CHROM", "qstart", "qend"),
    by.y = c("CHROM", "tstart", "tend"),
    type = "any", nomatch = 0L
  )
  ov_dt[, `:=`(
    ov_len  = pmax(0, pmin(qend, tend) - pmax(qstart, tstart)),
    int_len = pmax(0, qend - qstart)
  )]
  ov_dt[, prop := ov_len / int_len]
  ov_dt[is.na(prop) | is.infinite(prop), prop := 0]
  ov_best <- as.data.frame(ov_dt) %>%
    dplyr::select(CHROM, start, end, prop) %>%
    dplyr::group_by(CHROM, start, end) %>%
    slice_max(prop, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  annotated <- as.data.frame(annotated_dt) %>%
    dplyr::left_join(ov_best, by = c("CHROM", "start", "end")) %>%
    dplyr::mutate(
      hdstatus = case_when(
        is.na(prop) | prop < 0.5 ~ "non_HDR",
        TRUE                     ~ "HDR"
      )
    )
  return(annotated)
}

process_stat_annotation <- function(df, stat_type) {
  df_filt <- df %>% filter(stat_type == !!stat_type)
  annotated_dt <- annotate_overlaps(df_filt, region_rects2) %>%
    dplyr::mutate(
      domain = case_when(
        side == "Center" ~ "center",
        TRUE              ~ "arm"
      )
    ) %>%
    dplyr::filter(region != "Tip")
  annotated_df <- compute_hdr_status(annotated_dt)
  return(annotated_df)
}

stat_calc_hd <- function(path, region_name, stat_type) {
  file_path <- file.path(path, "chromosome_windows_diversity.csv")
  df <- read_csv(file_path, show_col_types = FALSE)
  annotated <- process_stat_annotation(df, stat_type)
  annotated <- annotated %>%
    dplyr::mutate(
      chrom_label = case_when(
        CHROM %in% c("I","II","III","IV","V") ~ "Autosomes",
        TRUE                                           ~ "X"
      )
    )
  summary_df <- annotated %>%
    dplyr::group_by(chrom_label, domain, hdstatus) %>%
    dplyr::summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(stat = stat_type, region = region_name)
  return(summary_df)
}

calc_all_stats_hd <- function(region_paths, stat_types = c("pi", "theta", "d")) {
  results <- purrr::map_df(stat_types, function(st) {
    purrr::imap_dfr(region_paths, function(path, region_name) {
      stat_calc_hd(path, region_name, st)
    })
  })
  return(results)
}

region_paths <- list(
  All            = "../../processed_data/pi_theta_d",
  Africa         = "../../processed_data/pi_theta_d_geo/Africa/",
  Caribbean      = "../../processed_data/pi_theta_d_geo/Caribbean/",
  `Central America` = "../../processed_data/pi_theta_d_geo/Central_America/",
  Hawaii         = "../../processed_data/pi_theta_d_geo/Hawaii/",
  Micronesia     = "../../processed_data/pi_theta_d_geo/Micronesia/",
  `South America`   = "../../processed_data/pi_theta_d_geo/South_America/",
  Taiwan         = "../../processed_data/pi_theta_d_geo/Taiwan/"
)

all_results <- calc_all_stats_hd(region_paths)

make_wide_table <- function(results_df, stat_type) {
  res <- results_df %>%
    dplyr::filter(stat == stat_type) %>%
    unite(col = "chrom_domain_hd", c(chrom_label, domain, hdstatus), sep = "_") %>%
    pivot_wider(
      id_cols     = region,
      names_from  = chrom_domain_hd,
      values_from = mean_value
    )
  res <- res %>%
    mutate(
      `fold change Autosomes HDR arm/center`      = `Autosomes_arm_HDR`     / `Autosomes_center_HDR`,
      `fold change ChromX HDR arm/center`         = `X_arm_HDR`             / `X_center_HDR`,
      `fold change Autosomes non_HDR arm/center`  = `Autosomes_arm_non_HDR` / `Autosomes_center_non_HDR`,
      `fold change ChromX non_HDR arm/center`     = `X_arm_non_HDR`         / `X_center_non_HDR`
    )
  res <- res %>%
    dplyr::rename(Region = region) %>%
    dplyr::mutate(Stat = stat_type)
  return(res)
}

round_final_table <- function(df) {
  fold_cols    <- grep("^fold change", names(df), value = TRUE)
  numeric_cols <- names(df)[sapply(df, is.numeric)]
  base_cols    <- setdiff(numeric_cols, fold_cols)
  
  dynamic_round <- function(x, base_digits) {
    sapply(x, function(value) {
      if (is.na(value)) return(NA_real_)
      digits  <- base_digits
      rounded <- round(value, digits)
      while (abs(rounded) == 0 && abs(value) > 0 && digits < 10) {
        digits  <- digits + 1
        rounded <- round(value, digits)
      }
      rounded
    })
  }
  
  if (length(base_cols) > 0) {
    df <- df %>% dplyr::mutate(across(all_of(base_cols), ~ dynamic_round(.x, 4)))
  }
  if (length(fold_cols) > 0) {
    df <- df %>% dplyr::mutate(across(all_of(fold_cols), ~ dynamic_round(.x, 2)))
  }
  df
}

wide_pi_table <- make_wide_table(all_results, "pi") %>% 
  round_final_table()
wide_theta_table <- make_wide_table(all_results, "theta") %>% 
  round_final_table()
wide_d_table <- make_wide_table(all_results, "d") %>% 
  round_final_table()

fold_cols <- grep("fold", colnames(wide_d_table), value = TRUE)
wide_d_table <- wide_d_table %>%
  dplyr::mutate(across(all_of(fold_cols), ~ "-"))


wide_pi_table_nonHDRs<-wide_pi_table %>% 
  dplyr::select(Region, Autosomes_arm_non_HDR,Autosomes_center_non_HDR,
                X_arm_non_HDR,X_center_non_HDR,
                `fold change Autosomes non_HDR arm/center`,
                `fold change ChromX non_HDR arm/center`,
                Stat)

wide_theta_table_nonHDRs<-wide_theta_table %>% 
  dplyr::select(Region, Autosomes_arm_non_HDR,Autosomes_center_non_HDR,
                X_arm_non_HDR,X_center_non_HDR,
                `fold change Autosomes non_HDR arm/center`,
                `fold change ChromX non_HDR arm/center`,
                Stat)

wide_d_table_nonHDRs<-wide_d_table %>% 
  dplyr::select(Region, Autosomes_arm_non_HDR,Autosomes_center_non_HDR,
                X_arm_non_HDR,X_center_non_HDR,
                `fold change Autosomes non_HDR arm/center`,
                `fold change ChromX non_HDR arm/center`,
                Stat)

output_table<-rbind(wide_pi_table_nonHDRs,
                    wide_theta_table_nonHDRs,
                    wide_d_table_nonHDRs)

write.table(output_table,
            "../../tables/TableS8_Geo_pi_theta_d_Autosomes_X_nonHDR.tsv", 
            quote = FALSE, row.names = FALSE,
            sep = '\t')

