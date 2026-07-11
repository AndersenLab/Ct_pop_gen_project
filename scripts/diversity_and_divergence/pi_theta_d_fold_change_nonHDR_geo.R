rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

format_4_or_first_nonzero <- function(x) {
  purrr::map_chr(
    x,
    function(value) {
      if (is.na(value)) {
        return(NA_character_)
      }
      
      if (value == 0) {
        return("0")
      }
      
      abs_value <- abs(value)
      
      if (abs_value >= 0.0001) {
        return(
          format(
            round(value, 4),
            nsmall = 4,
            scientific = FALSE,
            trim = TRUE
          )
        )
      }
      
      digits_to_keep <- ceiling(-log10(abs_value))
      value_to_write <- trunc(value * 10^digits_to_keep) / 10^digits_to_keep
      
      format(
        value_to_write,
        nsmall = digits_to_keep,
        scientific = FALSE,
        trim = TRUE
      )
    }
  )
}

format_fold_change <- function(x) {
  purrr::map_chr(
    x,
    function(value) {
      if (is.na(value) | is.infinite(value)) {
        return(NA_character_)
      }
      
      format(
        round(value, 2),
        nsmall = 2,
        scientific = FALSE,
        trim = TRUE
      )
    }
  )
}

read_pixy_stat <- function(file_path, stat_type, stat_col) {
  readr::read_tsv(file_path, show_col_types = FALSE) %>%
    dplyr::filter(chromosome != "MtDNA") %>%
    dplyr::transmute(
      Region = dplyr::case_when(
        pop == "GLOBAL" ~ "All",
        TRUE            ~ pop
      ),
      CHROM = chromosome,
      window_start = window_pos_1,
      window_stop = window_pos_2,
      stat_type = stat_type,
      stat = .data[[stat_col]]
    )
}

pixy_df <- dplyr::bind_rows(
  read_pixy_stat(
    "../../processed_data/pi_theta_d/results/Ct_GLOBAL_pi.txt",
    "pi",
    "avg_pi"
  ),
  read_pixy_stat(
    "../../processed_data/pi_theta_d/results/Ct_geo_pi.txt",
    "pi",
    "avg_pi"
  ),
  read_pixy_stat(
    "../../processed_data/pi_theta_d/results/Ct_GLOBAL_watterson_theta.txt",
    "theta",
    "avg_watterson_theta"
  ),
  read_pixy_stat(
    "../../processed_data/pi_theta_d/results/Ct_geo_watterson_theta.txt",
    "theta",
    "avg_watterson_theta"
  ),
  read_pixy_stat(
    "../../processed_data/pi_theta_d/results/Ct_GLOBAL_tajima_d.txt",
    "d",
    "tajima_d"
  ),
  read_pixy_stat(
    "../../processed_data/pi_theta_d/results/Ct_geo_tajima_d.txt",
    "d",
    "tajima_d"
  )
)

domains_raw <- readr::read_tsv("../../data/ct_ch_domains.tsv", show_col_types = FALSE)

domains_wide <- domains_raw %>%
  dplyr::filter(Location %in% c("left_arm", "right_arm")) %>%
  dplyr::mutate(
    arm = dplyr::case_when(
      Location == "left_arm"  ~ "left",
      Location == "right_arm" ~ "right",
      TRUE                    ~ NA_character_
    )
  ) %>%
  dplyr::select(CHROM, arm, start, end) %>%
  tidyr::pivot_wider(
    names_from  = arm,
    values_from = c(start, end),
    names_glue  = "{arm}_{.value}"
  ) %>%
  dplyr::ungroup()

region_rects <- dplyr::bind_rows(
  domains_wide %>%
    dplyr::transmute(
      CHROM,
      region = "Tip",
      xmin = 0,
      xmax = left_start
    ),
  domains_wide %>%
    dplyr::transmute(
      CHROM,
      region = "Arm",
      xmin = left_start,
      xmax = left_end
    ),
  domains_wide %>%
    dplyr::transmute(
      CHROM,
      region = "Center",
      xmin = left_end,
      xmax = right_start
    ),
  domains_wide %>%
    dplyr::transmute(
      CHROM,
      region = "Arm",
      xmin = right_start,
      xmax = right_end
    ),
  domains_wide %>%
    dplyr::transmute(
      CHROM,
      region = "Tip",
      xmin = right_end,
      xmax = Inf
    )
) %>%
  dplyr::mutate(
    xmin = xmin / 1e6,
    xmax = xmax / 1e6,
    domain = dplyr::case_when(
      region == "Center" ~ "center",
      region == "Arm"    ~ "arm",
      TRUE               ~ NA_character_
    )
  ) %>%
  dplyr::filter(region != "Tip")

hdrs_ct <- readr::read_tsv(
  "../../tables/TableS6_HDR_CT_allStrain_5kbclust_20251201.tsv",
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    CHROM,
    hdr_start = minStart / 1e6,
    hdr_end = maxEnd / 1e6
  ) %>%
  dplyr::distinct()

annotate_domain <- function(df) {
  df %>%
    dplyr::mutate(
      win_start = window_start / 1e6,
      win_end = window_stop / 1e6
    ) %>%
    dplyr::inner_join(region_rects, by = "CHROM") %>%
    dplyr::filter(
      win_start <= xmax,
      win_end >= xmin
    )
}

annotate_hdr_status <- function(df) {
  df %>%
    dplyr::left_join(hdrs_ct, by = "CHROM") %>%
    dplyr::mutate(
      ov_len = dplyr::case_when(
        is.na(hdr_start) | is.na(hdr_end) ~ 0,
        TRUE ~ pmax(0, pmin(win_end, hdr_end) - pmax(win_start, hdr_start))
      ),
      win_len = pmax(0, win_end - win_start),
      prop = dplyr::case_when(
        win_len == 0 ~ 0,
        TRUE ~ ov_len / win_len
      )
    ) %>%
    dplyr::group_by(
      Region,
      CHROM,
      window_start,
      window_stop,
      stat_type,
      stat,
      domain
    ) %>%
    dplyr::summarise(
      prop = max(prop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      prop = dplyr::case_when(
        is.na(prop) | is.infinite(prop) ~ 0,
        TRUE ~ prop
      ),
      hdstatus = dplyr::case_when(
        prop >= 0.5 ~ "HDR",
        TRUE        ~ "non_HDR"
      )
    )
}

stat_calc_hd <- function(df, stat_to_use, stat_label) {
  df %>%
    dplyr::filter(stat_type == stat_to_use) %>%
    annotate_domain() %>%
    annotate_hdr_status() %>%
    dplyr::mutate(
      chrom_label = dplyr::case_when(
        CHROM %in% c("I", "II", "III", "IV", "V") ~ "Autosomes",
        TRUE                                      ~ "X"
      )
    ) %>%
    dplyr::group_by(Region, chrom_label, domain, hdstatus) %>%
    dplyr::summarise(
      mean_value = mean(stat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Stat = stat_label)
}

all_results <- dplyr::bind_rows(
  stat_calc_hd(pixy_df, "pi", "pi"),
  stat_calc_hd(pixy_df, "theta", "theta"),
  stat_calc_hd(pixy_df, "d", "Tajima's D")
)

make_wide_table <- function(results_df, stat_label) {
  results_df %>%
    dplyr::filter(Stat == stat_label) %>%
    tidyr::unite(
      col = "chrom_domain_hd",
      c(chrom_label, domain, hdstatus),
      sep = "_"
    ) %>%
    tidyr::pivot_wider(
      id_cols     = c(Region, Stat),
      names_from  = chrom_domain_hd,
      values_from = mean_value
    ) %>%
    dplyr::mutate(
      `fold change Autosomes HDR arm/center`      = `Autosomes_arm_HDR`     / `Autosomes_center_HDR`,
      `fold change ChromX HDR arm/center`         = `X_arm_HDR`             / `X_center_HDR`,
      `fold change Autosomes non_HDR arm/center`  = `Autosomes_arm_non_HDR` / `Autosomes_center_non_HDR`,
      `fold change ChromX non_HDR arm/center`     = `X_arm_non_HDR`         / `X_center_non_HDR`
    )
}

wide_pi_table <- make_wide_table(all_results, "pi")

wide_theta_table <- make_wide_table(all_results, "theta")

wide_d_table <- make_wide_table(all_results, "Tajima's D")

output_table <- dplyr::bind_rows(
  wide_pi_table,
  wide_theta_table,
  wide_d_table
) %>%
  dplyr::select(
    Region,
    Autosomes_arm_non_HDR,
    Autosomes_center_non_HDR,
    X_arm_non_HDR,
    X_center_non_HDR,
    `fold change Autosomes non_HDR arm/center`,
    `fold change ChromX non_HDR arm/center`,
    Stat
  ) %>%
  dplyr::filter(
    Region %in% c(
      "All", "Africa", "Caribbean", "Central America",
      "Hawaii", "Micronesia", "South America", "Taiwan"
    )
  ) %>%
  dplyr::mutate(
    Stat = factor(
      Stat,
      levels = c("pi", "theta", "Tajima's D")
    ),
    Region = factor(
      Region,
      levels = c(
        "All", "Africa", "Caribbean", "Central America",
        "Hawaii", "Micronesia", "South America", "Taiwan"
      )
    )
  ) %>%
  dplyr::arrange(Stat, Region) %>%
  dplyr::mutate(
    Region = as.character(Region),
    Stat = as.character(Stat),
    Autosomes_arm_non_HDR = format_4_or_first_nonzero(Autosomes_arm_non_HDR),
    Autosomes_center_non_HDR = format_4_or_first_nonzero(Autosomes_center_non_HDR),
    X_arm_non_HDR = format_4_or_first_nonzero(X_arm_non_HDR),
    X_center_non_HDR = format_4_or_first_nonzero(X_center_non_HDR),
    `fold change Autosomes non_HDR arm/center` = dplyr::case_when(
      Stat == "Tajima's D" ~ "-",
      TRUE ~ format_fold_change(`fold change Autosomes non_HDR arm/center`)
    ),
    `fold change ChromX non_HDR arm/center` = dplyr::case_when(
      Stat == "Tajima's D" ~ "-",
      TRUE ~ format_fold_change(`fold change ChromX non_HDR arm/center`)
    )
  )

write.table(
  output_table,
  "../../tables/TableS8_Geo_pi_theta_d_Autosomes_X_nonHDR.tsv",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
