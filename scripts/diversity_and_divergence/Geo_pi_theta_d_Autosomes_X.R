rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(purrr)

sub_region_raw <- read.table("../../data/ct_ch_domains.tsv", sep = '\t')

sub_region_df <- sub_region_raw %>%
  dplyr::filter(V1 != "CHROM") %>%
  dplyr::mutate(
    V2 = as.numeric(V2),
    V3 = as.numeric(V3)
  ) %>%
  dplyr::rename(
    chrom = V1,
    region_start = V2,
    region_end = V3,
    subregion = V4
  )

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

read_pixy_file <- function(file_path) {
  readr::read_tsv(file_path, show_col_types = FALSE) %>%
    dplyr::filter(chromosome != "MtDNA") %>%
    dplyr::rename(
      chrom = chromosome,
      window_start = window_pos_1,
      window_stop = window_pos_2
    )
}

make_pixy_stat_table <- function(global_file, geo_file, stat_col, stat_name) {
  global_df <- read_pixy_file(global_file) %>%
    dplyr::mutate(Region = "All")
  
  geo_df <- read_pixy_file(geo_file) %>%
    dplyr::rename(Region = pop) %>%
    dplyr::mutate(
      Region = dplyr::case_when(
        Region == "Central America" ~ "Central America",
        Region == "South America" ~ "South America",
        TRUE ~ Region
      )
    )
  
  dplyr::bind_rows(global_df, geo_df) %>%
    dplyr::filter(
      Region %in% c(
        "All", "Africa", "Caribbean", "Central America",
        "Hawaii", "Micronesia", "South America", "Taiwan"
      )
    ) %>%
    dplyr::transmute(
      Region,
      chrom,
      window_start,
      window_stop,
      stat = as.numeric(.data[[stat_col]]),
      Stat = stat_name
    )
}

div_calc <- function(df) {
  process <- function(df, chrom_filter, chrom_label) {
    df %>%
      dplyr::filter(chrom %in% chrom_filter) %>%
      dplyr::mutate(mid = (window_start + window_stop) / 2) %>%
      dplyr::inner_join(sub_region_df, by = "chrom") %>%
      dplyr::filter(mid >= region_start, mid <= region_end) %>%
      dplyr::mutate(
        domain = dplyr::case_when(
          subregion == "center" ~ "center",
          subregion %in% c("left_arm", "right_arm") ~ "arm",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(domain)) %>%
      dplyr::group_by(Region, Stat, domain) %>%
      dplyr::summarise(mean_value = mean(stat, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(chrom = chrom_label)
  }
  
  autosomes <- process(df, c("I", "II", "III", "IV", "V"), "Autosomes")
  chrX <- process(df, c("X"), "X")
  
  dplyr::bind_rows(autosomes, chrX)
}

pi_table <- make_pixy_stat_table(
  "../../processed_data/pi_theta_d/results/Ct_GLOBAL_pi.txt",
  "../../processed_data/pi_theta_d/results/Ct_geo_pi.txt",
  "avg_pi",
  "pi"
)

theta_table <- make_pixy_stat_table(
  "../../processed_data/pi_theta_d/results/Ct_GLOBAL_watterson_theta.txt",
  "../../processed_data/pi_theta_d/results/Ct_geo_watterson_theta.txt",
  "avg_watterson_theta",
  "theta"
)

d_table <- make_pixy_stat_table(
  "../../processed_data/pi_theta_d/results/Ct_GLOBAL_tajima_d.txt",
  "../../processed_data/pi_theta_d/results/Ct_geo_tajima_d.txt",
  "tajima_d",
  "Tajima's D"
)

merged_table <- dplyr::bind_rows(pi_table, theta_table, d_table)

merged_wide_table <- div_calc(merged_table) %>%
  tidyr::unite(col = "chrom_domain", chrom, domain, sep = "_") %>%
  tidyr::pivot_wider(
    id_cols = c(Region, Stat),
    names_from = chrom_domain,
    values_from = mean_value
  ) %>%
  dplyr::rename(
    `Autosomes arm` = Autosomes_arm,
    `Autosomes center` = Autosomes_center,
    `ChromX arm` = X_arm,
    `ChromX center` = X_center
  )

geo_merged_wide_table <- merged_wide_table %>%
  dplyr::filter(
    Region %in% c(
      "All", "Africa", "Caribbean", "Central America",
      "Hawaii", "Micronesia", "South America", "Taiwan"
    )
  ) %>%
  dplyr::mutate(
    `fold change Autosomes arm/center` = round(`Autosomes arm` / `Autosomes center`, 2),
    `fold change ChromX arm/center` = round(`ChromX arm` / `ChromX center`, 2)
  ) %>%
  dplyr::mutate(
    `fold change Autosomes arm/center` = ifelse(Stat == "Tajima's D", "-", as.character(`fold change Autosomes arm/center`)),
    `fold change ChromX arm/center` = ifelse(Stat == "Tajima's D", "-", as.character(`fold change ChromX arm/center`))
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
    Stat = as.character(Stat),
    Region = as.character(Region),
    `Autosomes arm` = format_4_or_first_nonzero(`Autosomes arm`),
    `Autosomes center` = format_4_or_first_nonzero(`Autosomes center`),
    `ChromX arm` = format_4_or_first_nonzero(`ChromX arm`),
    `ChromX center` = format_4_or_first_nonzero(`ChromX center`)
  ) %>%
  dplyr::relocate(Stat, .after = dplyr::last_col())

write.table(
  geo_merged_wide_table,
  "../../tables/TableS5_Geo_pi_theta_d_Autosomes_X.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = '\t'
)
