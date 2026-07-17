rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)
library(purrr)


## Helpers
read_pixy_file <- function(file_path) {
  readr::read_tsv(
    file_path,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("NA", "NaN", "")
  )
}

make_old_format <- function(pixy_df, stat_type, stat_col) {
  pixy_df %>%
    dplyr::filter(chromosome != "MtDNA") %>%
    dplyr::transmute(
      chrom = chromosome,
      x = (window_pos_1 + window_pos_2) / 2,
      window_start = window_pos_1,
      window_stop = window_pos_2,
      stat_type = stat_type,
      stat = .data[[stat_col]]
    )
}

write_old_format <- function(chromosome_windows_diversity, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  readr::write_csv(
    chromosome_windows_diversity,
    file.path(out_dir, "chromosome_windows_diversity.csv"),
    na = ""
  )
  
  windows <- chromosome_windows_diversity %>%
    dplyr::distinct(chrom, x, window_start, window_stop) %>%
    dplyr::arrange(chrom, window_start)
  
  readr::write_csv(
    windows,
    file.path(out_dir, "windows.csv"),
    na = ""
  )
}

calculate_diversity <- function(file_path, diversity_type) {
  readr::read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::filter(stat_type == diversity_type) %>%
    dplyr::pull(stat) %>%
    mean(na.rm = TRUE)
}

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

## Convert Pixy GLOBAL output to old format
out_dir <- "../../processed_data/pi_theta_d"

global_pi <- read_pixy_file("../../processed_data/pi_theta_d/results/Ct_GLOBAL_pi.txt") %>%
  make_old_format(stat_type = "pi", stat_col = "avg_pi")

global_theta <- read_pixy_file("../../processed_data/pi_theta_d/results/Ct_GLOBAL_watterson_theta.txt") %>%
  make_old_format(stat_type = "theta", stat_col = "avg_watterson_theta")

global_d <- read_pixy_file("../../processed_data/pi_theta_d/results/Ct_GLOBAL_tajima_d.txt") %>%
  make_old_format(stat_type = "d", stat_col = "tajima_d")

global_diversity <- dplyr::bind_rows(global_theta, global_pi, global_d) %>%
  dplyr::arrange(chrom, window_start, stat_type)

write_old_format(global_diversity, out_dir)

## Convert Pixy geo output to old format
out_dir <- "../../processed_data/pi_theta_d_geo"

geo_pi <- read_pixy_file("../../processed_data/pi_theta_d/results/Ct_geo_pi.txt") %>%
  dplyr::mutate(region = gsub(" ", "_", pop))

geo_theta <- read_pixy_file("../../processed_data/pi_theta_d/results/Ct_geo_watterson_theta.txt") %>%
  dplyr::mutate(region = gsub(" ", "_", pop))

geo_d <- read_pixy_file("../../processed_data/pi_theta_d/results/Ct_geo_tajima_d.txt") %>%
  dplyr::mutate(region = gsub(" ", "_", pop))

regions <- geo_pi %>%
  dplyr::pull(region) %>%
  unique() %>%
  sort()

for (custom_region in regions) {
  region_pi <- geo_pi %>%
    dplyr::filter(region == custom_region) %>%
    make_old_format(stat_type = "pi", stat_col = "avg_pi")
  
  region_theta <- geo_theta %>%
    dplyr::filter(region == custom_region) %>%
    make_old_format(stat_type = "theta", stat_col = "avg_watterson_theta")
  
  region_d <- geo_d %>%
    dplyr::filter(region == custom_region) %>%
    make_old_format(stat_type = "d", stat_col = "tajima_d")
  
  region_diversity <- dplyr::bind_rows(region_theta, region_pi, region_d) %>%
    dplyr::arrange(chrom, window_start, stat_type)
  
  write_old_format(
    region_diversity,
    file.path(out_dir, custom_region)
  )
}

## Make Table S4 summary
species_summary <- tibble::tibble(
  Region = "All",
  pi = calculate_diversity(
    "../../processed_data/pi_theta_d/chromosome_windows_diversity.csv",
    "pi"
  ),
  theta = calculate_diversity(
    "../../processed_data/pi_theta_d/chromosome_windows_diversity.csv",
    "theta"
  ),
  d = calculate_diversity(
    "../../processed_data/pi_theta_d/chromosome_windows_diversity.csv",
    "d"
  )
)

geo_summary <- purrr::map_dfr(
  regions,
  function(custom_region) {
    tibble::tibble(
      Region = custom_region,
      pi = calculate_diversity(
        file.path("../../processed_data/pi_theta_d_geo", custom_region, "chromosome_windows_diversity.csv"),
        "pi"
      ),
      theta = calculate_diversity(
        file.path("../../processed_data/pi_theta_d_geo", custom_region, "chromosome_windows_diversity.csv"),
        "theta"
      ),
      d = calculate_diversity(
        file.path("../../processed_data/pi_theta_d_geo", custom_region, "chromosome_windows_diversity.csv"),
        "d"
      )
    )
  }
)

geo_diversity_result <- dplyr::bind_rows(species_summary, geo_summary) %>%
  dplyr::mutate(
    pi = format_4_or_first_nonzero(pi),
    theta = format_4_or_first_nonzero(theta),
    d = format_4_or_first_nonzero(d)
  ) %>%
  dplyr::select(Region, pi, theta, d) %>%
  dplyr::rename(
    "π value" = pi,
    "θ value" = theta,
    "Tajima's D" = d
  )

geo_freq_col <- readr::read_csv(
  "../../processed_data/geo_info/Ct_isotype_geo_freq.csv",
  show_col_types = FALSE
)

if ("freq" %in% colnames(geo_freq_col) && !("frequency" %in% colnames(geo_freq_col))) {
  geo_freq_col <- geo_freq_col %>%
    dplyr::rename(frequency = freq)
}

geo_freq_col <- geo_freq_col %>%
  dplyr::mutate(
    geo = dplyr::case_when(
      geo == "South America" ~ "South_America",
      geo == "Central America" ~ "Central_America",
      TRUE ~ geo
    )
  )

geo_freq_col <- dplyr::bind_rows(
  geo_freq_col,
  tibble::tibble(geo = "All", frequency = sum(geo_freq_col$frequency))
) %>%
  dplyr::mutate(
    geo = factor(
      geo,
      levels = c(
        "All", "Africa", "Australia", "Caribbean", "Central_America",
        "Hawaii", "South_America", "Indonesia", "Micronesia", "Taiwan"
      )
    )
  ) %>%
  dplyr::arrange(geo) %>%
  dplyr::rename(
    Region = geo,
    "Number of strains" = frequency
  )

all_result <- dplyr::left_join(
  geo_diversity_result,
  geo_freq_col,
  by = "Region"
)

dir.create("../../tables", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  all_result,
  "../../tables/TableS5_geo_pi_theta_d.csv",
  na = ""
)
