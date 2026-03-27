rm(list = ls())

library(extrafont)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
source("../utilities.R")

window_diversity<-read.csv("../../processed_data/pi_theta_exclude_HDRs/chromosome_windows_diversity.csv",
                           header = TRUE,
                           row.names = 1)
window_diversity<-window_diversity %>% 
  filter(!(is.na(stat)) & (stat != 0))

window_diversity_pi<-window_diversity%>%
  dplyr::filter(stat_type == "pi") %>% 
  na.omit()
mean_pi <- mean(window_diversity_pi$stat)
mean_pi

window_diversity_theta<-window_diversity%>%
  dplyr::filter(stat_type == "theta") %>% 
  na.omit()
mean_theta <- mean(window_diversity_theta$stat)
mean_theta

genome_domain_raw<-readr::read_tsv("../../data/ct_ch_domains.tsv")
genome_domain <- genome_domain_raw %>%
  rename(chrom=CHROM) %>% 
  mutate(
    category = case_when(
      grepl("tip$", Location)    ~ "Tip",
      grepl("arm$", Location)    ~ "Arm",
      Location == "center"       ~ "Center",
      TRUE                         ~ NA_character_
    ),
    xmin = start / 1e6,
    xmax = end  / 1e6
  ) %>%
  filter(!is.na(category)) %>%
  mutate(category = factor(category, levels = c("Tip", "Arm", "Center")))

windowed_div_stats_no_d <- function(windows_df){
  diversity <- windows_df %>%
    filter(stat_type %in% c("pi","theta"
                            )) %>%
    mutate(
      stat_type = case_when(
        stat_type == "pi"    ~ "Nucleotide diversity (\u03C0)",
        stat_type == "theta" ~ "Watterson's \u03B8"
      ),
      x_mb = x / 1e6
    )
  
  diversity$stat_type <- factor(
    diversity$stat_type,
    levels = c("Nucleotide diversity (\u03C0)", 
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
      color = "gray30", size = 0.05, alpha = 0.8, shape = 16
    ) +
    geom_smooth(
      data = diversity,
      aes(x = x_mb, y = stat),
      method = "loess", se = FALSE, span = 0.3, color = "lightgray"
    ) +
    facet_grid(stat_type ~ chrom, scales = "free") +
    scale_fill_manual(values = genome_domain_colors) +
    xlab("Physical genome position (Mb)") +
    ylab("Diversity statistic") +
    theme_bw() +
    theme(
      legend.position      = "none",
      panel.grid           = element_blank(),
      strip.background     = element_blank(),
      strip.text           = element_text(face = "bold", size = 7, color = "#525252"),
      axis.title           = element_text(face = "bold", size = 9),
      axis.text            = element_text(size = 6)
    )
}

result_plots_no_d <- windowed_div_stats_no_d(window_diversity)
ggsave("../../figures/FigureS35_pi_theta_nonHDRs.png", 
       plot = result_plots_no_d, 
       width = 7, height = 3, 
       units = "in", dpi = 600)

