rm(list = ls())

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(data.table)
library(plyr)

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

summary_nonoverlap_hdrs <- nonoverlap_hdrs %>%
  dplyr::summarise(
    n_nonoverlap_HDRs = dplyr::n(),
    mean_length_kb = mean(divSize) / 1e3,
    min_length_kb = min(divSize) / 1e3,
    max_length_kb = max(divSize) / 1e3,
    total_span_bp = sum(divSize),
    total_span_Mb = sum(divSize) / 1e6
  )

print(summary_nonoverlap_hdrs)



TAs <- data.frame(
  TA_id = paste0("TA", seq_len(9)),
  Chr   = c("ChrIII","ChrV","ChrII","ChrII","ChrIII","ChrV","ChrV","ChrII","ChrV"),
  Start = c(1.27e6,12.94e6,8009743,2836704,10426642,1289235,11935356,8077521,1260950),
  End   = c(1.51e6,13.32e6,8756516,8581670,10466762,1751145,13453526,8836894,1735109),
  stringsAsFactors = FALSE
)

norm_chr <- function(x) str_replace(x, "^Chr", "")

TAs2 <- TAs %>%
  dplyr::mutate(Chr_norm = norm_chr(Chr))

HDRs2 <- nonoverlap_hdrs %>%
  dplyr::mutate(
    CHROM_norm = norm_chr(CHROM),
    HDR_hit_id = paste0("nonoverlap_HDR_", dplyr::row_number())
  )



TA_HDR_hits <- TAs2 %>%
  dplyr::inner_join(HDRs2, by = c("Chr_norm" = "CHROM_norm")) %>%
  dplyr::filter(start <= End, end >= Start) %>%
  dplyr::mutate(
    overlap_start = pmax(start, Start),
    overlap_end   = pmin(end, End),
    overlap_bp    = overlap_end - overlap_start + 1
  )

TA_HDR_hits_long <- TA_HDR_hits %>%
  dplyr::select(
    TA_id,
    Chr,
    Start,
    End,
    CHROM,
    start,
    end,
    divSize,
    HDR_hit_id,
    overlap_start,
    overlap_end,
    overlap_bp,
    everything()
  )

TA_HDR_summary <- TA_HDR_hits_long %>%
  dplyr::group_by(TA_id, Chr, Start, End) %>%
  dplyr::summarise(
    n_HDRs = n_distinct(HDR_hit_id),
    HDRs   = paste(unique(HDR_hit_id), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(TA_id)

TA_HDR_summary_all <- TAs %>%
  dplyr::left_join(TA_HDR_summary, by = c("TA_id", "Chr", "Start", "End")) %>%
  dplyr::mutate(
    n_HDRs = replace_na(n_HDRs, 0L),
    HDRs   = replace_na(HDRs, "")
  )
print(TA_HDR_summary_all)

TA_HDR_export <- TA_HDR_summary_all %>%
  dplyr::select(-TA_id, -HDRs)

write_csv(TA_HDR_export,
          "../../tables/TableS13_TA_HDRs.csv")

HDRs_hit_coords <- TA_HDR_hits_long %>%
  dplyr::distinct(
    CHROM,
    start,
    end
  ) %>%
  dplyr::rename(
    chr = CHROM
  ) %>%
  dplyr::arrange(chr, start, end)

HDRs_hit_coords_with_TA <- TA_HDR_hits_long %>%
  dplyr::distinct(TA_id, CHROM, start, end) %>%
  dplyr::group_by(CHROM, start, end) %>%
  dplyr::summarise(
    hit_by_TAs = paste(sort(unique(TA_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::rename(chr = CHROM) %>%
  dplyr::arrange(chr, start, end)

plot_df <- HDRs_hit_coords_with_TA %>%
  tidyr::separate_rows(hit_by_TAs, sep = ";") %>%
  dplyr::rename(TA_id = hit_by_TAs) %>%
  dplyr::left_join(TAs, by = "TA_id") %>%
  dplyr::mutate(
    TA_chr = str_replace(Chr, "^Chr", ""),
    chr    = str_replace(chr, "^Chr", "")
  ) %>%
  dplyr::filter(chr == TA_chr) %>%
  dplyr::group_by(TA_id) %>%
  dplyr::mutate(chr = factor(chr, levels = unique(chr))) %>%
  dplyr::ungroup()

TA_df <- TAs %>%
  dplyr::mutate(
    TA_chr = str_replace(Chr, "^Chr", ""),
    TA_id  = factor(TA_id, levels = unique(plot_df$TA_id))
  )

TAs_plot <- ggplot(plot_df, aes(x = start, xend = end, y = chr, yend = chr)) +
  geom_segment(linewidth = 2, alpha = 0.85) +
  geom_segment(
    data = TA_df,
    aes(x = Start, xend = End, y = TA_chr, yend = TA_chr),
    inherit.aes = FALSE,
    linewidth = 6,
    alpha = 0.25
  ) +
  facet_wrap(
    ~ TA_id,
    scales = "free",
    labeller = as_labeller(function(x) {
      ta <- TA_df[TA_df$TA_id == x, ]
      paste0(ta$Chr, ": ", ta$Start, "-", ta$End)
    })
  ) +
  scale_x_continuous(
    labels = function(x) x / 1e6,
    name   = "Genomic position (Mb)"
  ) +
  scale_y_discrete(drop = TRUE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(
  TAs_plot,
  filename = "../../figures/FigureS36_TAs_overlap_nonoverlap_HDRs.pdf",
  width  = 7.5,
  height = 3,
  units  = "in"
)

