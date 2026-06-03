rm(list = ls())

library(dplyr)
library(readr)
library(plyr)
library(data.table)

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
          # dplyr::filter(check == FALSE & lag(check) == FALSE)
        
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
  "../../tables/TableS6_HDR_CT_allStrain_5kbclust_20251201.tsv",
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

# ============================================================
# If you want genome percentage, use the same bins file
# ============================================================

bins <- readr::read_tsv(
  "../../processed_data/genomic_bins/ONT_NIC58_1kb_bins.bed",
  col_names = c("CHROM", "binStart", "binEnd"),
  show_col_types = FALSE
) %>%
  dplyr::filter(CHROM != "MtDNA")

genome_span_bp <- nrow(bins) * 1e3

summary_nonoverlap_hdrs <- summary_nonoverlap_hdrs %>%
  dplyr::mutate(
    genome_span_bp = genome_span_bp,
    percent_reference_genome = total_span_bp / genome_span_bp * 100
  )

print(summary_nonoverlap_hdrs)


# ============================================================
# Additional HDR summaries 
# Questions:
# 1. How many genes does each HDR contain?
# 2. What is the average and range of gene number per HDR?
# 3. How many HDR SNVs are in exons, introns, and intergenic regions?
# 4. How many HDR SNVs are in chromosome domains?
# ============================================================

###### Add HDR ID ######
nonoverlap_hdrs <- nonoverlap_hdrs %>%
  dplyr::arrange(CHROM, start, end) %>%
  dplyr::mutate(HDR_id = paste0("HDR_", dplyr::row_number())) %>%
  dplyr::select(HDR_id, CHROM, start, end, divSize)

head(nonoverlap_hdrs)
dim(nonoverlap_hdrs)


library(rtracklayer)
gff_gr <- rtracklayer::import("../../processed_data/gene_enrichment/c_tropicalis.NIC58_20251002.csq.longest.gff3")
gff_raw <- as.data.frame(gff_gr)

head(gff_raw)
colnames(gff_raw)
table(gff_raw$type)

gene_regions <- gff_raw %>%
  dplyr::filter(type == "gene") %>%
  dplyr::filter(seqnames != "MtDNA") %>%
  dplyr::mutate(
    CHROM = as.character(seqnames),
    gene_id = ID,
    gene_start = start - 1,
    gene_end = end
  ) %>%
  dplyr::select(CHROM, gene_start, gene_end, gene_id) %>%
  dplyr::arrange(CHROM, gene_start, gene_end)

head(gene_regions)
dim(gene_regions)


###### Exon intervals from NIC58 GFF ######

exon_regions <- gff_raw %>%
  dplyr::filter(type == "exon") %>%
  dplyr::filter(seqnames != "MtDNA") %>%
  dplyr::mutate(
    CHROM = as.character(seqnames),
    exon_start = start - 1,
    exon_end = end
  ) %>%
  dplyr::select(CHROM, exon_start, exon_end) %>%
  dplyr::distinct() %>%
  dplyr::arrange(CHROM, exon_start, exon_end)

head(exon_regions)
dim(exon_regions)

# ============================================================
# 1. Count genes in each non-overlapping HDR
# ============================================================

hdr_for_overlap <- nonoverlap_hdrs %>%
  dplyr::select(HDR_id, CHROM, start, end, divSize)

gene_for_overlap <- gene_regions %>%
  dplyr::rename(start = gene_start,
                end = gene_end)

hdr_dt <- as.data.table(hdr_for_overlap)
gene_dt <- as.data.table(gene_for_overlap)

setkey(hdr_dt, CHROM, start, end)
setkey(gene_dt, CHROM, start, end)

###### Overlap genes with HDRs ######
hdr_gene_overlap <- data.table::foverlaps(gene_dt,
                                          hdr_dt,
                                          by.x = c("CHROM","start","end"),
                                          by.y = c("CHROM","start","end"),
                                          type = "any",
                                          nomatch = 0)

head(hdr_gene_overlap)
dim(hdr_gene_overlap)


###### Gene count per HDR ######
gene_count_per_HDR <- hdr_gene_overlap %>%
  as.data.frame() %>%
  dplyr::select(HDR_id, gene_id) %>%
  dplyr::distinct() %>%
  dplyr::group_by(HDR_id) %>%
  dplyr::summarise(n_genes = dplyr::n()) %>%
  dplyr::ungroup()

gene_count_per_HDR <- nonoverlap_hdrs %>%
  dplyr::select(HDR_id, CHROM, start, end, divSize) %>%
  dplyr::left_join(gene_count_per_HDR, by = "HDR_id") %>%
  dplyr::mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes))

head(gene_count_per_HDR)
summary(gene_count_per_HDR$n_genes)


###### Summary of gene number per HDR ######
summary_gene_count_per_HDR <- gene_count_per_HDR %>%
  dplyr::summarise(
    n_HDRs = dplyr::n(),
    total_gene_overlaps = sum(n_genes),
    mean_genes_per_HDR = mean(n_genes),
    median_genes_per_HDR = median(n_genes),
    min_genes_per_HDR = min(n_genes),
    max_genes_per_HDR = max(n_genes)
  )

print(summary_gene_count_per_HDR)


# ============================================================
# 2. Classify HDR SNVs into exon, intron, and intergenic
# ============================================================

###### Load SNV file with HDR/non-HDR labels ######

snv_data <- read.table("../../processed_data/HDR_stats/Ct_soft_filtered_hets.het_counts_per_site.HDR_labeled.tsv",
                       header = TRUE)

head(snv_data)
dim(snv_data)
table(snv_data$region)

hdr_snv_data <- snv_data %>%
  dplyr::filter(region == "HDR") %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(
    snv_id = paste(CHROM, POS, sep = ":"),
    snv_start = POS - 1,
    snv_end = POS
  ) %>%
  dplyr::select(snv_id, CHROM, POS, snv_start, snv_end,
                n_het, n_called, n_total, has_het, region)

head(hdr_snv_data)
dim(hdr_snv_data)

snv_for_overlap <- hdr_snv_data %>%
  dplyr::rename(start = snv_start,
                end = snv_end)

snv_dt <- as.data.table(snv_for_overlap)

setkey(snv_dt, CHROM, start, end)

exon_for_overlap <- exon_regions %>%
  dplyr::rename(start = exon_start,
                end = exon_end)

exon_dt <- as.data.table(exon_for_overlap)

setkey(exon_dt, CHROM, start, end)

snv_exon_overlap <- data.table::foverlaps(snv_dt,
                                          exon_dt,
                                          by.x = c("CHROM","start","end"),
                                          by.y = c("CHROM","start","end"),
                                          type = "any",
                                          nomatch = 0)

exon_snv_ids <- unique(snv_exon_overlap$snv_id)
length(exon_snv_ids)


gene_snv_for_overlap <- gene_regions %>%
  dplyr::select(CHROM, gene_start, gene_end) %>%
  dplyr::distinct() %>%
  dplyr::rename(start = gene_start,
                end = gene_end)

gene_snv_dt <- as.data.table(gene_snv_for_overlap)

setkey(gene_snv_dt, CHROM, start, end)

snv_gene_overlap <- data.table::foverlaps(snv_dt,
                                          gene_snv_dt,
                                          by.x = c("CHROM","start","end"),
                                          by.y = c("CHROM","start","end"),
                                          type = "any",
                                          nomatch = 0)

gene_snv_ids <- unique(snv_gene_overlap$snv_id)
length(gene_snv_ids)

hdr_snv_classified <- hdr_snv_data %>%
  dplyr::mutate(
    feature_class = dplyr::case_when(
      snv_id %in% exon_snv_ids ~ "exon",
      snv_id %in% gene_snv_ids ~ "intron",
      TRUE ~ "intergenic"
    )
  )

table(hdr_snv_classified$feature_class)

summary_HDR_SNV_feature_class <- hdr_snv_classified %>%
  dplyr::group_by(feature_class) %>%
  dplyr::summarise(
    n_SNVs = dplyr::n()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    percent_SNVs = n_SNVs / sum(n_SNVs) * 100
  )

print(summary_HDR_SNV_feature_class)


###### Summary: genic / intergenic ######
summary_HDR_SNV_genic_intergenic <- hdr_snv_classified %>%
  dplyr::mutate(
    genic_status = ifelse(feature_class %in% c("exon","intron"),
                          "genic",
                          "intergenic")
  ) %>%
  dplyr::group_by(genic_status) %>%
  dplyr::summarise(
    n_SNVs = dplyr::n()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    percent_SNVs = n_SNVs / sum(n_SNVs) * 100
  )

print(summary_HDR_SNV_genic_intergenic)
