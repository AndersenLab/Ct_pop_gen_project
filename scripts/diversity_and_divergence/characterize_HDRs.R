library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(ape) 
library(circlize)


lineage_colors <- c(
  LAC = "#FFCA10",
  Tw1 = "#161308",
  Tw2 = "#FF66CC",
  Indo1 = "#3390FF",
  Tw3 = "#D1007A",
  Indo2 = "#4EA3C8",
  Tw4 = "#AA3399",
  Mic1 = "#7A541F",
  Mic2 = "#A99020",
  HC = "#D9EF8B",
  Hw1 = "#1B4D1B",
  Tw5 = "#6600CC",
  Au = "#002FA7",
  Tw6 = "#A366FF",
  Hw2 = "#4C8C4C",
  Af = "#F4B6B6",
  Indo3 = "#C8E6F0",
  Tw7 = "#D9B3FF",
  Hw3 = "#00FF00"
)


hdrs <- readr::read_tsv("../../tables/TableS6_HDR_CT_allStrain_5kbclust_20251201.tsv")  %>% 
  dplyr::filter(divSize>=5e3)
bins <- readr::read_tsv("../../processed_data/genomic_bins/ONT_NIC58_1kb_bins.bed",col_names = c("CHROM","binStart","binEnd")) 
all_variants <- readr::read_tsv("../../processed_data/variant_counts/CT_all_strain_vc.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN")) %>% dplyr::mutate(source="QX1410") %>% dplyr::filter(!STRAIN=="QX1410")
lineages <- readr::read_csv("../../processed_data/geo_info/geo_and_lineage.csv") 

dt_variants <- as.data.table(all_variants)
dt_hdrs <- as.data.table(hdrs)

setnames(dt_variants, old = c("START_BIN", "END_BIN"), new = c("start", "end"))
setnames(dt_hdrs, old = c("minStart", "maxEnd"), new = c("start", "end"))

setkey(dt_variants, STRAIN, CHROM, start, end)
setkey(dt_hdrs, STRAIN, CHROM,  start, end)

overlap_result <- foverlaps(dt_variants, dt_hdrs, type = "any", nomatch = 0L)

dt_variants[, in_hd := FALSE]
overlap_bins <- overlap_result[, .(STRAIN, CHROM, start = i.start, end = i.end)]
dt_variants[overlap_bins, in_hd := TRUE, on = .(STRAIN, CHROM, start, end)]

countSummary_nr <- as.data.frame(dt_variants) %>%
  dplyr::group_by(STRAIN,in_hd) %>%
  dplyr::summarise(tot_var=sum(COUNT))

variantSummary_nr <- countSummary_nr %>%
  tidyr::pivot_wider(names_from = in_hd, values_from = tot_var, names_prefix = "tot_var_") %>%
  dplyr::rename(non_divergent_variants=tot_var_FALSE,divergent_variants=tot_var_TRUE) %>%
  dplyr::mutate(genome_wide_variants=non_divergent_variants+divergent_variants,divergent_variants_prop=divergent_variants/genome_wide_variants) 

covSpans_nr <- as.data.frame(dt_variants) %>%
  dplyr::group_by(STRAIN,in_hd) %>%
  dplyr::summarise(span=n()*1e3) %>% 
  dplyr::ungroup()

spanSummary_nr <- covSpans_nr %>%
  tidyr::pivot_wider(names_from = in_hd, values_from = span, names_prefix = "tot_bases_") %>%
  dplyr::rename(non_divergent_span=tot_bases_FALSE,divergent_span=tot_bases_TRUE) %>%
  dplyr::mutate(genome_span=non_divergent_span+divergent_span,divergent_span_prop=divergent_span/genome_span) 

allSummary_nr <- variantSummary_nr %>%
  dplyr::left_join(spanSummary_nr,by="STRAIN") %>%
  dplyr::left_join(lineages,by=c("STRAIN"="strain")) %>%
  dplyr::mutate(lineage=factor(lineage,levels=sort(unique(lineages$lineage))))

meanRGSummary_nr <- allSummary_nr %>%
  dplyr::mutate(source="NIC58") %>%
  dplyr::group_by(lineage) %>%
  dplyr::summarise(mean_var_prop = mean(divergent_variants_prop),
                   mean_var_prop_sd=sd(divergent_variants_prop),
                   min_var_prop=min(divergent_variants_prop),
                   max_var_prop=max(divergent_variants_prop),
                   mean_span_prop=mean(divergent_span_prop),
                   mean_span_prop_sd=sd(divergent_span_prop),
                   min_span_prop=min(divergent_span_prop),
                   max_span_prop=max(divergent_span_prop)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lineage=factor(lineage,levels=sort(unique(lineages$lineage))))

meanRGSummary_tot <- allSummary_nr %>%
  dplyr::mutate(source="NIC58") %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(mean_var_prop = mean(divergent_variants_prop),
                   mean_var_prop_sd=sd(divergent_variants_prop),
                   min_var_prop=min(divergent_variants_prop),
                   max_var_prop=max(divergent_variants_prop),
                   mean_span_prop=mean(divergent_span_prop),
                   mean_span_prop_sd=sd(divergent_span_prop),
                   min_span_prop=min(divergent_span_prop),
                   max_span_prop=max(divergent_span_prop)) %>%
  dplyr::ungroup() 

meanRGSummary_nr <- rbind(meanRGSummary_nr,meanRGSummary_tot %>% dplyr::rename(lineage=source) %>% dplyr::mutate(lineage=ifelse(lineage=="NIC58","All",lineage)))

p1_nr <- ggplot() +
  geom_point(data=allSummary_nr, aes(x = divergent_span_prop * 100, y = divergent_variants_prop * 100, color = lineage),shape=21,stroke=0.8) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
    # axis.text=element_blank(),
    # axis.ticks=element_blank(),
    plot.margin = margin(5.5, 1, 1, 1)
  ) +
  scale_color_manual(values = lineage_colors,
                     breaks = names(lineage_colors)) +
  scale_y_continuous(breaks = seq(0, 90, 10),limits = c(0,90)) +
  scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,24)) +
  labs(color = "Relatedness group") +
  xlab("Percent genome span of hyper-divergent regions") +
  ylab("Percent of variants in hyper-divergent regions")


p2_nr <- ggplot() +
  geom_point(data=meanRGSummary_nr,
             aes(x = mean_span_prop * 100, y = mean_var_prop * 100, color = lineage)) +
  geom_errorbar(data=meanRGSummary_nr,aes(x = mean_span_prop * 100,
                   ymin = (mean_var_prop - mean_var_prop_sd) * 100,
                   ymax = (mean_var_prop + mean_var_prop_sd) * 100,
                   color = lineage),
               width = 0) + 
  geom_errorbarh(data=meanRGSummary_nr,aes(y = mean_var_prop * 100,
                     xmin = (mean_span_prop - mean_span_prop_sd) * 100,
                     xmax = (mean_span_prop + mean_span_prop_sd) * 100,
                     color = lineage),
                 height = 0) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin = margin(5.5, 1, 1, 1)) +
  scale_color_manual(values = lineage_colors,
                     breaks = names(lineage_colors)) +
  scale_y_continuous(breaks = seq(0, 90, 10),limits = c(60,90)) +
  scale_x_continuous(breaks = seq(0, 100, 2),limits = c(0,24)) +
  labs(color = "Relatedness\ngroup") +
  xlab("Percent genome span of hyper-divergent regions") +
  ylab("Percent of variants in hyper-divergent regions")


legend <- cowplot::get_legend(p2_nr)

p12_nr <- cowplot::plot_grid(p2_nr + theme(legend.position="none"),p1_nr,rel_heights = c(0.4,1),nrow=2)
padded_plot <- cowplot::plot_grid(NULL,p12_nr,legend,rel_widths = c(0.035,1,0.2),nrow=1)

full_plot <- cowplot::ggdraw(padded_plot) +
  #draw_label("Percent genome span of hyper-divergent regions", x = 0.5, y = 0.01, vjust = 0, size = 12) +
  draw_label("Percent of variants in hyper-divergent regions", x = 0.005, y = 0.5, angle = 90, vjust = 1, size = 12)
#full_plot

ggsave(plot = full_plot, filename = "../../figures/FigureS28_propVC_20260323.png",width = 7,height = 6,dpi = 600,device = 'png',bg = "white")

write.table(meanRGSummary_nr,file = "../../tables/TableS7_HDRsummaryStats_20251205.tsv",sep = "\t",quote = F,row.names = F)

##########
conc <- readr::read_tsv("../../data/gtcheck.txt") %>%
  dplyr::filter(i=="NIC58" | j =="NIC58") %>%
  dplyr::mutate(isotype=ifelse(i=="NIC58",j,i)) %>%
  dplyr::mutate(conc=(sites-discordance)/sites) %>%
  dplyr::select(isotype,conc,discordance)

allSummary_nr_wConc <- allSummary_nr %>%
  dplyr::left_join(conc,by="isotype")

dat <- allSummary_nr_wConc %>%
  dplyr::filter(divergent_variants > 20)

t1 <- ggplot() +
  geom_point(data = dat,
             aes(x = divergent_variants_prop * 100,
                 y = conc*100,
                 color = lineage),
             shape = 21, stroke = 0.8) +
  theme_bw() +
  theme(
    panel.grid     = element_blank(),
    axis.title.y   = element_blank(),
    legend.position = 'none',
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    plot.margin    = margin(5.5, 1, 2, 1)
  ) +
  scale_color_manual(values = lineage_colors,
                     breaks = names(lineage_colors)) +
  labs(color = "Relatedness group") +
  ylab("Percent genetic similarity with NIC58") +
  xlab("Percent of variants in HDRs")

t2 <- ggplot() +
  geom_point(data = dat,
             aes(x = divergent_span_prop * 100,
                 y = conc*100,
                 color = lineage),
             shape = 21, stroke = 0.8) +
  theme_bw() +
  theme(
    panel.grid   = element_blank(),
    legend.position = 'none',
    plot.margin  = margin(5.5, 1, 2, 2)
  ) +
  scale_color_manual(values = lineage_colors,
                     breaks = names(lineage_colors)) +
  labs(color = "Relatedness group") +
  xlab("Percent genome span of HDRs") +
  ylab("Percent genetic similarity with NIC58")

t2_bar <- ggplot(dat,
                 aes(x = divergent_span_prop * 100, fill = lineage)) +
  geom_histogram(binwidth = 1, position = "stack", color = NA) +
  scale_fill_manual(values = lineage_colors,
                    breaks = names(lineage_colors)) +
  theme_bw() +
  theme(
    panel.grid    = element_blank(),
    axis.title.x    = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    legend.position = "none",
    plot.margin   = margin(5.5, 1, 0, 2)
  ) +
  ylab("Number of\nisotypes")

t1_bar <- ggplot(dat,
                 aes(x = divergent_variants_prop * 100, fill = lineage)) +
  geom_histogram(binwidth = 1, position = "stack", color = NA) +
  scale_fill_manual(values = lineage_colors,
                    breaks = names(lineage_colors)) +
  theme_bw() +
  theme(
    panel.grid    = element_blank(),
    axis.title    = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    legend.position = "none",
    plot.margin   = margin(5.5, 1, 0, 1)
  ) +
  ylab("Number of\nisotypes")

p2_combined <- plot_grid(t2_bar, t2, ncol = 1, rel_heights = c(0.3, 1),align = "v",axis="lr")
p1_combined <- plot_grid(t1_bar, t1, ncol = 1, rel_heights = c(0.3, 1),align = "v",axis="lr")
#p3_combined <- plot_grid(t3_bar, t3, ncol = 1, rel_heights = c(0.3, 1),align = "v",axis="lr")

simi_metrics <- cowplot::plot_grid(p2_combined, p1_combined,legend,
                   nrow = 1,
                   rel_widths = c(1, 0.9, 0.4))
ggsave(plot = simi_metrics, filename = "../../figures/FigureS29_similarity_hdrs_20260323.png",width = 7,height = 6,dpi = 600,device = 'png',bg = "white")

fold_average <- allSummary_nr_wConc %>%
  dplyr::mutate(divergent_variant_density=divergent_variants/divergent_span,genome_wide_variant_density=genome_wide_variants/genome_span) %>%
  dplyr::mutate(fold_change=divergent_variant_density/genome_wide_variant_density) 

size_summary <- hdrs %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarize(total_divSize = sum(divSize, na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(total_divSize)) %>%
  dplyr::mutate(prop=total_divSize/(nrow(bins)*1e3))

meanSize <- mean(hdrs$divSize)/1e3
minSize <- min(hdrs$divSize)/1e3
maxSize <- max(hdrs$divSize)/1e3

getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      #print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>% 
        dplyr::arrange(CHROM,minStart) %>%
        dplyr::mutate(check=ifelse(lead(minStart) <= maxEnd,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
      if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
        print("NO MORE INTERSECTS")
        k=0
      } else {
        
        temp <- checkIntersect %>%
          dplyr::mutate(gid=data.table::rleid(check)) %>%
          dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))
        
        collapse <- temp %>%
          dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart=min(minStart)) %>%
          dplyr::mutate(newEnd=max(maxEnd)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(minStart=newStart,maxEnd=newEnd) %>%
          dplyr::select(-newEnd,-newStart)
        
        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print("WRITING TO LIST")
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

collapsed <- plyr::ldply(getRegFreq(hdrs %>% dplyr::group_split(CHROM)),data.frame) %>% dplyr::mutate(divSize=maxEnd-minStart)
cov_wg <- sum(collapsed$divSize) / (nrow(bins)*1e3)
wg_tot <- nrow(collapsed)
wg_mean <- mean(collapsed$divSize) / 1e3
wg_min <- min(collapsed$divSize) / 1e3
wg_max <- max(collapsed$divSize) / 1e3

bins_dt <- as.data.table(bins)
setnames(bins_dt, c("binStart", "binEnd"), c("start", "end"))
bins_dt[, id := .I]  # optional: keep track of bins

hdrs_dt <- as.data.table(hdrs)
setnames(hdrs_dt, c("minStart", "maxEnd"), c("start", "end"))

setkey(bins_dt, CHROM, start, end)
setkey(hdrs_dt, CHROM, start, end)

overlaps <- foverlaps(hdrs_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(STRAIN)), by = .(CHROM, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("CHROM", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq=n_strains/695)

hdrs_ordered <- hdrs %>% 
  dplyr::left_join(lineages %>% dplyr::select(strain,lineage),by=c("STRAIN"="strain")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() 


x_breaks_by_chrom <- bins_wFreq %>%
  group_by(CHROM) %>%
  summarise(x_breaks = list(seq(
    floor(min((start) / 1e6)),
    ceiling(max((end) / 1e6)),
    by = 5
  )))

p1 <- ggplot(hdrs_ordered) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45)) + 
  facet_wrap(~CHROM,scales = 'free_x',nrow=1) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text = element_blank())  +
  ylab("622 Isotype strains") +
  # scale_fill_manual(values = lineage_colors,
  #                   breaks = names(lineage_colors)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 5),expand = c(0, 0))

p2 <- ggplot(bins_wFreq %>% dplyr::filter(CHROM!="MtDNA")) +
  geom_point(aes(x=(start+500)/1e6,y=freq),size=0.3,stroke = 0) +
  facet_wrap(~CHROM,scales = 'free_x',nrow=1) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        strip.text = element_blank())  +
  ylab("Frequency")+
  xlab("Physical position (Mb)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x+(500/1e6))), by = 5),expand = c(0, 0))

#p21 <- cowplot::plot_grid(p1,p2,nrow=2,ncol=1,rel_heights = c(1,0.3),align = "v", axis = "l", labels = c("a","b"))

diversity <- readr::read_csv("../../processed_data/pi_theta_d/chromosome_windows_diversity.csv",col_select = c(-1))

domains_raw <- readr::read_tsv("../../data/chromosome_windows_diversity.csv") 

region_rects2 <- domains_raw %>% 
  dplyr::mutate(region=ifelse(grepl("tip",Location),"Tip",ifelse(grepl("arm",Location),"Arm","Center"))) %>%
  dplyr::mutate(ymin=-Inf,ymax=Inf,xmin=start/1e6,xmax=end/1e6)
region_colors <- c("Tip" = "#5E3C99", "Center" = "#FDB863", "Arm" = "#4393C3")

pitheta_plot <- ggplot() + 
  geom_rect(data = region_rects2, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(data=diversity %>% 
               dplyr::filter(stat_type!="d" & stat_type!="mis") %>%
               dplyr::rename(CHROM=chrom) %>%
               dplyr::mutate(stat_type = factor(ifelse(stat_type == "pi", "π", "Watterson's θ"), 
                                                levels = c("π", "Watterson's θ"))),
             aes(x=x/1e6,y=stat),
             size=0.05) +
  geom_smooth(data=diversity %>% 
                dplyr::filter(stat_type!="d" & stat_type!="mis") %>%
                dplyr::rename(CHROM=chrom) %>%
                dplyr::mutate(stat_type = factor(ifelse(stat_type == "pi", "π", "Watterson's θ"), 
                                                 levels = c("π", "Watterson's θ"))),
              aes(x=x/1e6,y=stat),
              method="loess", se=FALSE, color="grey85", linewidth=1,span=0.3) + 
  scale_fill_manual(values = region_colors) +
  facet_grid(stat_type~CHROM,scales = "free_x") +
  theme_bw() +
  theme(panel.border = element_rect(fill=NA),
        strip.background.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        plot.margin = margin(0, 0, 5, 10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=9),
        strip.background = element_blank(),
        strip.text.y=element_text(size=10)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill="Domain")

div_plot <- cowplot::plot_grid(pitheta_plot,p1,p2,nrow=3,align = "v",axis = "lr",rel_heights = c(0.5,0.9,0.3),labels = c("a","b"))
ggsave(div_plot,filename = "../../figures/Figure3_div_hd_main.png",width = 7,height = 8,dpi = 600,device = 'png')

hdrs_ordered_priv_iv <- hdrs %>% 
  dplyr::left_join(lineages %>% dplyr::select(strain,lineage),by=c("STRAIN"="strain")) %>%
  dplyr::filter(lineage %in% c("Hw3","LAC"),CHROM=="IV") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ncalls,STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup()

# identify 20 evenly spaced rleID values
rle_ids_to_keep <- hdrs_ordered_priv_iv %>%
  filter(lineage == "LAC") %>%
  distinct(rleID) %>%
  mutate(idx = row_number()) %>%
  filter(idx %in% round(seq(1, n(), length.out = 20))) %>%
  pull(rleID)

# now keep ALL rows with those rleID values
lac_subsample_20 <- hdrs_ordered_priv_iv %>%
  dplyr::filter(lineage == "LAC", rleID %in% rle_ids_to_keep) %>%
  dplyr::group_by(rleID) %>%
  dplyr::mutate(rleID_new = cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::select(-rleID) %>%
  dplyr::rename(rleID=rleID_new)# re-enumerate rleID starting at 1

lac20_hw3 <- rbind(lac_subsample_20,hdrs_ordered_priv_iv %>%
                     dplyr::filter(lineage == "Hw3"))


strain_labs <- lac20_hw3 %>%
  distinct(rleID, STRAIN) %>%
  arrange(rleID) %>%
  tibble::deframe()

NICvECA <- readr::read_tsv("../../processed_data/genome_alignments/CT_all_nucmer_aln.tsv",
                           col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")) %>%
  dplyr::filter(STRAIN=="ECA1307")

p_priv <- ggplot(lac20_hw3) + 
  geom_rect(aes(xmin=minStart/1e6,
                xmax=maxEnd/1e6,
                ymin=rleID-0.45,
                ymax=rleID+0.45)) + 
  facet_grid(lineage~CHROM, scales='free_y', space="free_y") + 
  scale_y_continuous(breaks = unique(lac20_hw3$rleID),
                     labels = strain_labs,
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 5),
                     expand = c(0, 0),
                     limits = c(0, max(lac20_hw3$maxEnd)/1e6))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA)) +
  ylab("") +
  xlab("NIC58 phyisical position (Mb)")

priv_aln <- ggplot(NICvECA %>% dplyr::filter(REF == "IV")) +
  geom_rect(aes(xmin=S1/1e6, xmax=E1/1e6, ymin=IDY+0.1, ymax=IDY-0.1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic"),
    strip.background = element_blank()
  ) +
  scale_x_continuous(expand = c(0,0),
                     limits= c(0,max(lac20_hw3$maxEnd)/1e6)) +
  ylab("Percent identity") +
  xlab("") +
  annotate(
    "text",
    x = -Inf, y = -Inf,                    # bottom-left corner
    label = "NIC58 (LAC) vs. ECA1307 (Hw3) - Chromosome IV",
    hjust = -0.02, vjust = -1,            # nudges text inward
    size = 3
  ) +
  xlab("NIC58 phyisical position (Mb)") 

p_range <- ggplot_build(priv_aln)$layout$panel_params[[1]]
x_min <- p_range$x_range$range[1]
x_max <- p_range$x_range$range[2]

fix_diff <- cowplot::plot_grid(p_priv,priv_aln,nrow=2,align = "v",axis = "lr",labels = c("a","b"),rel_heights = c(1,0.6))
ggsave(fix_diff,filename = "../../figures/FigureS30_fixed_differences.png",width = 7,height = 6.5,dpi = 600,device = 'png')
