library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ape)
library(data.table)
library(valr)
library(stringr)

# This function intersects genome-to-genome (g2g) long-read alignments with genome bins to extract bin-level coverage and identity 
# a miriad of bin metrics are calculated (see mutate cluster)
alignmentPartition <- function(df) {
  a2a_bed <- df %>% dplyr::select(REF,S1,E1,IDY,STRAIN) %>%
    dplyr::rename(chrom=REF,start=S1,end=E1,Identity=IDY)
  
  df_LRbins <- valr::bed_intersect(bins_1kb_CE_stripped,a2a_bed) %>%
    dplyr::rename(CHROM=chrom, START_BIN=start.x, END_BIN=end.x, 
                  start=start.y, stop=end.y, coverage=.overlap, Identity=Identity.y)
  
  df_LRbins_unclass <- bins_1kb_CE_stripped %>% 
    dplyr::rename(CHROM=chrom, START_BIN=start, END_BIN=end) %>%
    dplyr::left_join(df_LRbins,by=c("CHROM"="CHROM","START_BIN"="START_BIN","END_BIN"="END_BIN")) %>%
    dplyr::arrange(CHROM,START_BIN) %>% 
    dplyr::mutate(coverage=as.numeric(coverage)) %>%
    dplyr::mutate(coverage=ifelse(is.na(coverage),0,coverage)) %>%
    dplyr::mutate(STRAIN.y=ifelse(is.na(STRAIN.y),unique(a2a_bed$STRAIN),STRAIN.y)) %>%
    dplyr::mutate(full_cov = ifelse(coverage == 1e3, T, F)) %>% 
    dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
    dplyr::mutate(group_size=n()) %>%
    dplyr::mutate(group_id=cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(group_id) %>%
    dplyr::mutate(maxIDY=max(Identity)) %>%
    #dplyr::mutate(minIDY=min(Identity)) %>%
    #dplyr::mutate(minCov=min(coverage)) %>%
    dplyr::mutate(maxCov=max(coverage)) %>%
    dplyr::mutate(bin_IDY=maxIDY) %>%
    dplyr::mutate(bin_COV=maxCov) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(CHROM,START_BIN) %>% 
    dplyr::select(CHROM,START_BIN,END_BIN,Identity,coverage,bin_IDY,bin_COV,group_id,group_size,STRAIN.y) %>%
    dplyr::rename(STRAIN=STRAIN.y) %>%
    dplyr::distinct(group_id,.keep_all = T) 
  
  return(df_LRbins_unclass)
}

# given a set of coverage and identity thresholds, this function classifies bins as divergent or non-divergent
# divergent bins have sub-classifications based on the variant determinants 
classifyPartition <- function(df,cf,idy) {
  df_bins_clasif_LR <- df %>%
    dplyr::mutate(div=ifelse(coverage<(cf*10),"C",ifelse(!(is.na(Identity)) & Identity<=idy,"I","nondiv"))) %>%
    dplyr::arrange(CHROM,START_BIN) %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(div_class = ifelse((lag(div) == "C" | lag(div) == "I") & (lead(div) == "C" | lead(div) == "I") & div=="nondiv","G", div)) %>%
    dplyr::mutate(div_class = ifelse(div_class=="C" & coverage == 0,"Z",div_class)) %>%
    dplyr::mutate(div_gf=ifelse(div_class=="nondiv","nondiv","div")) %>%
    dplyr::ungroup()
  
  return(df_bins_clasif_LR)
}

## this function estimates the frequency at which each bin is classified as hyperdivergent in the sampled population
## low frequency () bins are discarded to avoid super-clustering of divergent regions
frequencyEstimates <- function(df) {
  denom<-length(unique(df$STRAIN))
  df <- all_stats %>% 
    dplyr::group_by(CHROM,START_BIN) %>%
    dplyr::mutate(binCt=stringr::str_count(div_gf,"nondiv")) %>%
    dplyr::mutate(binCt=ifelse(is.na(binCt),1,binCt)) %>%
    dplyr::mutate(binCt=ifelse(binCt==1,0,1)) %>%
    dplyr::mutate(binCtSum=sum(binCt)) %>%
    dplyr::mutate(binFreq=(binCtSum/denom)*100)
  return(df)
}

## This function clusters contiguous div bins into div regions, 
## and assembles a divergent 'footprint' from the string of bin classifications
clusterBins <- function(df,mode) {
  
  if(mode=="SRF") {
    temp <- df %>% 
      dplyr::arrange(CHROM,START_BIN) %>%
      dplyr::mutate(div_gf=ifelse(div_gf =="div" & binFreq < 5,"nondiv",div_gf)) %>%
      dplyr::select(-binCt,-binCtSum,-binFreq)
  } else {
    temp <- df %>% 
      dplyr::arrange(CHROM,START_BIN) 
  }
  
  temp$enum <- sequence(rle(as.character(temp$div_gf))$lengths)
  
  div_bins <- temp %>%
    dplyr::arrange(CHROM,START_BIN) %>%
    dplyr::group_by(CHROM,data.table::rleid(div_gf)) %>%
    dplyr::mutate(gid=cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(clusTips=ifelse(div_gf == "div" & lead(enum) == 2 & lead(div_gf)=="div","clust_start", 
                                  ifelse(div_gf == "div" & lead(enum) > 1 & lead(div_gf)=="div","clust_center",
                                         ifelse(div_gf == "div" & enum > 1 & lead(div_gf)=="nondiv","clust_end","unclust")))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!(clusTips=="unclust"))
  
  if (mode == "LR") {
    div_regions <- div_bins %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(bin_covz=sum(stringr::str_detect(div_class,"Z"))/n()) %>%
      dplyr::mutate(bin_lowcov=sum(stringr::str_detect(div_class,"C"))/n()) %>%
      dplyr::mutate(bin_gf=sum(stringr::str_detect(div_class,"G"))/n()) %>%
      dplyr::mutate(bin_idy=sum(stringr::str_detect(div_class,"I"))/n()) %>%
      #dplyr::mutate(reg_idy=mean(Identity)) %>%
      dplyr::mutate(meanIDY=mean(bin_IDY,na.rm=TRUE),meanCov=mean(bin_COV),minStart=min(START_BIN),maxEnd=max(END_BIN),divSize=maxEnd-minStart) %>%
      dplyr::mutate(bin_foot=paste0(div_class,collapse = "")) %>%
      dplyr::distinct(gid,.keep_all = T) %>%
      dplyr::ungroup() %>%
      dplyr::select(CHROM,minStart,maxEnd,divSize,meanIDY,meanCov,bin_covz,bin_lowcov,bin_gf,bin_idy,bin_foot)
  } else {
    div_regions <- div_bins %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(bin_covz=sum(stringr::str_detect(div_class,"Z"))/n()) %>%
      dplyr::mutate(bin_lowcov=sum(stringr::str_detect(div_class,"C"))/n()) %>%
      dplyr::mutate(bin_gf=sum(stringr::str_detect(div_class,"G"))/n()) %>%
      dplyr::mutate(bin_idy=sum(stringr::str_detect(div_class,"I"))/n()) %>%
      #dplyr::mutate(reg_idy=mean(Identity)) %>%
      dplyr::mutate(meanVC=mean(COUNT),meanCF=mean(pc1X),minStart=min(START_BIN),maxEnd=max(END_BIN),divSize=maxEnd-minStart) %>%
      dplyr::mutate(bin_foot=paste0(div_class,collapse = "")) %>%
      dplyr::distinct(gid,.keep_all = T) %>%
      dplyr::ungroup() %>%
      dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_covz,bin_lowcov,bin_gf,bin_idy,bin_foot)
  }
  
  out <- list(div_bins,div_regions)
  return(out)
}

#get nucmer g2g coords
transformed_coords <- readr::read_tsv("../../processed_data/genome_alignments/CT_all_nucmer_aln.tsv",col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")) 

#list of strains with long-read genomes
div_str <- unique(transformed_coords$STRAIN)

#get 1kb windows (bedtools makewindows)
bins_1kb_CE_stripped <- readr::read_tsv("../../processed_data/genomic_bins/ONT_NIC58_1kb_bins.bed",col_names = F) %>%
  dplyr::rename(chrom=X1,start=X2,end=X3)

#Filter out alignments smaller than 1kb
#tcords is now a list of datframes that contain the g2g alignment coordinates for each strain
tcords <- transformed_coords %>%
  dplyr::filter(L2>1000) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::group_split()

#use alignementPartition() to extract coverage and idy stats from g2g coordinate sets for each bin
strain_df <- list()
for (i in 1:length(tcords)) {
  print(i)
  tmp <- alignmentPartition(tcords[[i]])
  strain_df[[i]]  <- tmp
}
all_LR_stats <-ldply(strain_df,data.frame)

####### CALL LR BASED HDRs ##########
#set a range of identities to test HDR calling using LR data
cf <- 60
idthresh <- c(seq(90,99,1))

#iterate through identity range and classify partitions previously generated
strain_dv <- list()
clasiBins <- list()
for (i in 1:length(strain_df)) {
  print(i)
  strain_varID <- list()
  for (k in 1:length(idthresh)) {
    strain_varID[[k]] <- classifyPartition(strain_df[[i]],cf,idthresh[k])
    strain_varID[[k]]$threshIDY <- idthresh[k]
  }
  strain_dv[[i]] <- strain_varID
  clasiBins[[i]] <- ldply(strain_varID,data.frame)
}

# this df holds bin data after HDR classification, mainly for diagnostic purposes, can be omitted
# all_bins_LR <- ldply(clasiBins,data.frame) %>%
#   dplyr::filter(threshIDY==95 & !(CHROM=="MtDNA")) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(gwIDY=mean(bin_IDY,na.rm=T),gwIDY_SD=sd(bin_IDY,na.rm=T)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(swIDY=mean(gwIDY,na.rm=T),swIDY_SD=sd(gwIDY,na.rm=T)) %>%
#   dplyr::group_by(STRAIN,CHROM) %>%
#   dplyr::mutate(chrIDY=mean(bin_IDY,na.rm=T),chrIDY_SD=sd(bin_IDY,na.rm=T)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(CHROM) %>%
#   dplyr::mutate(chrIDY_sw=mean(chrIDY,na.rm=T),chrIDY_sw_SD=sd(chrIDY,na.rm=T)) %>%
#   dplyr::ungroup()

#iterate through bin classifications at various idy thresholds to cluster bins into regions using clusterBins()
dvReg <- list()
dvReg_part <- list()
for (i in 1:length(strain_df)) {
  print(i)
  strain_declust <- data.frame(CHROM=character(),minStart=integer(),maxEnd=integer(),divSize=integer(),meanIDY=numeric(),treshIDY=integer())
  bin_declust <- data.frame(CHROM=character(),START_BIN=integer(),END_BIN=integer(),bin_IDY=double(),bin_COV=double(),group_size=integer(),div_class=character(),gid=integer(),clusTips=character())
  strain <- unique(strain_dv[[i]][[1]]$STRAIN)
  
  for (k in 1:length(idthresh)) {
    temp <- clusterBins(strain_dv[[i]][[k]],"LR")
    temp_bins <- temp[[1]] %>%
      dplyr::select(CHROM,START_BIN,END_BIN,bin_IDY,bin_COV,group_size,div_class,gid,clusTips)
    temp_div <- temp[[2]]
    temp_div$threshIDY <- idthresh[k]
    temp_bins$threshIDY <- idthresh[k]
    strain_declust <- rbind(strain_declust,temp_div)
    bin_declust <- rbind(bin_declust,temp_bins)
  }
  strain_declust$STRAIN <- strain
  bin_declust$STRAIN <- strain
  
  dvReg[[i]] <- strain_declust %>%
    dplyr::mutate(regionClass=ifelse(bin_covz>0.95,"coverage gap",
                                     ifelse(bin_idy>0.5, "identity call", "coverage call")))
  dvReg_part[[i]] <- bin_declust
}
all_calls_LR <- ldply(dvReg, data.frame) %>% dplyr::filter(CHROM!="MtDNA")

#write.table(all_calls_LR, "processed_data/LR_calls/LR_HDRs_allThresh.tsv",quote = F,sep = "\t",row.names = F)
#write.table(all_LR_stats, "processed_data/LR_calls/LR_HDRs_allThresh_stats.tsv",quote = F,sep = "\t",row.names = F)

s1 <- ggplot(all_calls_LR %>% dplyr::filter(divSize > 200e3)) + 
  geom_histogram(aes(x=divSize/1e3),binwidth = 5) + 
  facet_grid(CHROM~threshIDY, scales = "free_x") +
  theme_classic() +
  theme(strip.text.x=element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1)) +
  xlab("HDR size (kb)") +
  scale_y_continuous(breaks = seq(0, max(all_calls_LR$divSize/1e3, na.rm = TRUE), by = 5))

s2 <- ggplot(all_calls_LR %>% dplyr::filter(divSize > 50e3 & divSize <=200e3)) + 
  geom_histogram(aes(x=divSize/1e3),binwidth = 5) + 
  facet_grid(CHROM~threshIDY, scales = "free_x")  +
  theme_classic() +
  theme(strip.text.x=element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1),
        axis.title.x = element_blank()) +
  xlab("HDR size (kb)") +
  scale_y_continuous(breaks = seq(0, max(all_calls_LR$divSize/1e3, na.rm = TRUE), by = 25))

s3 <- ggplot(all_calls_LR %>% dplyr::filter(divSize >= 5e3 & divSize <=50e3)) + 
  geom_histogram(aes(x=divSize/1e3),binwidth = 5) + 
  facet_grid(CHROM~threshIDY, scales = "free_x") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1))

idy_bins <- all_LR_stats %>% dplyr::filter(!is.na(bin_IDY)) %>% dplyr::filter(CHROM!="MtDNA")

strain_thresh <- idy_bins %>%
  dplyr::filter(bin_COV>=600) %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::summarise(
    sample_mean = mean(bin_IDY, na.rm = TRUE),
    sample_sd   = sd(bin_IDY, na.rm = TRUE),
    thresh      = sample_mean - 2 * sample_sd,
    upper_q95   = quantile(bin_IDY, probs = 0.05, na.rm = TRUE),
    .groups = "drop")


strain_thresh_all <- idy_bins %>%
  dplyr::filter(bin_COV>=600) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::summarise(
    sample_mean = mean(bin_IDY, na.rm = TRUE),
    sample_sd   = sd(bin_IDY, na.rm = TRUE),
    thresh      = sample_mean - 2 * sample_sd,
    upper_q95   = quantile(bin_IDY, probs = 0.05, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(CHROM="All")

order_levels <- strain_thresh_all %>%
  dplyr::arrange(thresh) %>%
  pull(STRAIN)

plot_df <- rbind(strain_thresh, strain_thresh_all) %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = order_levels))

s4 <- ggplot(plot_df,
             aes(x = STRAIN, y = thresh, group = CHROM, color = CHROM)
) +
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(axis.text = element_text(angle = 45, hjust = 1)) +
  ylab("Bin identity\nlower 5% quantile") +
  xlab("Strain") +
  scale_color_manual(
    values = c("I"   = "#1b9e77",
               "II"  = "#d95f02",
               "III" = "#7570b3",
               "IV"  = "#e7298a",
               "V"   = "#66a61e",
               "X"   = "#e6ab02",
               "All" = "black")
  ) +
  geom_hline(yintercept = 97, linetype="dashed")

comp_reg<-cowplot::plot_grid(s3,s2,s1,s4,ncol=1, rel_heights = c(1.1,1,1,1),align = "v",axis = "lr",labels=c("a","b","c","d"))

ggsave(plot = comp_reg, filename = "../../figures/FigureS28_LR_STATS_95idy_20251201.png",width = 7.5,height = 9,dpi = 600,device = 'png')

# 
# plot_strain_identity <- function(strain_OI, chromlist, rect_colour="blue", rect_linewidth=0.3,xlim=c(NA,NA)) {
#   
#   idy_sub <- idy_bins %>%
#     dplyr::filter(CHROM %in% chromlist & STRAIN==strain_OI)
#   
#   calls_sub <- all_calls_LR %>%
#     dplyr::filter(CHROM %in% chromlist & STRAIN==strain_OI)
#   
#   # shared X limit
#   idy_x <- (idy_sub$START_BIN/1e6) + ((idy_sub$END_BIN/1e6)-(idy_sub$START_BIN/1e6))/2
#   call_x <- calls_sub$maxEnd/1e6
#   x_max <- max(c(idy_x, call_x), na.rm=TRUE)
#   
#   idy_plot <- ggplot()+
#     geom_point(data=idy_sub, aes(x=idy_x, y=bin_IDY)) +
#     geom_hline(yintercept = 95,linetype="dashed")+
#     facet_wrap(~CHROM, nrow=2, scales='free_x') +
#     theme_bw() +
#     theme(axis.title.x=element_blank()) +
#     scale_y_continuous(breaks=seq(floor(min(idy_sub$bin_IDY,na.rm=TRUE)), ceiling(max(idy_sub$bin_IDY,na.rm=TRUE)), by=1)) +
#     scale_x_continuous(limits=c(0,x_max), expand=c(0.01,0)) +
#     ylab("Bin identity") +
#     ggtitle(paste("Chromosome-wide identity -", strain_OI)) +
#     coord_cartesian(xlim = xlim)
#   
#   call_plot <- ggplot()+
#     geom_rect(data=calls_sub, aes(xmin=minStart/1e6, xmax=maxEnd/1e6, ymin=threshIDY-0.2, ymax=threshIDY+0.2),
#               colour=rect_colour, linewidth=rect_linewidth) +
#     facet_wrap(~CHROM, nrow=2, scales='free_x') +
#     scale_y_continuous(breaks=seq(floor(min(calls_sub$threshIDY,na.rm=TRUE)), ceiling(max(calls_sub$threshIDY,na.rm=TRUE)), by=1)) +
#     scale_x_continuous(limits=c(0,x_max), expand=c(0.01,0)) +
#     theme_bw() +
#     xlab("Physical position") +
#     ylab("Identity threshold") +
#     ggtitle(paste("HDRs -", strain_OI)) +
#     coord_cartesian(xlim = xlim) 
#   
#   cowplot::plot_grid(idy_plot, call_plot, nrow=2, align="v", axis="lr")
# }
# 
# q5 <- idy_bins |>
#   dplyr::group_by(CHROM) |>
#   dplyr::summarize(q05 = stats::quantile(bin_IDY, 0.05))
# 
q1sd <- idy_bins |>
  dplyr::group_by(CHROM) |>
  dplyr::summarize(
    cutoff = base::mean(bin_IDY) - 2*stats::sd(bin_IDY))

ggplot2::ggplot(idy_bins, ggplot2::aes(x = bin_IDY, color = CHROM)) +
  ggplot2::geom_freqpoly(binwidth = 0.1, linewidth = 1) +
  ggplot2::geom_vline(
    data = q1sd,
    ggplot2::aes(xintercept = cutoff, color = CHROM),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    breaks = base::seq(
      from = 75,
      to   = 100,
      by   = 1
    )
  ) +
  ggplot2::theme_classic()
# 
# ggplot2::ggplot(idy_bins, ggplot2::aes(x = bin_IDY, color = CHROM)) +
#   ggplot2::geom_freqpoly(binwidth = 0.01, linewidth = 1) +
#   ggplot2::geom_vline(
#     data = q5,
#     ggplot2::aes(xintercept = q05, color = CHROM),
#     linetype = "dashed",
#     linewidth = 0.8,
#     show.legend = FALSE
#   ) +
#   ggplot2::theme_classic() +
#   ggplot2::labs(x = "bin_IDY", y = "Frequency") +
#   ggplot2::scale_x_continuous(
#     breaks = base::seq(
#       from = 75,
#       to   = 100,
#       by   = 1
#     )
#   ) +
#   ggplot2::scale_color_brewer(palette = "Dark2")
# 
# plot_strain_identity(strain_OI="QG2899",chromlist=c("III"))
# plot_strain_identity(strain_OI="JU1373",chromlist=c("IV"))
# plot_strain_identity(strain_OI="ECA1307",chromlist=c("X"),xlim=c(4,5))
####### CALL SR BASED HDRS (FOR LR STRAINS) ##########
covfrac <- c(seq(0.05,0.9,0.05))
varct <- c(seq(5,25,1))

#all strain coverage data
coverage_df <- readr::read_table("../../processed_data/mosdepth_coverage/thresh/CT_all_thresh_cov.tsv",col_names = c("CHROM","START_BIN","END_BIN","NAME","c1X","c2X","c5X","STRAIN"))
#all strain variant count data
varct_df <- readr::read_table("../../processed_data/variant_counts/CT_all_strain_vc.tsv", col_names = c("CHROM","START_BIN","END_BIN","COUNT","STRAIN"))

#join coverage and variant count data
SR_stats_all <- varct_df %>%
  dplyr::left_join(coverage_df,by=c("STRAIN","CHROM","START_BIN","END_BIN")) %>%
  dplyr::select(-NAME) %>%
  dplyr::filter(!(CHROM=="MtDNA")) %>%
  dplyr::mutate(c1X=ifelse(is.na(c1X),0,c1X)) %>%
  dplyr::mutate(c2X=ifelse(is.na(c2X),0,c2X)) %>%
  dplyr::mutate(c5X=ifelse(is.na(c5X),0,c5X)) %>%
  dplyr::mutate(pc1X=c1X/1e3,pc2X=c2X/1e3,pc5X=c5X/1e3) %>%
  dplyr::filter(STRAIN %in% div_str)

#iterate through range of coverage and variant count thresholds
#at each iteration, classify bins and cluster adjacent bins
SR_list <- list()
ct=1
for (i in 1:length(covfrac)) {
  for (k in 1:length(varct)){
    print(paste0("cf:",covfrac[i]," / vc:",varct[k]))
    #classify
    all_stats <- SR_stats_all %>%
      dplyr::mutate(div_class=ifelse(pc1X < covfrac[i],"C", 
                                     ifelse(COUNT >= varct[k],"I","R"))) %>% 
      dplyr::mutate(div=ifelse(div_class=="C" | div_class=="I","div","nondiv")) %>%               
      dplyr::group_by(STRAIN,CHROM) %>%
      dplyr::mutate(div_gf=ifelse(div=="nondiv" & lead(div)=="div" & lag(div)=="div","div",div)) %>%
      dplyr::mutate(div_class=ifelse(div==div_gf,div_class,"G")) %>%
      dplyr::ungroup() 
    
    #cluster
    regList <- list()
    #binList <- list()
    for (j in 1:length(div_str)) {
      temp <- all_stats %>% dplyr::filter(STRAIN==div_str[j])
      div_call <- clusterBins(temp,"SR")
      # div_call[[1]]$STRAIN <- div_str[[j]]
      # div_call[[1]]$CFT <- covfrac[[i]]
      # div_call[[1]]$VCT <- varct[[k]]
      div_call[[2]]$STRAIN <- div_str[[j]]
      div_call[[2]]$CFT <- covfrac[[i]]
      div_call[[2]]$VCT <- varct[[k]]
      #binList[[i]] <- div_call[[1]]
      regList[[j]] <- div_call[[2]]
    }
    
    SR_list[[ct]] <- ldply(regList, data.frame) 
    ct = ct + 1
  }
}
all_SR_meta_calls <- ldply(SR_list,data.frame)



#based on Daehan's paper, we use 95% IDY LR calls as the truth set
t97_LR_calls <- all_calls_LR %>%
  dplyr::filter(threshIDY==97) %>%
  dplyr::select(CHROM,minStart,maxEnd,STRAIN,regionClass) %>%
  dplyr::rename(chrom=CHROM,start=minStart,end=maxEnd)

#iterate through range of coverage and variant count thresholds
#at each iteration, we identify intersects between LR and SR HDR calls
#we estimate overlap fraction, excess fraction, recall, and precision
#we correct for 1:Many overlaps by grouping intersects by the LR call, and aggregating their overlap and excess fractions
#for Many:1, we keep the SR call with the longest overlap
interList_JC <- list()
ct=1
for (i in 1:length(covfrac)) {
  for (k in 1:length(varct)){
    print(paste0("cf:",covfrac[i]," / vc:",varct[k]))
    temp <- all_SR_meta_calls %>%
      dplyr::filter(VCT==varct[k] & CFT==covfrac[i])
    
    for (j in 1:length(div_str)) {
      #get SR calls
      temp_SR <- temp %>%
        dplyr::filter(STRAIN==div_str[j]) %>%
        dplyr::select(CHROM, minStart,maxEnd,STRAIN,VCT,CFT) %>%
        dplyr::rename(chrom=CHROM,start=minStart,end=maxEnd) %>%
        dplyr::filter(end-start >= 5e3)
      
      #get LR calls
      temp_LR <- t97_LR_calls  %>% 
        dplyr::filter(STRAIN==div_str[j]) %>%
        dplyr::filter(end-start >= 5e3)
      
      #call intersects
      interSet <- valr::bed_intersect(temp_SR,temp_LR) %>%
        dplyr::mutate(LR_size=end.y-start.y) %>%
        dplyr::mutate(SR_size=end.x-start.x) %>%
        dplyr::mutate(overlap_fraction=.overlap/LR_size) %>% #overlap fraction (OF)
        dplyr::mutate(exc_fraction=(SR_size-.overlap)/SR_size) %>% # excess fraction (EF)
        dplyr::mutate(LRID=paste(chrom,start.y,end.y,sep = "_")) %>% #generate a unique ID for each LR call
        dplyr::mutate(SRID=paste(chrom,start.x,end.x,sep = "_")) %>% #generate a unique ID for each SR call
        dplyr::group_by(LRID) %>% #group by LR call
        dplyr::mutate(LRgsize=n(),
                      agg_OF=sum(.overlap)/LR_size,
                      agg_EF=(sum(SR_size)-sum(.overlap))/sum(SR_size),
                      locmax_OF=ifelse(.overlap==max(.overlap),T,F)) %>% #correct OF and EF, flag longest overlap
        dplyr::ungroup() %>%
        dplyr::filter(locmax_OF==T) %>% #keep one entry per group (longest overlap)
        dplyr::group_by(SRID) %>% #group by SR call
        dplyr::mutate(SRgsize=n(), 
                      locmax_SRo=ifelse(.overlap==max(.overlap),T,F)) %>% #flag longest overlap
        dplyr::ungroup() %>%
        dplyr::filter(locmax_SRo==T) #keep one entry per group (longest overlap)
      
      interSet$noverRec <- nrow(interSet)/nrow(temp_LR) #recall
      interSet$noverPre <- nrow(interSet)/nrow(temp_SR) #precision
      interSet$noverF1 <- 2*((nrow(interSet)/nrow(temp_SR) * nrow(interSet)/nrow(temp_LR))/(nrow(interSet)/nrow(temp_SR) + nrow(interSet)/nrow(temp_LR)))
      
      interList_JC[[ct]] <-interSet
      
      ct=ct+1
    }
  }
}

#aggregate all intersect data
all_intersects <- ldply(interList_JC,data.frame) %>%
  dplyr::mutate(PP=ifelse(VCT.x <10,paste0("0",VCT.x,":",CFT.x),paste0(VCT.x,":",CFT.x))) %>% #set parameter pair (PP) ID
  dplyr::filter(overlap_fraction > 0) #placeholder, could be used to filter low quality overlaps

############ PARAMETER PAIR SELECTION ##########
statSumm <- all_intersects %>%
  dplyr::group_by(STRAIN.x,PP) %>%
  dplyr::mutate(meanOF=mean(agg_OF)) %>%
  dplyr::mutate(meanEF=mean(agg_EF)) %>%
  dplyr::mutate(gid=cur_group_id()) %>%
  dplyr::ungroup()

maxOF <- statSumm %>%
  dplyr::group_by(gid) %>%
  dplyr::arrange(desc(meanOF)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gid,.keep_all = T)

MOF <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=meanOF,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Mean overlap fraction") +
  geom_line(aes(x=VCT.x,y=meanOF,color=as.factor(CFT.x*100),group=CFT.x))

#meanEF
MEF <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=meanEF,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Mean excess fraction") +
  geom_line(aes(x=VCT.x,y=meanEF,color=as.factor(CFT.x*100),group=CFT.x)) 

REC <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=noverRec,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Recall") +
  geom_line(aes(x=VCT.x,y=noverRec,color=as.factor(CFT.x*100),group=CFT.x)) 

PRE <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=noverPre,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  #geom_vline(xintercept = 14) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("Precision") +
  geom_line(aes(x=VCT.x,y=noverPre,color=as.factor(CFT.x*100),group=CFT.x)) 

F1 <- ggplot(maxOF) + geom_point(aes(x=VCT.x,y=noverF1,color=as.factor(CFT.x*100)),size=0.5)  +
  facet_wrap(~STRAIN.x) +
  #geom_vline(xintercept = 14) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill=NA)) +
  guides(color=guide_legend(title="Percent\nbases\ncovered")) +
  xlab("Variant count") +
  ylab("F1") +
  geom_line(aes(x=VCT.x,y=noverF1,color=as.factor(CFT.x*100),group=CFT.x)) 

topHits <- statSumm %>%
  dplyr::group_by(gid) %>%
  dplyr::arrange(desc(noverF1)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gid,.keep_all = T) %>%
  dplyr::group_by(STRAIN.x) %>%
  dplyr::arrange(desc(noverF1)) %>%
  dplyr::slice_max(order_by = noverF1, n = 19) #%>% #n is increased from N=1 manually to find consensus optimal

bestHit <- topHits %>%
  dplyr::group_by(PP) %>%
  dplyr::mutate(freqPP=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(maxfreq=max(freqPP)) %>%
  dplyr::filter(freqPP==max(freqPP)) %>%
  dplyr::mutate(class="HIFREQ") %>%
  dplyr::select(-maxfreq)

worseHits <- topHits %>%
  dplyr::filter(!(gid %in% bestHit$gid)) %>%
  dplyr::mutate(class="TOPHIT") %>%
  dplyr::group_by(STRAIN.x) %>%
  dplyr::mutate(class=ifelse(noverF1==max(noverF1),"Strain\noptimal",class)) %>%
  dplyr::ungroup()

BEST <- ggplot() + 
  geom_point(data = worseHits %>% filter(class == "Strain\noptimal"),
             aes(x = VCT.x, y = CFT.x, color = meanOF, shape = class), size = 2) +
  geom_point(data = bestHit %>% mutate(class = "Consensus\noptimal\n"),
             aes(x = VCT.x, y = CFT.x, color = meanOF, shape = class), size = 2) +
  facet_wrap(~STRAIN.x) +
  xlab("Variant count") + ylab("Percent bases covered") +
  scale_color_gradient(name = "Mean\noverlap\nfraction", low = "orange", high = "blue") +
  scale_shape_manual(values = c("Strain\noptimal" = 15, "Consensus\noptimal\n" = 17)) +
  guides(shape = guide_legend(title = NULL)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))


####################### CALL SR HDRS SPECIES-WIDE ###############################
#get list of WI (generated from VCF)
strainL <- readLines("processed_data/sample_list/CT_all_samples.txt")
COV_thresh = 0.9
VC_thresh = 9

#merge bin-based variant counts and coverage data
SR_stats_WI <- varct_df %>%
  dplyr::left_join(coverage_df,by=c("STRAIN","CHROM","START_BIN","END_BIN")) %>%
  dplyr::select(-NAME) %>%
  dplyr::filter(!(CHROM=="MtDNA")) %>%
  dplyr::mutate(c1X=ifelse(is.na(c1X),0,c1X)) %>%
  dplyr::mutate(c2X=ifelse(is.na(c2X),0,c2X)) %>%
  dplyr::mutate(c5X=ifelse(is.na(c5X),0,c5X)) %>%
  dplyr::mutate(pc1X=c1X/1e3,pc2X=c2X/1e3,pc5X=c5X/1e3) 


#classify bins
all_stats <- SR_stats_WI %>%
  dplyr::mutate(div_class=ifelse(COUNT >= VC_thresh,"I", #over variant count thresh
                                 ifelse(pc1X < COV_thresh,"C","R"))) %>%  #under coverage thresh
  dplyr::mutate(div=ifelse(div_class=="C" | div_class=="I","div","nondiv")) %>%               
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(div_gf=ifelse(div=="nondiv" & lead(div)=="div" & lag(div)=="div","div",div)) %>% #fill 1kb gaps
  dplyr::mutate(div_class=ifelse(div==div_gf,div_class,"G")) %>% #assign a classification to gaps
  dplyr::ungroup()

#cluster bins (no frequency filter)
regList <- list()
binList <- list()
for (i in 1:length(strainL)) {
  print((i/length(strainL))*100)
  temp <- all_stats %>% dplyr::filter(STRAIN==strainL[i])
  div_call <- clusterBins(temp,"SR")
  #div_call[[1]]$STRAIN <- strainL[[i]]
  div_call[[2]]$STRAIN <- strainL[[i]]
  #binList[[i]] <- div_call[[1]]
  regList[[i]] <- div_call[[2]]
}

all_calls_SR <- ldply(regList, data.frame)

#flag adjacent regions that are 5kb apart
gap_clust <- all_calls_SR %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()

#get flagged and merged them
joinClust <- gap_clust %>% 
  dplyr::filter(dec3g=="join") %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
  tidyr::fill(segbreak,.direction = 'up') %>%
  dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  #dplyr::mutate(gid2=cur_group_id()) #%>%
  dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gapFoot=paste0(rep("R",forGapSize/1e3),collapse = "")) %>%
  dplyr::mutate(new_foot=ifelse(flag3g=="clust",paste0(bin_foot,gapFoot),bin_foot)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(clust_foot=paste0(new_foot,collapse = "")) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(newMeanVC=mean(meanVC)) %>%
  dplyr::mutate(newMeanCF=mean(meanCF)) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize,-meanVC,-meanCF,-bin_foot) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize,meanVC=newMeanVC,meanCF=newMeanCF,bin_foot=clust_foot) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN,nclust)

#keep unflagged 
nojoin <- gap_clust %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,meanVC,meanCF,bin_foot,STRAIN) %>%
  dplyr::mutate(nclust=1)

#bind unflagged and merged calls, arrange strains by number of HDRs
all_calls_SR_clustered <- rbind(joinClust,nojoin) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() 

############# PLOT and WRITE

p1 <- ggplot(all_calls_SR_clustered %>% dplyr::filter(divSize >= 5e3 & !(CHROM=="MtDNA"))) + 
  geom_rect(aes(xmin=minStart/1e6,xmax=maxEnd/1e6,ymin=rleID-0.45,ymax=rleID+0.45)) + 
  facet_wrap(~CHROM,scales = 'free',nrow = 1,ncol = 6) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank())  +
  xlab("Physical position (Mb)") +
  ylab("622 Isotype strains") +
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 5),expand = c(0, 0))

ggsave(plot = p1, filename = "../../figures/FigureSx_HDR_CT_allStrain_5kbclust_20251201.png",width = 8.5,height = 3.5,dpi = 600,device = 'png')
ggsave(plot = MOF, filename = "../../figures/FigureS29_MOF_HDR_CT_20251201.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = MEF, filename = "../../figures/FigureS30_MEF_HDR_CT_20251201.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = REC, filename = "../../figures/FigureS31_REC_HDR_CT_20251201.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = PRE, filename = "../../figures/FigureS32_PRE_HDR_CT_20251201.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = F1, filename = "../../figures/FigureS33_F1_HDR_CT_20251201.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')
ggsave(plot = BEST, filename = "../../figures/FigureS34_BEST_HDR_CT_20251201.png",width = 7.5,height = 6.5,dpi = 600,device = 'png')

options(scipen=10)
write.table(all_calls_SR_clustered, file="../../tables/HDR_CT_allStrain_5kbclust_20251201.tsv",row.names = F,quote = F,sep = '\t')
options(scipen=0)

#save.image("Trop_HDR_checkpoint_20251201.Rda")
