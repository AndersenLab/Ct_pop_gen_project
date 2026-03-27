rm(list=ls())

library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(viridis)
library(pheatmap)
library(gplots)
library(gridExtra)
library(grid)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(phytools)
library(geosphere)
library(igraph)
library(ggraph)
library(ggforce)
library(ComplexHeatmap)
library(circlize)
library(admixtools)

source("../utilities.R")

geo_colors <- geo.colours
df_colors <- data.frame(unname(geo_colors),names(geo_colors)) %>% dplyr::rename(color=`unname.geo_colors.`, geo=`names.geo_colors.`)

lineages <- read.csv("../../processed_data/geo_info/Ct_lineage_all.csv") %>% 
  rename(isotype=sample) %>% 
  rename(Lineage=lineage)

geo_info <- readr::read_csv(file="../../processed_data/Geo_info/Ct_indep_isotype_info_geo.csv") %>%
  dplyr::left_join(df_colors,by=c("geo")) %>%
  dplyr::left_join(lineages %>% dplyr::select(isotype, Lineage), by="isotype") %>%
  dplyr::filter(!is.na(Lineage)) %>%
  dplyr::mutate(abslat=abs(lat))

cvmat <- readr::read_tsv(file="../../processed_data/Ct_admixture_best_k/cv_matrix.tsv")

cvmat_long <- cvmat %>%
  tidyr::pivot_longer(
    cols = starts_with("S"),
    names_to = "seed",
    values_to = "cv_error",
    values_drop_na = T) %>%
  dplyr::mutate(seed = as.integer(sub("^S", "", seed)))

cv <- ggplot() +
  geom_jitter(data=cvmat_long,aes(x=K,y=cv_error,group=K),alpha=0.7,size=0.5, width = 0.2, height = 0) +
  geom_boxplot(data=cvmat_long,aes(x=K,y=cv_error,group=K),fill=NA,outliers = F) +
  labs(color="Seed") +
  scale_x_continuous(breaks = seq(min(cvmat_long$K), max(cvmat_long$K), by = 1), expand = c(0,0)) +
  ylim(0,NA)+
  theme_bw() +
  ylab("CV")

kmin=2
kmax=30
fraclist <- list()
qlist <- list()
for (i in kmin:kmax) {
  
  if (i>26) {
    colnames = c()
    for (i in 27:i){
      symbol=paste0(LETTERS[floor((i+1)/26)],LETTERS[i-(26*floor((i+1)/26))])
      colnames=c(colnames,symbol)
    }
    colnames=c(LETTERS[1:26],colnames)
  } else {
    colnames=LETTERS[1:i]
    symbol=LETTERS[i]
  }
  
  tmpQ <- readr::read_tsv(paste0("../../processed_data/Ct_admixture_best_k/concat_Qfiles_K",i,".tsv"),col_names=c("isotype",colnames,"run_ID"))
  
  qlist[[i-(kmin-1)]] <- tmpQ %>%
    tidyr::pivot_longer(
      cols = A:!!sym(symbol),
      names_to = "subpop",
      values_to = "fraction") %>%
    tidyr::separate(run_ID, into=c("LD_tag","LD_val","K_val","seed"),sep = "_") %>%
    dplyr::select(-LD_tag) %>%
    dplyr::left_join(geo_info,by="isotype")
  
  
  fraclist[[i-(kmin-1)]] <- qlist[[i-(kmin-1)]] %>%
    dplyr::filter(fraction >= 0.999) #get non-admixed individuals
}

ctlist<- list()
rgm_perseed <- list()
for (i in kmin:kmax) {
  RGassignment<-fraclist[[i-(kmin-1)]] %>%
    dplyr::group_by(K_val,seed) %>%
    dplyr::mutate(nlin_assigned=length(unique(Lineage))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(K_val,seed,subpop) %>%
    dplyr::mutate(nlin=length(unique(Lineage))) %>%
    dplyr::distinct(Lineage,.keep_all = T) %>%
    dplyr::select(-isotype) %>%
    dplyr::ungroup()
  
  mergedRG<- RGassignment %>%  dplyr::filter(nlin>1)
  
  ctlist[[i-(kmin-1)]] <- c(i,length(unique(mergedRG$seed)),max(mergedRG$nlin_assigned),max(mergedRG$nlin))
  
  per_seed <- RGassignment %>%
    dplyr::group_by(K_val,seed) %>%
    dplyr::filter(nlin==max(nlin)) %>%
    dplyr::distinct(seed,.keep_all = T)
  
  rgm_perseed[[i-(kmin-1)]] <- per_seed$nlin
}
df <- as.data.frame(do.call(rbind, ctlist)) %>% dplyr::filter(V1>=10)
df_perseed <- do.call(rbind, lapply(seq_along(rgm_perseed), function(i) {data.frame(K   = i + (kmin - 1), num = rgm_perseed[[i]])})) %>% dplyr::filter(K>=10)
colnames(df) <- c("K", "nrun_with_merged_RG", "max_nonadmix_RGs", "max_nonadmix_RG_merged")

rg1 <- ggplot()+geom_point(data=df,aes(x=K,y=nrun_with_merged_RG),size=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7)) +
  ylab("Number of independent runs\nwith multiple RGs per subpopulation") +
  scale_x_continuous(breaks = seq(min(df$K), max(df$K), by = 1)) +
  scale_y_continuous(breaks = seq(min(df$nrun_with_merged_RG), max(df$nrun_with_merged_RG), by = 1))

rg2 <- ggplot()+geom_point(data=df,aes(x=K,y=max_nonadmix_RGs),size=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7)) +
  ylab("Maximum number of  RGs\nwith non-admixed individuals") +
  scale_x_continuous(breaks = seq(min(df$K), max(df$K), by = 1)) +
  scale_y_continuous(breaks = seq(min(df$max_nonadmix_RGs), max(df$max_nonadmix_RGs), by = 1))

df_counts <- df_perseed %>%
  count(K, num)

rg3<-ggplot(df_counts, aes(x = K, y = num, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "black",size=2.5) +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5),
        axis.text = element_text(size=6),
        axis.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7),
        legend.justification.right = "top",
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(0,0,0,0))+
  labs(fill = "Number of\nindependent\nruns") +
  ylab("Maximum number of RGs\nper subpopulation") +
  scale_x_continuous(breaks = seq(min(df_perseed$K), max(df_perseed$K), by = 1)) +
  scale_y_continuous(breaks = seq(min(df_perseed$num), max(df_perseed$num), by = 1))

comprg <- cowplot::plot_grid(rg2,rg1,rg3,
                             nrow=1,
                             ncol=3,
                             rel_widths = c(0.9,0.9,1.2), 
                             align="h",
                             axis = "tb", 
                             labels=c("b","c","d"))
rg_cv <- cowplot::plot_grid(cv + theme(panel.grid.major = element_line(color="grey80"),panel.grid.minor = element_blank()),comprg,nrow=2,rel_heights = c(1.2,1), align = "v",axis="r",labels=c("a",NA))
ggsave(rg_cv,file="../../figures/FigureS8_admixture_CV_error.pdf",
       width = 7,height = 6,
       units = "in",device = 'pdf',
       bg="white",dpi=600)

## Set 28 as the best K
### lowest CV error at K=28
cvmat %>% filter(K==28) %>% min()
# [1] 0.05033

cvmat %>% filter(K == 28) %>%   
  select(where(~ .x == cvmat %>% filter(K==28) %>% min())) %>%
  names()

