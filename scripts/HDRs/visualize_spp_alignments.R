library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


brig_ex <- readr::read_tsv("../../processed_data/genome_alignments/other_species_divergent/ECA2666xQX1410.nucmer.alignments.tsv",
                           col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI"))%>%
  dplyr::mutate(STRAIN="ECA2666",sp="C. briggsae", comp="QX1410 x ECA2666")

ele_ex <- readr::read_tsv("../../processed_data/genome_alignments/other_species_divergent/XZ1516xN2.nucmer.alignments.tsv",
                          col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI"))%>%
  dplyr::mutate(STRAIN="XZ1516",sp="C. elegans", comp="N2 x XZ1516")

trop_ex <- readr::read_tsv("../../processed_data/genome_alignments/CT_all_nucmer_aln.tsv",
                           col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")) %>%
  dplyr::filter(STRAIN=="ECA1307") %>%
  dplyr::mutate(sp="C. tropicalis", comp="NIC58 x ECA1307")

all_ex <- rbind(trop_ex,ele_ex,brig_ex) %>%
  dplyr::filter(L2>1e3)

compplot<-ggplot(all_ex %>% dplyr::filter(REF == "V")) +
  geom_rect(aes(xmin=S1/1e6, xmax=E1/1e6, ymin=IDY+0.2, ymax=IDY-0.2)) +
  facet_wrap(sp + comp ~ ., ncol = 1,
             scales = "free_x",
             strip.position = "right") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "italic"),
        strip.background = element_blank()) +
  scale_x_continuous(expand = c(0,0))+
  ylab("Genome alignment percent identity") +
  xlab("Reference genome physical position (Mb)") +
  ggtitle("Chromosome V")


ggsave(compplot,filename = "../../figures/FigureS26_REFxDIV_V_alnidy.png",width = 7,height = 5.5,dpi = 600,device = 'png')
