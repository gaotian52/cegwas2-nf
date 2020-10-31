#!/usr/bin/env Rscript
library(tidyverse)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# permutation significance threshold

fdr_threshold <- args[3]

fdr_threshold <- round(as.numeric(fdr_threshold), digits = 3)


chr_pos <- read.table(args[1], quote="\"", comment.char="", stringsAsFactors=FALSE)



QTL_peaks <- read.delim(args[2], header=FALSE, stringsAsFactors=FALSE)


t2g_all <- read.delim(args[4], stringsAsFactors=FALSE)




trait_pos_empir <- chr_pos %>%
  dplyr::mutate(V5=V3-4717,V7=V4+4717) %>%
  dplyr::mutate(V6=ifelse(V5<0,0,V5)) %>%
  dplyr::rename(trait=V1,Chr=V2,start=V6,end=V7) %>%
  dplyr::select(-V3,-V4,-V5)


qtl_peak2 <- QTL_peaks %>%
  dplyr::rename(trait=V1,chrom=V2,peak=V4) %>% 
  left_join(trait_pos_empir) %>%
  dplyr::mutate(eQTL_classification=ifelse(Chr==chrom & V3< start & V5>end, "local",
                                ifelse(Chr==chrom & V3<start & V5<end & V5>start, "local",
                                       ifelse(Chr==chrom & V3>start & V3<end & V5>end, "local",
                                              ifelse(Chr==chrom & V3>start & V3<end & V5<end, "local","distant"))))) %>% # a local eQTL defined as an eQTL whose confidence interval includes the gene it influences at genome-wide
  #dplyr::filter(!(Chr=="MtDNA" | chrom=="MtDNA")) %>%
  dplyr::select(trait,Chr,start,end,e_chr=chrom,e_start=V3,e_peak=peak,e_end=V5,logP=V6,var_exp=V7,eQTL_classification) %>% 
  dplyr::mutate(e_length=e_end-e_start) %>% 
  left_join(t2g_all,by=c("trait"="transcript"))



write.table(qtl_peak2, paste("eQTL_Peaks_",fdr_threshold,".tsv",sep=""), sep = "\t", row.names = F, quote = FALSE)



qtl_peak_exp05 <- qtl_peak2 %>% dplyr::filter(var_exp>=0.05)

write.table(qtl_peak_exp05, paste("eQTL_Peaks_",fdr_threshold,"_exp05.tsv",sep=""), sep = "\t", row.names = F, quote = FALSE)




qtl_peak <- qtl_peak2 %>% 
  dplyr::filter(!(Chr=="MtDNA" | e_chr=="MtDNA")) 



n_traits <- length(unique(qtl_peak$trait)) 

n_local <- nrow(dplyr::filter(qtl_peak,eQTL_classification=="local")) 

n_distant <- nrow(dplyr::filter(qtl_peak,eQTL_classification=="distant")) 

qtl_peak$Chr_pos<- factor(qtl_peak$Chr,levels = c("X","V","IV","III","II","I"))

plt_eqtlmap <- ggplot(data=qtl_peak)  + 
  geom_rect( aes(xmin = e_start/1E6, xmax = e_end/1E6,
                 ymin = (start/1E6+0.1), ymax = (end/1E6-0.1),
                 fill = eQTL_classification)) +   
  geom_point(aes(x=e_peak/1E6,y=start/1E6),color="red",size=1.2) +
  scale_fill_manual(values = c("distant"="gray69","local"="gold2")) +
  facet_grid(cols=vars(e_chr), rows=vars(Chr_pos), 
             scales = "free", switch="both") +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "grey", fill = NA, size = 1),
        panel.spacing = unit(0,"line"),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        axis.title =  element_text(size=12),
        strip.text = element_text(size=12, vjust = 1),
        strip.background = element_blank())  +
  ggtitle(paste("Traits: ",n_traits,"   ","Local_all: ",n_local,"   ","Distant_all: ",n_distant, " ","Pemutated_5%FDR: ", fdr_threshold, sep="")) +
  ylab("Gene position") + 
  xlab("eQTL position")


ggsave(plt_eqtlmap, filename = paste("eQTL_map_",fdr_threshold,".png",sep = ""), height = 8, width = 12)


