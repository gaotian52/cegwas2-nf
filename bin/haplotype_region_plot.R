
library(tidyverse)


NJ_label_order <- read.delim("NJ_403isotype_order.tsv", stringsAsFactors=FALSE)

swept_chr_count <- read.delim("swept_chr_count.tsv", stringsAsFactors=FALSE)


# load arguments
args <- commandArgs(trailingOnly = TRUE)

processed_mapping <- read.delim(args[1], stringsAsFactors=FALSE)


trait <- unique(processed_mapping$trait)



QTL_peaks <- na.omit(processed_mapping) %>% dplyr::distinct(CHROM,startPOS,endPOS)

write.table(QTL_peaks, paste("QTL_peaks.bed",sep=""), sep = "\t", row.names = F, col.names = F, quote = FALSE)


system(glue::glue("bedtools intersect -wo -a QTL_peaks.bed -b swept_haplotype_1123.bed > QTL_peaks_swept_region.bed"))





QTL_peaks_swept_region <- read.delim(paste("QTL_peaks_swept_region.bed",sep = ""), header=FALSE, stringsAsFactors=FALSE)

colnames(QTL_peaks_swept_region) <- c("CHROM","startPOS","endPOS",
                                      "chromosome","start","stop","haplotype","isotype","sat", "plotpoint","segment","cvalue","hap_length","color","chrom_haplotype_sum",           
                                      "swept_haplotype","swept_haplotype_name","isotype_has_swept_haplotype","isotypes_w_haplotype","isotype_swept_haplotype_length","max_swept_haplotype_length",    
                                      "max_haplotype_shared","filtered_swept_haplotype_len","filtered_sweep_len","filtered_sweep_ratio","is_swept",  
                                      "nnn")





qtl_swept_order <- na.omit(processed_mapping) %>%
  dplyr::distinct(marker, isotype=strain, allele, .keep_all = T)  %>%
  dplyr::select(isotype,CHROM,allele,startPOS,peakPOS,endPOS) %>%
  left_join(NJ_label_order) %>% 
  dplyr::group_by(CHROM,allele) %>%
  mutate(my_ranks = order(order(plotpoint_hi))) %>%
  left_join(QTL_peaks_swept_region) %>%
  dplyr::mutate(start_s=ifelse(is.na(start),startPOS,
                               ifelse(start<startPOS,startPOS,start)))%>%
  dplyr::mutate(stop_s=ifelse(is.na(stop),endPOS,
                              ifelse(stop>endPOS,endPOS,stop))) %>%
  dplyr::mutate(plotpoint_hi=as.numeric(my_ranks)) %>%
  dplyr::mutate(ALL=ifelse(allele==-1,"REF","ALT")) %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>% 
  left_join(swept_chr_count) %>% 
  dplyr::select(isotype,allele,ALL,start_s,stop_s,plotpoint_hi,swept_haplotype,CHROM,peakPOS,facet_marker,genotype) 


write.table(qtl_swept_order, paste(trait,"_haplotype.tsv",sep=""), sep = "\t", row.names = F, quote = FALSE)



for (m in unique(qtl_swept_order$facet_marker)) {

signle_marker <- qtl_swept_order %>% dplyr::filter(facet_marker==m)

signle_marker$alleles <- factor(signle_marker$ALL,levels = c("REF","ALT"))

plt <- ggplot(signle_marker,
       aes(xmin = start_s/1E6, xmax = stop_s/1E6,
           ymin = -(plotpoint_hi + 0.5), ymax = -(plotpoint_hi - 0.5))) +
  geom_rect(show.legend = FALSE,
            aes(fill = swept_haplotype)) +
  scale_fill_manual(values = c("FALSE"="gray","TRUE"="red")) +
  scale_y_continuous(breaks = unique(signle_marker$plotpoint_hi),
                     expand = c(0, 0)) +
  geom_point(aes(x=(peakPOS/1E6),y=-plotpoint_hi,color=genotype),size=0.5)+
  scale_color_manual(values=c("swept"="gold2","divergent"="plum4")) + 
  xlab("Genomic position (Mb)") +
  theme_bw() +
  facet_grid(alleles~facet_marker, scales="free") +
  theme( axis.text.y =  element_blank(),
         axis.text.x =  element_text(size=12),
        axis.ticks.y = element_blank(),
        axis.title.x =  element_text(size=12),
        axis.title.y =  element_blank(),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        panel.grid = ggplot2::element_blank(),
        text=element_text(family="Helvetica"),
        plot.margin = unit(c(0, 2, 0, 10), "mm"),
        panel.spacing = unit(0.3, "lines")) +
  xlab("Genomic position (Mb)")


ggsave(plt, filename = paste(trait,"_", m, "_qtl_haplotype.png",sep=""), units = "mm",height = 100, width = 60)

}




