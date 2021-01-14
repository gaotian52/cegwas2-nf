#!/usr/bin/env Rscript
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)


 med_analysis <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE)
 

 med_analysis_head <- med_analysis %>% 
   #dplyr::filter(V15=="Threshold")%>% 
   dplyr::filter(V1=="metabolite")


 colnames(med_analysis) <- med_analysis_head[1,]
 
 
 pr_med_df <- med_analysis %>% 
   dplyr::filter(!metabolite=="metabolite") %>% 
   dplyr::mutate(abs_est=as.numeric(S),
                 prob=as.numeric(p),
                 m_start=as.numeric(m_start),
                 m_peak=as.numeric(m_peak),
                 m_end=as.numeric(m_end),
                 m_log10p=as.numeric(m_log10p),
                 m_var_exp=as.numeric(m_var_exp),
                 gwpeak=as.numeric(gwpeak)) %>% 
   dplyr::select(-S,-p) %>% 
   dplyr::mutate(gwas_qtl=paste(gwtrait,m_chr,gwpeak,sep="_"))
 


med <- pr_med_df

# save mapping data set
readr::write_tsv(med, 
                 path = glue::glue("mediation_analysis_med.tsv"),
                 col_names = T)




med_transcript <- med 

med_plt <- ggplot(med_transcript) +
  aes(x = m_peak/1e6, y = abs_est)+#,color=eQTL_classification) +
  geom_point(aes(alpha = prob)) +
  scale_alpha_continuous(range = c(1, 0.1)) +  
 # scale_color_manual(values = c("distant"="gray69","local"="gold2"))+
  theme_bw(10) + 
   facet_grid(gwtrait~m_chr, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.1, size=12),
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        axis.title =  element_text(size=12),
        strip.text = element_text(size=12, vjust = 1),
        strip.background = element_blank(),
        axis.text = element_text(size=12) ) +
  labs(x = "mQTL peak position (Mb)", y = "Mediation estimate") 

ggsave(med_plt,
       filename = ("mediation_est_plt_metabolites.png"))






med05 <- med %>% dplyr::filter(prob<0.1)

if (nrow(med05) > 0) {
med_plt_05 <- ggplot(med05) +
  aes( y=gwtrait, x = abs_est,color=metabolite) +
  geom_point(aes(alpha = prob)) +
 # geom_text_repel(aes(label = mediator_gene), fontface="italic")+
  scale_alpha_continuous(range = c(1, 0.1)) +  
#  scale_color_manual(values = c("gray69","gold2"))+
  theme_bw(10) + 
  facet_grid(.~m_chr, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.1, size=12),
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        axis.title =  element_text(size=12),
        strip.text = element_text(size=12, vjust = 1),
        strip.background = element_blank(),
        axis.text.y = element_text(size=10) ,
        axis.text.x =  element_text(size=10 ,angle = 90, hjust = 1),
        axis.title.y=element_blank()) +
  labs( x = "Mediation estimates")




ggsave(med_plt_05,
       filename = ("mediation_est_plt_05.png"))
}