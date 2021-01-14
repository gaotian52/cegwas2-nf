#!/usr/bin/env Rscript
library(tidyverse)

# load arguments
args <- commandArgs(trailingOnly = TRUE)


gtrait = args[1]

# load phenotpe data and scale
phenotype_data <- readr::read_tsv(args[2]) %>%
  na.omit() %>%
  as.data.frame() %>% 
  dplyr::select(strain,ph=gtrait) %>%
  dplyr::mutate(newpheno = (ph - mean(ph, na.rm = T)) / sd(ph, na.rm = T)) %>%
  dplyr::select(strain, tr = newpheno)

# save mapping data set
readr::write_tsv(phenotype_data, 
                 path = glue::glue("{gtrait}_scaled_mapping.tsv"),
                 col_names = T)


gwas_pchr = args[3] 

gwas_start = args[4] %>% as.numeric()

gwas_stop = args[5] %>% as.numeric()

gwas_p = args[6] %>% as.numeric()

#eqtl_int <- read.delim(args[7], stringsAsFactors=FALSE)

qtl_metabolite <- read.delim(args[7], stringsAsFactors=FALSE)

# eQTL overlap GWAS peak
meta_qtl <- qtl_metabolite %>% 
  #bind_rows(eqtl_int,eqtl_transcript) %>% 
  dplyr::filter(m_chr==gwas_pchr) %>% 
  dplyr::filter(!m_end<=(gwas_start-1e6)) %>% 
  dplyr::filter(!m_start>=(gwas_stop+1e6)) %>% 
  dplyr::mutate(gwtrait=gtrait, gwpeak=gwas_p) #%>% dplyr::filter(Threshold=="FDR_BF") %>% dplyr::filter(logP>60)



if( nrow(meta_qtl) > 0 ){

# save mapping data set
readr::write_tsv(meta_qtl, 
                 path = glue::glue("{gtrait}_{gwas_pchr}_{gwas_p}_mqtl.tsv"),
                 col_names = T)

}

