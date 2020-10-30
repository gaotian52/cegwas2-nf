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

eqtl_transcript <- read.delim(args[7], stringsAsFactors=FALSE)

# eQTL overlap GWAS peak
gene_qtl <- eqtl_transcript %>% 
  #bind_rows(eqtl_int,eqtl_transcript) %>% 
  dplyr::filter(e_chr==gwas_pchr) %>% 
  dplyr::filter(!e_end<=(gwas_start-1e6)) %>% 
  dplyr::filter(!e_start>=(gwas_stop+1e6)) %>% 
  dplyr::mutate(gwtrait=gtrait, gwpeak=gwas_p) #%>% dplyr::filter(Threshold=="FDR_BF") %>% dplyr::filter(logP>60)

# save mapping data set
readr::write_tsv(gene_qtl, 
                 path = glue::glue("{gtrait}_{gwas_pchr}_{gwas_p}_eqtl.tsv"),
                 col_names = T)



gene_qtl_genelist <- gene_qtl %>% 
  dplyr::select(gwtrait,e_chr,gwpeak,trait) %>% 
  distinct()


# save gene list 
readr::write_tsv(gene_qtl_genelist, 
                 path = glue::glue("{gtrait}_{gwas_pchr}_{gwas_p}_elist.tsv"),
                 col_names = F)
