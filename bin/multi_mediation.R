.libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths() ))

#!/usr/bin/env Rscript

library(tidyverse)
library(MultiMed)
# load arguments
args <- commandArgs(trailingOnly = TRUE)


# load genotype matrix
Genotype_Matrix <- readr::read_tsv(args[1]) %>%
  na.omit()

# load pheno data
trait_phenotype <- read.delim(args[3], stringsAsFactors=FALSE) 




# load expression data
# gene level 
 # expression_pheno_raw <- read.delim(args[2], stringsAsFactors=FALSE)
  
 # expression_pheno <- expression_pheno_raw %>% 
 #   dplyr::filter(strain %in% trait_phenotype$strain)
  

# transcript level
  
  texpression_pheno_raw <- read.delim(args[2], stringsAsFactors=FALSE)
  
  texpression_pheno <- texpression_pheno_raw %>% 
    gather(trait2,value,-strain) %>% 
    dplyr::mutate(trait=sub("(^X)(.*)","\\2",trait2)) %>% 
    dplyr::select(strain,trait,value) %>% spread(trait,value) %>% 
    dplyr::filter(strain %in% trait_phenotype$strain)
  


  # processed pheno data
  trait_pheno <- trait_phenotype%>% 
    dplyr::rename(trait=tr) %>% 
    dplyr::filter(strain %in% texpression_pheno$strain)




# GWAS qtl infor
gwas_intchr = args[4]

gwas_peak = args[5] %>% as.numeric()

gwtrait = args[6]

# load eqtl data
eqtl_infor <- read.delim(args[7], stringsAsFactors=FALSE)

gene_qtl_genelist <- eqtl_infor %>% 
  dplyr::select(gwtrait,e_chr,gwpeak,trait) %>% 
  distinct()


gene_list <- gene_qtl_genelist %>% 
  dplyr::filter(grepl("WBGene",trait))


transcript_list <- gene_qtl_genelist %>% 
  dplyr::filter(!grepl("WBGene",trait))



# get the genotype at the peak marker
gwas_g <- Genotype_Matrix %>% 
  dplyr::filter(CHROM==gwas_intchr & POS == gwas_peak) %>% 
  dplyr::select(-(1:4)) %>% 
  tidyr::gather(strain,geno) %>% 
  dplyr::filter(strain %in% texpression_pheno$strain) %>% 
  dplyr::arrange(strain)



#head(df_multi_gene2)

# get the transcript expression data
t_lgmtpm_gwas <- texpression_pheno %>% 
  dplyr::select(c(strain,unique(transcript_list$trait))) %>% 
  dplyr::arrange(strain)  %>% 
  dplyr::select(-strain)

if( length(unique(gwas_g$geno)) == 2 ){
  
#pick transcripts with variation in expression
t_lgmtpm_gwas_vari <- t_lgmtpm_gwas[vapply(t_lgmtpm_gwas, function(x) length(unique(x)) > 1, logical(1L))]


exp_matr_transcript <- as.matrix(t_lgmtpm_gwas_vari)


mt_multi_transcript <- medTest(gwas_g$geno, exp_matr_transcript, trait_pheno$trait, nperm = 1000)

df_multi_transcript <- as.data.frame(mt_multi_transcript) 

row.names(df_multi_transcript) <- (colnames(exp_matr_transcript))

df_multi_transcript2 <- df_multi_transcript %>% 
  rownames_to_column(var="gene")  %>% 
  dplyr::arrange(p)


df_multi_med <- df_multi_transcript2 %>% 
  left_join(eqtl_infor,by=c("gene"="trait"))


# save mapping data set
readr::write_tsv(df_multi_med, 
                 path = glue::glue("{gwtrait}_{gwas_intchr}_{gwas_peak}_med.tsv"),
                 col_names = T)

}
