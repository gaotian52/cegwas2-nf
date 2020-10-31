#!/usr/bin/env Rscript
library(tidyverse)


# load arguments
args <- commandArgs(trailingOnly = TRUE)


# load genotype matrix
pheno_ori <- read.delim(args[1], stringsAsFactors=FALSE) %>%
  na.omit()

phenotype_data <- pheno_ori

trait <- colnames(phenotype_data[2])

for (i in 1:200) {
  sss<-base::sample(phenotype_data[,2], length(phenotype_data[,2]))
  
  ss_h<-phenotype_data$strain
  
  ss_n<-glue::glue("{trait}_p{i}")
  
  ss_per <- data.frame(strian=ss_h,ss=sss) 
  
  colnames(ss_per) <- c("strain", ss_n)
  

#write.table(ss_per, paste("pr_",trait,"p",i,".tsv",sep=""), sep = "\t", row.names = F, quote = FALSE)
  
write.table(ss_per, 
            file = glue::glue("pr_{ss_n}.tsv"),
            quote = F, col.names = T, row.names = F, sep="\t")

}

