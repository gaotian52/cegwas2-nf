#!/usr/bin/env Rscript
library(tidyverse)
library(rrBLUP)

args <- commandArgs(trailingOnly = TRUE)


save_name <- gsub(".LD.tsv", "", args[4])

cores_avail <- as.numeric(args[5])

complete_matrix <- readr::read_tsv(args[1]) %>%
  na.omit()

kinship_matrix <- rrBLUP::A.mat(t(complete_matrix[,5:ncol(complete_matrix)]), 
                                n.core = cores_avail)

ROI_matrix <- readr::read_tsv(args[2]) %>%
  na.omit()

pr_map <- data.table::fread(args[3])

map_pheno <- na.omit(pr_map) %>%
  dplyr::distinct(strain, value) %>%
  as.data.frame()

map_pheno$value <- as.numeric(map_pheno$value)

ROI_matrix <- ROI_matrix[colnames(ROI_matrix)%in%c("CHROM","POS","REF","ALT",map_pheno$strain)]

kinship_matrix <- kinship_matrix[row.names(kinship_matrix)%in%map_pheno$strain,
                                 colnames(kinship_matrix)%in%map_pheno$strain]


gwa_mapping <- function (data, 
                         cores = cores_avail, 
                         kin_matrix = kinship_matrix, 
                         snpset = NULL, 
                         min.MAF = 0.05,
                         p3d = as.logical(args[6])) {
  x <- data
  
  y <- snpset %>% dplyr::mutate(marker = paste0(CHROM, "_", POS)) %>% 
    dplyr::select(marker, everything(), -REF, -ALT) %>% 
    as.data.frame()
  
  kin <- as.matrix(kin_matrix)
  
  pmap <- rrBLUP::GWAS(pheno = x, 
                       geno = y, 
                       K = kin, 
                       min.MAF = min.MAF, 
                       n.core = cores, 
                       P3D = as.logical(p3d), 
                       plot = FALSE)
  
  return(pmap)
}

roi_mapping <- gwa_mapping(data = map_pheno,
                           snpset = ROI_matrix,
                           kin_matrix = kinship_matrix)

roi_ld <- readr::read_tsv(args[4])

head(roi_mapping)
head(roi_ld)

peakp <- unique(roi_ld$BP_A)

map_peaks <- na.omit(roi_mapping) %>%
  dplyr::filter(POS == unique(roi_ld$BP_A))

map_peaks$POS <- as.numeric(map_peaks$POS)

pr_roi_ld <- roi_ld %>%
  dplyr::mutate(peak_marker = gsub("_", ":", unique(map_peaks$marker))) %>%
  dplyr::mutate(marker = gsub(":", "_", SNP_B)) %>%
  dplyr::select(peak_marker, peak_maf = MAF_A, marker, maf_marker_b = MAF_B, ld_r2 = R2) %>%
  dplyr::left_join(roi_mapping,., by = "marker") %>%
  dplyr::filter(value > 0)

readr::write_tsv(pr_roi_ld, path = glue::glue("{save_name}_prLD_df.tsv"))

