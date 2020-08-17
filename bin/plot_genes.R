#!/usr/bin/env Rscript
library(tidyverse)
library(cegwas2)

# input arguments
# 1 = LD file
# 2 = phenotype file
# 3 = gene file
# 4 = vcf


args <- commandArgs(trailingOnly = TRUE)


pr_trait_ld <- data.table::fread(args[1])
phenotypes <- readr::read_tsv(args[2])

load(file = args[3])

analysis_trait <- colnames(phenotypes)[2]

colnames(phenotypes) <- c("strain", "Phenotype_Value")

query_regions <- pr_trait_ld %>%
  dplyr::select(CHROM, start_pos, end_pos)%>%
  dplyr::distinct()

query_regions

# verify that provided VCF has annotation field
test_out <- system(glue::glue("bcftools view {args[4]} | head -10000 | grep ANN"), intern=T)

if(length(test_out) > 0) {
    q_vcf <- args[4]
} else {
    system("echo 'Provided VCF does not have an ANN column, using CeNDR default'")
    if(grepl("20180527", args[5])){
        q_vcf <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/{args[5]}/variation/WI.{args[5]}.soft-filter.vcf.gz")
    } else {
        q_vcf <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/{args[5]}/variation/WI.{args[5]}.soft-filtered.vcf.gz")
    }
}

snpeff_out <- list()
for(r in 1:nrow(query_regions)){
  cq <- query_regions$CHROM[r]
  sq <- query_regions$start_pos[r]
  eq <- query_regions$end_pos[r]

  snpeff_out[[r]] <- cegwas2::query_vcf(glue::glue("{cq}:{sq}-{eq}"),
                                        impact = "ALL",
                                        samples = unique(phenotypes$strain),
                                        vcf = q_vcf) %>%
    dplyr::select(CHROM:ALT, strain = SAMPLE, allele = a1, effect:transcript_biotype, nt_change:aa_change)%>%
    dplyr::distinct(CHROM, POS, strain, REF, ALT, .keep_all = T) %>%
    tidyr::unite(marker, CHROM, POS, sep = "_")
}

snpeff_df <- dplyr::bind_rows(snpeff_out) %>%
  dplyr::left_join(pr_trait_ld, ., by = "marker")

genes_in_region <- gene_ref_flat %>%
  dplyr::filter(wbgene %in% snpeff_df$gene_id) %>%
  dplyr::select(gene_id = wbgene, strand, txstart, txend, feature_id = gene) %>%
  dplyr::arrange(txstart, feature_id)%>%
  dplyr::distinct(gene_id, feature_id, .keep_all = TRUE)

ugly_genes_in_region <- genes_in_region%>%
  dplyr::left_join(snpeff_df, ., by = "gene_id") %>%
  dplyr::distinct(marker, CHROM, POS, log10p, peak_marker, strain, impact, .keep_all = T) %>%
  dplyr::left_join(., phenotypes, by = "strain")

tidy_genes_in_region <- genes_in_region%>%
  dplyr::left_join(snpeff_df, ., by = "gene_id") %>%
  dplyr::distinct(marker, CHROM, POS, log10p, peak_marker, strain, impact, .keep_all = T) %>%
  dplyr::left_join(., phenotypes, by = "strain") %>%
  dplyr::select(MARKER = marker, CHROM, POS, REF, ALT, MAF_variant = maf_marker_b,
                GENE_NAME = gene_name, WBGeneID = gene_id, WBFeature_TYPE = feature_type,
                WBFeature_ID = feature_id.x, TRANSCRIPT_BIOTYPE = transcript_biotype, VARIANT_IMPACT = impact,
                NUCLEOTIDE_CHANGE = nt_change, AMINO_ACID_CHANGE = aa_change,
                STRAND = strand, TRANSCRIPTION_START_POS = txstart, TRANSCRIPTION_END_POS = txend,
                PEAK_MARKER = peak_marker, PEAK_MAF = peak_maf, TRAIT = trait,
                QTL_INTERVAL_START = start_pos, QTL_INTERVAL_END = end_pos,
                VARIANT_LD_WITH_PEAK_MARKER = ld_r2, VARIANT_LOG10p = log10p,
                STRAIN = strain, STRAIN_GENOTYPE = allele, Phenotype_Value)

write_tsv(tidy_genes_in_region,
          path = glue::glue("{analysis_trait}_snpeff_genes.tsv"))


