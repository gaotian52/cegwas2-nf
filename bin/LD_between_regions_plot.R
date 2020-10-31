.libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths() ))
library(genetics) 
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
gm <- read.table(file = args[1], header = T)
processed_mapping <- read.delim(args[2], stringsAsFactors=FALSE)
TRAIT <- args[3]

snp_df <- processed_mapping %>% na.omit()

ld_snps <- dplyr::filter(gm, CHROM %in% snp_df$CHROM, POS %in% snp_df$POS)


if ( nrow(ld_snps) > 1 ) {
    
    ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS,
                                         sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
    
    sn <- list()
    
    for (i in 1:nrow(ld_snps)) {
        sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
                                                        gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
    }
    
    test <- data.frame(sn)
    colnames(test) <- (ld_snps$snp_id)
    ldcalc <- t(genetics::LD(test)[[4]])^2
    diag(ldcalc) <- 1
    
    write.table(ldcalc, paste0(TRAIT, "_LD_between_QTL_regions.tsv"), quote=F, row.names = T, col.names = NA, sep="\t")
    
    
    LDs <- tbl_df(data.frame(ldcalc) %>% dplyr::add_rownames(var = "SNP1")) %>%
      tidyr::gather(SNP2, r2, -SNP1) %>% dplyr::arrange(SNP1) %>%
      tidyr::separate(SNP1, sep = "_", into = c("CHROM1","POS1"), remove = F) %>%
      dplyr::arrange(CHROM1, as.numeric(POS1)) 
    

    
    LD_plot <- ggplot2::ggplot(LDs) +
      ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), y = factor(SNP2, levels = unique(SNP1), ordered = T)) +
      ggplot2::geom_tile(ggplot2::aes(fill = as.numeric(r2))) +
      ggplot2::geom_text(ggplot2::aes(label = signif(r2, 3)),  size = 12*5/14) +
      ggplot2::theme(text = ggplot2::element_text(size = 12, color = "black"),
                     axis.text.x = ggplot2::element_text(size = 8, color = "black"),
                     axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     legend.text = element_text(size=12,  color = "black"),
                     legend.title =   element_text(size=12,  color = "black"),
                     legend.position = "left"
      ) +
      scale_x_discrete(labels = function(x) { gsub("_", ":", x) }, expand = c(0, 0)) + 
      scale_y_discrete(position = "right", labels = function(x) { gsub("_", ":", x) }, expand = c(0, 0)) + 
      scale_fill_continuous(high = "red", low = "white", na.value = "white") +
      labs(fill=expression(R^{2}))
    
    ggsave(LD_plot, filename = paste(TRAIT,"_LD_plot.png",sep=""), units = "mm",height = 100, width = 100)
    
    
}
