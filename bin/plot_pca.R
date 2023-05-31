#!/usr/bin/env Rscript

# DESeq2 QC on various count tables produced by iraiosub/riboseq pipeline
# Author: Ira Iosub
# Usage:plot_pca.R --featurecounts x --cds y --cds_window z --transcript_info w

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("", "--featurecounts"), action = "store", type = "character", default=NA, help = "path to multi-sample featurecounts table"),
                    make_option(c("", "--cds"), action = "store", type = "character", default=NA, help = "path to multi-sample CDS coverage table"),
                    make_option(c("", "--cds_window"), action = "store", type = "character", default=NA, help = "path to multi-sample CDS window coverage table"),
                    make_option(c("", "--transcript_info"), action = "store", type = "character", default=NA, help = "transcript info table"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



# A function that takes a count table with features as rownames, and generates PCA data using rlog normalised counts
get_rlog_pca <- function(count_data) {
  
  # Create metadata
  sample.ls <- colnames(count_data)
  
  meta.df <- data.table(sample = sample.ls, condition = 1)
  # Create dds oject
  dds <- DESeqDataSetFromMatrix(countData = round(count_data), 
                                colData = meta.df, 
                                design = ~1)
  dds <- estimateSizeFactors(dds)
  # idx <- rowSums( counts(dds, normalized = TRUE) >= 2 ) >= 1
  # dds <- dds[idx,]
  results.dds <- DESeq(dds)
  
  # PCA using rlog
  rlog <- rlogTransformation(results.dds)
  
  count_data.pca <- plotPCA(rlog, intgroup = c("sample"), returnData = FALSE)
  
  return(list(plot = count_data.pca, rlog = rlog))
  
}


# =========
# RPF gene-level counts from FeatureCounts
# =========

featurecounts.df <- fread(opt$featurecounts)

featurecounts.df <- featurecounts.df %>%
      rename_with(~str_remove(., '.genome.dedup.sorted.bam')) %>%
      rename_with(~str_remove(., '.Aligned.sortedByCoord.out.bam')) 

# If there is only one sample, there is no point in running the analysis

if (ncol(featurecounts.df) < 3 ) {
  
  message("There is only one sample provided. This analysis is only valid for 2 or more samples.") # Not enough samples in counts file for PCA.
  # file.create("pca.pdf")
  
} else {
  
  tx_info.df <- fread(opt$transcript_info)
  colours <- colorRampPalette(c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2","#91D1C2B2", "#DC0000B2", "#7E6148B2", "#f47942","#C39BD3","#fbb04e","#AAB7B8"))(ncol(featurecounts.df))
  
  
  # Reformat df and PCA of top 500 genes
  featurecounts.df <- featurecounts.df %>%
    remove_rownames() %>% 
    column_to_rownames(var = "Geneid") 
  
  featurecounts.pca.gg <- get_rlog_pca(featurecounts.df)$plot +
    geom_point(aes(color = sample)) +
    ggtitle("Gene-level counts", "FeatureCounts") +
    labs(caption = "*top 500 most variable CDS") +
    theme_cowplot() +
    scale_fill_manual(values = colours) +
    scale_color_manual(values = colours) +
    geom_text(hjust=0, vjust=0.1)

  fwrite(get_rlog_pca(featurecounts.df)$rlog, "featurecounts.rlog.tsv.gz", sep = "\t")
  
  
  if (!is.na(opt$cds)) {
    
      # =========
      # P-site CDS counts from riboWaltz
      # =========
      
      cds.df <- fread(opt$cds)
      
      # Select the longest CDS transcript based on transcript info table
      cds_longest.df <- semi_join(cds.df, tx_info.df, by = c("transcript" = "transcript_id")) %>%
        remove_rownames() %>% 
        column_to_rownames(var = "transcript")  %>%
        dplyr::select(-length_cds)
      
      cds.pca.gg <- get_rlog_pca(cds_longest.df)$plot +
        geom_point(aes(color = sample)) +
        ggtitle("Gene-level CDS occupancy", "P-sites") +
        labs(caption = "*top 500 most variable CDS") +
        theme_cowplot() +
        scale_fill_manual(values = colours) +
        scale_color_manual(values = colours) +
        geom_text(hjust=0.1, vjust=0.1)

      fwrite(get_rlog_pca(cds_longest.df)$rlog, "psite_cds_coverage.rlog.tsv.gz", sep = "\t")
      
      
      # =========
      # P-site CDS window counts from riboWaltz
      # =========

      cds_window.df <- fread(opt$cds_window)
      
      # Select the longest CDS transcript based on transcript info table
      cds_window_longest.df <- semi_join(cds_window.df, tx_info.df, by = c("transcript" = "transcript_id")) %>%
        remove_rownames() %>% 
        column_to_rownames(var = "transcript")  %>%
        dplyr::select(-length_selection, -length_cds)
      
      cds_window.pca.gg <- get_rlog_pca(cds_window_longest.df)$plot +
        geom_point(aes(color = sample)) +
        ggtitle("Gene-level CDS (+15th codon to -10th codon) occupancy", "P-sites") +
        labs(caption = "*top 500 most variable CDS") +
        theme_cowplot() +
        scale_fill_manual(values = colours) +
        scale_color_manual(values = colours) +
        geom_text(hjust=0.1, vjust=0.1)

      fwrite(get_rlog_pca(cds_window_longest.df)$rlog, "psite_cds_window_coverage.rlog.tsv.gz", sep = "\t")
      
      
      # =========
      # Merge PCA plots
      # =========
      
      pca.gg <- cowplot::plot_grid(featurecounts.pca.gg, cds.pca.gg, cds_window.pca.gg, rows = 3)
      
      ggsave("pca.pdf", pca.gg, dpi = 600, height = 30, width = 12)
      
      # save longest CDS tables
      # fwrite(semi_join(cds.df, tx_info.df, by = c("transcript" = "transcript_id")), "longest_cds_coverage_psite.tsv.gz", sep = "\t")
      # fwrite(semi_join(cds_window.df, tx_info.df, by = c("transcript" = "transcript_id")), "longest_cds_window_coverage_psite.tsv.gz", sep = "\t")
      
      
    } else {
      
      ggsave("pca.pdf", featurecounts.pca.gg, dpi = 600, height = 10, width = 12)
      
    }
}
  