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
  idx <- rowSums( counts(dds, normalized = TRUE) >= 2 ) >= 3
  dds <- dds[idx,]
  results.dds <- DESeq(dds)
  
  # PCA using rlog
  rlog <- rlogTransformation(results.dds)
  rlog.df <- data.frame(assay(rlog))
  
  # Plot
  count_data.pca <- plotPCA(rlog, intgroup = c("sample"), returnData = FALSE)

  # Data for the plot
  count_data.pca_data <- plotPCA(rlog, intgroup = c("sample"), returnData = TRUE) %>%
    dplyr::select(sample, PC1, PC2) %>%
    column_to_rownames(var = "sample")
  
  return(list(plot = count_data.pca, data = count_data.pca_data, rlog = rlog.df))
  
}


# =========
# RPF gene-level counts from FeatureCounts
# =========

# Create empty file if plot cannot be made 
file.create("pca_mqc.png")

featurecounts.df <- fread(opt$featurecounts)

featurecounts.df <- featurecounts.df %>%
      rename_with(~str_remove(., '.genome.dedup.sorted.bam')) %>%
      rename_with(~str_remove(., '.Aligned.sortedByCoord.out.bam')) 

# If there is only one sample, there is no point in running the analysis

if (ncol(featurecounts.df) < 4 ) {
  
  message("There aren't enough samples provided. This analysis is only valid for 3 or more samples.") # Not enough samples in counts file for PCA.
  # file.create("pca.pdf")
  
} else {
  
  tx_info.df <- fread(opt$transcript_info)
  
  # Reformat df and PCA of top 500 genes
  featurecounts.df <- featurecounts.df %>%
    remove_rownames() %>% 
    column_to_rownames(var = "Gene") 
  
  featurecounts.pca.gg <- get_rlog_pca(featurecounts.df)$plot +
    geom_point(color = "606060") +
    ggtitle("Gene-level counts", "featureCounts (rlog-normalised counts)") +
    labs(caption = "*top 500 most variable genes") +
    theme_cowplot() +
    # ggrepel::geom_text_repel(aes(label = sample))
    geom_text(aes(label = sample), size = 3, hjust = -0.1, vjust = 0.8) +
    theme(legend.position = "none") +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    coord_cartesian(clip = "off")

  fwrite(get_rlog_pca(featurecounts.df)$rlog, "featurecounts.rlog.tsv.gz", sep = "\t")
  fwrite(get_rlog_pca(featurecounts.df)$data, "featurecounts_pca_mqc.tsv", row.names = TRUE, sep = "\t")
  
  
  if (!is.na(opt$cds)) {
    
      # =========
      # P-site CDS counts from riboWaltz
      # =========
      
      cds.df <- fread(opt$cds)

      if (ncol(cds.df) < 5) {
        message("There aren't enough samples provided. This analysis is only valid for 3 or more samples.")
        # Create empty plot to specify the analysis is not available
        cds.pca.gg <- ggplot() + theme_void() + ggtitle("CDS occupancy", "P-sites (rlog-normalised counts)") + geom_text(aes(0,0,label='N/A'))
      } else {

        # Select the longest CDS transcript based on transcript info table
        cds_longest.df <- semi_join(cds.df, tx_info.df, by = c("transcript" = "transcript_id")) %>%
          remove_rownames() %>% 
          column_to_rownames(var = "transcript")  %>%
          dplyr::select(-length_cds)
        
        cds.pca.gg <- get_rlog_pca(cds_longest.df)$plot +
          geom_point(color = "606060") +
          ggtitle("CDS occupancy", "P-sites (rlog-normalised counts)") +
          labs(caption = "*top 500 most variable CDS") +
          theme_cowplot() +
          # ggrepel::geom_text_repel(aes(label = sample))
          geom_text(aes(label = sample), size = 3, hjust = -0.1, vjust = 0.8)+
          theme(legend.position = "none") +
          theme(plot.margin = unit(c(1,1,1,1), "cm")) +
          coord_cartesian(clip = "off")

        fwrite(get_rlog_pca(cds_longest.df)$rlog, "psite_cds_coverage.rlog.tsv.gz", sep = "\t")
        fwrite(get_rlog_pca(cds_longest.df)$data, "psite_pca_mqc.tsv", row.names = TRUE, sep = "\t")

      }
      
      # =========
      # P-site CDS window counts from riboWaltz
      # =========

      cds_window.df <- fread(opt$cds_window)

      if (ncol(cds_window.df) < 5) {
        message("There aren't enough samples provided. This analysis is only valid for 3 or more samples.")

        # Create empty plot to specify the analysis is not available
        cds_window.pca.gg <- ggplot() + theme_void() + ggtitle("CDS occupancy", "P-sites (rlog-normalised counts)") + geom_text(aes(0,0,label='N/A'))
      
      } else {
      
        # Select the longest CDS transcript based on transcript info table
        cds_window_longest.df <- semi_join(cds_window.df, tx_info.df, by = c("transcript" = "transcript_id")) %>%
          remove_rownames() %>% 
          column_to_rownames(var = "transcript")  %>%
          dplyr::select(-length_selection, -length_cds)
        
        cds_window.pca.gg <- get_rlog_pca(cds_window_longest.df)$plot +
          geom_point(color = "606060") +
          ggtitle("CDS (+15th codon to -10th codon) occupancy", "P-sites (rlog-normalised counts)") +
          labs(caption = "*top 500 most variable CDS") +
          theme_cowplot() +
        # ggrepel::geom_text_repel(aes(label = sample))
          geom_text(aes(label = sample), size = 3, hjust = -0.1, vjust = 0.8) +
          theme(legend.position = "none") +
          theme(plot.margin = unit(c(1,1,1,1), "cm")) +
          coord_cartesian(clip = "off")

        fwrite(get_rlog_pca(cds_window_longest.df)$rlog, "psite_cds_window_coverage.rlog.tsv.gz", sep = "\t")
        fwrite(get_rlog_pca(cds_window_longest.df)$data, "psite_cds_window_pca_mqc.tsv", row.names = TRUE, sep = "\t")

      }
      
      # =========
      # Merge PCA plots
      # =========
      
      pca.gg <- cowplot::plot_grid(featurecounts.pca.gg, cds.pca.gg, cds_window.pca.gg, nrow = 3)
      ggsave("pca.pdf", pca.gg, dpi = 600, height = 30, width = 12)

      pca_mqc.gg <- cowplot::plot_grid(featurecounts.pca.gg, cds.pca.gg, cds_window.pca.gg, ncol = 3)
      ggsave("pca_mqc.png", pca_mqc.gg, dpi = 600, height = 6, width = 21)
      
      # save longest CDS tables
      # fwrite(semi_join(cds.df, tx_info.df, by = c("transcript" = "transcript_id")), "longest_cds_coverage_psite.tsv.gz", sep = "\t")
      # fwrite(semi_join(cds_window.df, tx_info.df, by = c("transcript" = "transcript_id")), "longest_cds_window_coverage_psite.tsv.gz", sep = "\t")
      
      
    } else {
      
      ggsave("pca.pdf", featurecounts.pca.gg, dpi = 600, height = 10, width = 10)
      gsave("pca_mqc.png", pca.gg, dpi = 600, height = 6, width = 7)
      
    }
}
  