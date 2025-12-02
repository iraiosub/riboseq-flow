#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--summary_list"), action = "store", type = "character", default=NA, help = "list of comma separated qc summary tables"),
                    make_option(c("-l", "--read_len_list"), action = "store", type = "character", default=NA, help = "list of tab separated read length distribution"),
                    make_option(c("-u", "--useful_len_list"), action = "store", type = "character", default=NA, help = "list of tab separated useful read length distribution"),
                    make_option(c("-r", "--region_counts_list"), action = "store", type = "character", default=NA, help = "list of tab separated region counts"),
                    make_option(c("", "--start_dist_list"), action = "store", type = "character", default=NA, help = "list of tab separated useful read count density around start codon"),
                    make_option(c("-m", "--mapping_counts_list"), action = "store", type = "character", default=NA, help = "list of tab separated mapping counts"),
                    make_option(c("-f", "--frame_counts_list"), action = "store", type = "character", default=NA, help = "list of tab separated frame counts"),
                    make_option(c("", "--length_filter_list"), action = "store", type = "character", default=NA, help = "list of tab length-filtered counts"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Collate summary results
# =========

# Summary files produced by riboseq_qc.R
summary_df.ls <- as.list(strsplit(opt$summary_list, ",")[[1]])

# summary_df.ls <- list(opt$input_list)
full_summary.df <- rbindlist(lapply(summary_df.ls, fread), use.names = TRUE)

# If no UMIs, the duplication column has NAs, in that case, the plot should be skipped

if (is.na(unique(full_summary.df$duplication))) {

  dup <- ggplot() + theme_void() + ggtitle("Duplication %") + geom_text(aes(0,0,label='N/A'))

} else {

  dup <- ggplot(full_summary.df, aes(x = name, y = duplication)) + # , fill = name
    geom_bar(stat="identity") +
    ggtitle("Duplication %") +
    ylab("%") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
    # ggeasy::easy_rotate_x_labels()

}

exp <- ggplot(full_summary.df , aes(x = name, y = percent_expected_length)) + 
  geom_bar(stat="identity") +
  ggtitle( paste0("% Expected length: ", unique(full_summary.df$expected_length)), "Low % indicates over-digestion") +
  ylab("%") +
  theme_classic() +
  # ggeasy::easy_rotate_x_labels() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

use <- ggplot(full_summary.df, aes(x = name, y = y)) +
  geom_bar(stat="identity") +
  ggtitle("% Reads that are useful", "Reads mapped uniquely to longest CDS transcripts") +
  ylab("%") +
  theme_classic() +
  # ggeasy::easy_rotate_x_labels() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

use + exp + dup + 
  plot_layout(nrow = 3)

# fwrite(full_summary.df, "qc_summary.tsv.gz", sep = "\t")

# Set dimensions for plot based on the number of samples
if (length(unique(full_summary.df$name)) >= 3) {

  plot.width <- length(unique(full_summary.df$name))
} else {

  plot.width <- 6
}

ggsave("qc_summary.pdf", dpi = 300, height = 12, width = plot.width)
# ggsave("qc_summary_mqc.png", dpi = 300, height = 12, width = plot.width)



# =========
# MultiQC tsv preparation
# =========

pcoding_percentage_mqc.df <- full_summary.df %>%
  dplyr::select(name, y) %>%
  dplyr::rename(sample = name, pcoding_percentage = y) %>%
  dplyr::arrange(sample)
fwrite(pcoding_percentage_mqc.df, "pcoding_percentage_mqc.tsv", row.names = FALSE, sep = "\t")

expected_length_mqc.df <- full_summary.df %>%
  dplyr::select(name, percent_expected_length) %>%
  dplyr::rename(sample = name, expected_length_percentage = percent_expected_length) %>%
  dplyr::arrange(sample)

fwrite(expected_length_mqc.df, "expected_length_mqc.tsv", row.names = FALSE, sep = "\t")

duplication_mqc.df <- full_summary.df %>%
  dplyr::select(name, duplication) %>%
  dplyr::rename(sample = name) %>%
  dplyr::arrange(sample)
fwrite(duplication_mqc.df, "duplication_mqc.tsv", row.names = FALSE, sep = "\t")

# FASTQ Read length files produced by riboseq_qc.R (starting read length distribution)
read_length.ls <- as.list(strsplit(opt$read_len_list, ",")[[1]])
read_length.df <- rbindlist(lapply(read_length.ls , fread), use.names = TRUE, fill = TRUE)

# Strip nt from colnames to allow linegraph
read_length.df <- read_length.df %>%
  rename_with(~str_remove(., 'nt')) %>%
  dplyr::arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)
fwrite(read_length.df, "starting_length_mqc.tsv", row.names = TRUE, sep = "\t")

# Useful read length files produced by riboseq_qc.R
useful_length.ls <- as.list(strsplit(opt$useful_len_list, ",")[[1]])
useful_length.df <- rbindlist(lapply(useful_length.ls , fread), use.names = TRUE, fill = TRUE)

# Strip nt from colnames to allow linegraph
useful_length.df <- useful_length.df %>%
  rename_with(~str_remove(., 'nt')) %>%
  dplyr::arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)
fwrite(useful_length.df, "useful_length_mqc.tsv", row.names = TRUE, sep = "\t")

# Distance from start files produced by riboseq_qc.R
start_dist.ls <- as.list(strsplit(opt$start_dist_list, ",")[[1]])
start_dist.df <- rbindlist(lapply(start_dist.ls , fread), use.names = TRUE, fill = TRUE)

# Strip nt from colnames to allow linegraph
start_dist.df <- start_dist.df %>%
  rename_with(~str_remove(., 'nt')) %>%
  dplyr::arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)
fwrite(start_dist.df, "start_dist_mqc.tsv", row.names = TRUE, sep = "\t")


# Barplots data for multiqc
# Region counts
region_counts.ls <- as.list(strsplit(opt$region_counts_list, ",")[[1]])
region_counts.df <- rbindlist(lapply(region_counts.ls , fread), use.names = TRUE, fill = TRUE) %>%
  dplyr::arrange(sample)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(region_counts.df) == nrow(region_counts.df) + 1) {

  region_counts.df <- region_counts.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}


fwrite(region_counts.df, "region_mqc.tsv", row.names = FALSE, sep = "\t")

# Mapping counts
mapping_counts.ls <- as.list(strsplit(opt$mapping_counts_list, ",")[[1]])
mapping_counts.df <- rbindlist(lapply(mapping_counts.ls , fread), use.names = TRUE) %>%
  dplyr::arrange(sample)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(mapping_counts.df) == nrow(mapping_counts.df) + 1) {

  mapping_counts.df <- mapping_counts.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}

fwrite(mapping_counts.df, "mapping_mqc.tsv", row.names = FALSE, sep = "\t")

# Frame counts
frame_counts.ls <- as.list(strsplit(opt$frame_counts_list, ",")[[1]])
frame_counts.df <- rbindlist(lapply(frame_counts.ls , fread), use.names = TRUE) %>%
  dplyr::arrange(sample)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(frame_counts.df) == nrow(frame_counts.df) + 1) {

  frame_counts.df <- frame_counts.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}

fwrite(frame_counts.df, "frame_mqc.tsv", row.names = FALSE, sep = "\t")

# Length filter counts
length_filter.ls <- as.list(strsplit(opt$length_filter_list, ",")[[1]])
length_filter.df <- rbindlist(lapply(length_filter.ls, fread), use.names = TRUE) %>%
  dplyr::arrange(sample) 

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(length_filter.df) == nrow(length_filter.df) + 1) {

  length_filter.df <- length_filter.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}

fwrite(length_filter.df, "length_filter_mqc.tsv", row.names = FALSE, sep = "\t")