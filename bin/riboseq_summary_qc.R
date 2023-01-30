#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--input_list"), action = "store", type = "character", default=NA, help = "qc summary directory"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Collate results
# =========

# Summary files produced by riboseq_qc.R
summary_df.ls <- list(opt$input_list)
full_summary.df <- rbindlist(lapply(summary_df.ls, fread), use.names = TRUE)


dup <- ggplot(full_summary.df, aes(x = name, y = duplication, fill = name)) +
  geom_bar(stat="identity") +
  ggtitle("Duplication %") +
  ylab("%") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
  # ggeasy::easy_rotate_x_labels()

exp <- ggplot(full_summary.df , aes(x = name, y = percent_expected_length, fill = name)) +
  geom_bar(stat="identity") +
  ggtitle("% Expected length", "Low % indicates over-digestion") +
  ylab("%") +
  theme_classic() +
  # ggeasy::easy_rotate_x_labels() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

use <- ggplot(full_summary.df, aes(x = name, y = y, fill = name)) +
  geom_bar(stat="identity") +
  ggtitle("% Reads that are useful") +
  ylab("%") +
  theme_classic() +
  # ggeasy::easy_rotate_x_labels() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

use | exp | dup


fwrite(full_summary.df, "qc_summary.tsv.gz", sep = "\t")
ggsave("qc_summary.pdf", dpi = 300, height = 12,
  width = 8*length(unique(full_summary.df$name)))
