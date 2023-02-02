#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--input_list"), action = "store", type = "character", default=NA, help = "list of comma separated count tables tables"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Collate results
# =========

# Count tables priduced by e.g. featurecounts
counts.ls <- as.list(strsplit(opt$input_list, ",")[[1]])

counts.df.ls <- lapply(counts.df.ls, fread)





fwrite(counts.df, "qc_summary.tsv.gz", sep = "\t")
