#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))

load_count_tables <- function(table_path) {
  
  data <- fread(table_path)
  data <- data %>%
    dplyr::select(Geneid, contains(".genome.dedup.sorted.bam")) %>% 
    rename_with(~str_remove(., '.genome.dedup.sorted.bam'))
  return(data)
}


option_list <- list(make_option(c("-i", "--input_list"), action = "store", type = "character", default=NA, help = "List of comma separated count tables tables"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Raw counts
# =========

# Count tables produced by e.g. feature_counts
raw_counts.ls <- as.list(strsplit(opt$input_list, ",")[[1]])
raw_counts.df.ls <- lapply(raw_counts.df.ls, load_count_tables)
raw_counts.df <- purrr::reduce(list(x,y,z), dplyr::left_join, by = 'Geneid')

fwrite(raw_counts.df, "featurecounts.tsv.gz", sep = "\t")