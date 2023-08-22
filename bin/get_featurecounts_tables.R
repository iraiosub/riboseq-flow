#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))


# A function that loads quantification files depending on whether they had been produced from deduplicated bam files or not
load_count_tables <- function(table_path) {
  
  data <- fread(table_path)
  
  data <- data %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length)
  
  if (TRUE %in% str_detect(colnames(data), ".genome.dedup.sorted.bam")) {
    data <- data %>%
      dplyr::select(Geneid, contains(".genome.dedup.sorted.bam")) %>% 
      rename_with(~str_remove(., '.genome.dedup.sorted.bam'))
    
  } else if (!(TRUE %in% str_detect(colnames(data), ".genome.dedup.sorted.bam")) & TRUE %in% str_detect(colnames(data), ".Aligned.sortedByCoord.out.bam")) {
    
    data <- data %>%
      dplyr::select(Geneid, contains(".Aligned.sortedByCoord.out.bam")) %>% 
      rename_with(~str_remove(., '.Aligned.sortedByCoord.out.bam')) 
    
  }
  
  return(data)
}


option_list <- list(make_option(c("-i", "--input_list"), action = "store", type = "character", default=NA, help = "List of comma separated count tables"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Raw counts
# =========

# Count tables produced by e.g. feature_counts
raw_counts.ls <- as.list(strsplit(opt$input_list, ",")[[1]])

raw_counts.df.ls <- lapply(raw_counts.ls, load_count_tables)
raw_counts.df <- purrr::reduce(raw_counts.df.ls, dplyr::left_join, by = 'Geneid')

# Reorder columns in alphabetical order, but keep Gene as the 1st column
raw_counts.df <- raw_counts.df %>%
  dplyr::select(Geneid, sort(colnames(.))) %>%
  dplyr::rename(Gene = Geneid)

fwrite(raw_counts.df, "merged.featureCounts.tsv.gz", sep = "\t")