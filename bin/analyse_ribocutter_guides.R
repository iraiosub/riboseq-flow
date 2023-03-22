#!/usr/bin/env Rscript

# Script for ribocuttter QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 22.03.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-g", "--guides"), action = "store", type = "character", default=NA, help = "list of comma separated ribocutter guides csv files"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


get_ribocutter_table <- function(ribocutter_csv) {
  
  actual_name <- str_remove(basename(ribocutter_csv), ".csv")
  
  df <- read_csv(ribocutter_csv) %>%
    mutate(name = str_split(actual_name, ".min")[[1]][1]) %>%
    mutate(min_length = case_when(str_detect(actual_name, ".min") ~ as.integer(str_split(actual_name, ".min")[[1]][2]),
                                  TRUE ~ 1L))
           
  return(df)
           
}

# Load data
ribocutter.ls <- as.list(strsplit(opt$guides, ",")[[1]])

ribocutter.df <- rbindlist(lapply(ribocutter.ls, get_ribocutter_table)) %>%
  dplyr::select(name, min_length, total_targeted) %>%
  distinct()

ribocutter.gg <- ggplot(ribocutter.df, aes(x = name, y = total_targeted, fill = factor(min_length))) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("% of reads targeted") +
  ggtitle("% of reads targeted by 50 guides with ribocutter",
          "Higher is better")

ribocutter.gg


ggsave("ribocutter.pdf", ribocutter.gg, dpi = 300)