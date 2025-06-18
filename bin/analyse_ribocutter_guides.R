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
    mutate(name = str_split(actual_name, "\\.min")[[1]][1]) %>%
    mutate(min_length = case_when(str_detect(actual_name, "\\.min") ~ as.integer(str_split(actual_name, "\\.min")[[1]][2]),
                                  TRUE ~ 1L))
  
  df <- mutate(df, guide_number = nrow(df))
           
  return(df)
           
}

# Load data
ribocutter.ls <- as.list(strsplit(opt$guides, ",")[[1]])

ribocutter.df <- rbindlist(lapply(ribocutter.ls, get_ribocutter_table)) %>%
  dplyr::select(name, min_length, guide_number, total_library_fraction_targeted) %>%
  distinct()

guide_number <- unique(ribocutter.df$guide_number)

ribocutter.gg <- ggplot(ribocutter.df, aes(x = name, y = total_library_fraction_targeted, fill = factor(min_length))) +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("% of reads targeted") +
  ggtitle(paste0("% of reads targeted by ", guide_number," guides with ribocutter"))

ggsave("ribocutter.pdf", ribocutter.gg, dpi = 300)

# Reformat for MultiQC report plotting 
ribocutter_mqc.df <- ribocutter.df %>%
  distinct(name, min_length, total_library_fraction_targeted) %>%
  mutate(min_length = paste0(min_length, " nt")) %>%
  pivot_wider(names_from = min_length, values_from = total_library_fraction_targeted) %>%
  dplyr::rename(sample = name) 
  

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(ribocutter_mqc.df) == nrow(ribocutter_mqc.df) + 1) {

  ribocutter_mqc.df <- ribocutter_mqc.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}
  
  

fwrite(ribocutter_mqc.df, "ribocutter_mqc.tsv", sep = "\t")
