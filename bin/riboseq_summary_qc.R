#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(patchwork))

# =========
# Collate results
# =========

# full_summary produced by riboseq_qc.R

if(sample_name == sample_names[1]) {
  full_summary <- summary_df
} else {
  full_summary <- bind_rows(full_summary, summary_df)
}



full_summary2 <- full_summary %>%
  mutate(size = case_when(str_detect(name, "10cm") ~ "10cm",
                          str_detect(name, "_6well") ~ "6well",
                          str_detect(name, "24well") ~ "24well",
                          str_detect(name, "96well") ~ "96well"))


dup <- ggplot(full_summary2, aes(x = name, y = duplication, fill = size)) +
  geom_bar(stat="identity") +
  ggtitle("Duplication %") +
  ylab("%") +
  theme_classic() +
  ggeasy::easy_rotate_x_labels()

exp <- ggplot(full_summary2, aes(x = name, y = percent_expected_length, fill = size)) +
  geom_bar(stat="identity") +
  ggtitle("% Expected length", "Low % indicates over-digestion") +
  ylab("%") +
  theme_classic() +
  ggeasy::easy_rotate_x_labels() +
  ggeasy::easy_remove_legend()

use <- ggplot(full_summary2, aes(x = name, y = y, fill = size)) +
  geom_bar(stat="identity") +
  ggtitle("% reads that are useful") +
  ylab("%") +
  theme_classic() +
  ggeasy::easy_rotate_x_labels() +
  ggeasy::easy_remove_legend()

use | exp | dup

ggsave("summary_of_results.pdf")
