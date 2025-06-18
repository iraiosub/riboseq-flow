#!/usr/bin/env Rscript

# Author: Oscar Wilkins

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))


df <- read_csv("~/Downloads/test_multi.csv")

df2 <- df %>%
  pivot_longer(cols = contains("mean_weight_frame"))

# Violin plot with all data, grouped visually at x = 1
ggplot(df2, aes(x = name, y = value)) +
  geom_violin(fill = "lightblue", alpha = 0.5, scale = 'width') +
  
  # Add the highlighted ORF as a red dot
  geom_point(data = subset(df2, orf_id == "ENSMUST00000028377.14_289_571_1"),
             color = "red", size = 3) +
  
  ylab("Estimated fraction") +
  
  # Improve appearance
  theme_minimal() 

scn <- read_tsv("~/Downloads/scn_orfs.tsv.gz") %>%
  distinct(orf_id)

# Violin plot with all data, grouped visually at x = 1
ggplot(df2, aes(x = name, y = value)) +
  geom_violin(fill = "lightblue", alpha = 0.5, scale = 'width') +
  
  # Add the highlighted ORF as a red dot
  geom_point(data = subset(df2, orf_id == "ENSMUST00000028377.14_691_733_1"),
             color = "red", size = 3) +
  
  ylab("Estimated fraction") +
  
  # Improve appearance
  theme_minimal() 