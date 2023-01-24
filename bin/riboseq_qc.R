#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(patchwork))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-b", "--bam"), action = "store", type = "character", default=NA, help = "bam file"),
                    make_option(c("-t", "--transcript_info"), action = "store", type = "character", default=NA, help = "longest protein coding transcript details"),
                    make_option(c("", "--after_premap"), action = "store", type = "character", default=NA, help = "after_premap length analysis csv file"),
                    make_option(c("", "--before_dedup"), action = "store", type = "character", default=NA, help = "before_dedup length analysis csv file"),
                    make_option(c("", "--after_dedup"), action = "store", type = "character", default=NA, help = "after_dedup length analysis csv file"),)
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


get_info_from_bam <- function(bam, info) {

    # add sample name
    bam.df <- data.frame(scanBam(bam)) %>%
        mutate(rl = str_length(seq)) %>%
        dplyr::rename(transcript_id = rname) %>%
        inner_join(info) %>%
        mutate(distance_from_start = pos - cds_start,
          distance_from_end = pos - cds_end) %>%
        mutate(frame = distance_from_start %% 3)

    frame.df <- bam.df %>% 
      group_by(frame, rl) %>%
      mutate(n = n()) %>%
      select(frame, read_length = rl, n) %>%
      unique() %>%
      ungroup()

    start_dist.df <- bam.df %>%
      group_by(rl, distance_from_start) %>%
      mutate(n = n()) %>%
      select(rl, distance_from_start, n) %>%
      unique() %>%
      ungroup()

    end_dist.df <- bam.df %>%
      group_by(rl, distance_from_end) %>%
      mutate(n = n()) %>%
      select(rl, distance_from_end, n) %>%
      unique() %>%
      ungroup()

    return(list(bam = bam.df, frame = frame.df, start_dist = start_dist.df, end_dist = end_dist.df))

}

# =========
# Frame, distance from start and end codon analyses
# =========

riboseq_info <- get_info_from_bam(opt$bam, opt$transcript_info)

p1 <- ggplot(riboseq_info$frame, aes(x = read_length, y = n, fill=factor(frame))) +
    geom_bar(stat="identity", position="dodge") +
    theme_classic() +
    ggeasy::easy_add_legend_title("Frame") +
    xlim(NA, 40)

p2 <- ggplot(riboseq_info$start_dist %>% filter(rl > 18 & rl < 40), 
         aes(x = distance_from_start, y = rl, fill = n)) +
    geom_raster() +
    xlim(-50, 50) +
    scale_fill_gradient(low = "white", high="black") +
    theme_classic() +
    ylab("Read length") +
    xlab("Distance of start of read from start codon") +
    ggtitle("Reads near start codon")

p3 <- ggplot(riboseq_info$end_dist %>% filter(rl > 18 & rl < 40), 
         aes(x = distance_from_end, y = rl, fill = n)) +
    geom_raster() +
    xlim(-80, 20) +
    scale_fill_gradient(low = "white", high="black") +
    theme_classic() +
    ylab("Read length") +
    xlab("Distance of start of read from stop codon") +
    ggtitle("Reads near stop codon")



# # Make a list of bam files from --bam 
# bams <- Sys.glob("C:/Users/ogw/Downloads/cristina_tso/*.bam")

# sample_data <- read_csv("Copy of PM22312.csv")
# info <- read_tsv("C:/Users/ogw/Downloads/longest_proteincoding_transcript_hs_details.txt")

# sample_names <- str_replace(word(word(bams, -1, -1, sep="/"), 1, sep="\\."), "Aligned", "")

# for(sample_name in sample_names){
#   actual_name <- sample_data$`Sample Name`[which(sample_data$`Sample limsid`==word(sample_name, 1, sep="_"))]
#   print(actual_name)  
#   df <- data.frame(scanBam(paste0("C:/Users/ogw/Downloads/cristina_tso/",
#                                   sample_name,
#                                   "Aligned.sortedByCoord.out.bam")))



# =========
# Length analysis
# =========

  original_fq <- read_csv(opt$before_dedup) %>%
    dplyr::rename(length2 = length) %>%
    mutate(length = length2 - 3) %>%
    dplyr::select(-length2) %>%
    dplyr::select(length, original_n = n)
  
  before_dedup <- read_csv(opt$before_dedup) %>%
    dplyr::select(length, before_dedup_bam = n_bam)
  
  after_premap <- read_csv(opt$after_premap) %>%
    dplyr::select(length, after_premap_n = n)
  
  after_dedup <- read_csv(opt$after_dedup) %>%
    dplyr::select(length, after_dedup_bam = n)

  rRNA_df <- inner_join(original_fq, after_premap) %>%
    mutate(perc_rRNA = 100-100*after_premap_n/original_n)

  mapping_df <- inner_join(original_fq, after_premap) %>%
    inner_join(before_dedup) %>%
    mutate(total = 100*before_dedup_bam/original_n,
           ignoring_rRNA = 100*before_dedup_bam/after_premap_n) %>%
    dplyr::select(length, total, ignoring_rRNA) %>%
    pivot_longer(-length)

  mapping_df2 <- inner_join(original_fq, after_premap) %>%
    inner_join(before_dedup) %>%
    pivot_longer(-length)
  
  duplicate_df <- inner_join(before_dedup, after_dedup) %>%
    mutate(perc_duplicates = 100*(before_dedup_bam-after_dedup_bam)/before_dedup_bam)


  useful_read_perc <- 100*nrow(riboseq_info$bam) / sum(original_fq$original_n)
  useful_df <- data.frame(x = "", y = useful_read_perc)
  
  length_plot <- ggplot(original_fq, aes(x = length, y = original_n)) +
    geom_bar(stat="identity") +
    xlim(0,70) +
    ggpubr::theme_pubr() +
    ylab("N reads in original fastq") +
    ggtitle("Read length distribution")  

  
  rRNA_plot <- ggplot(rRNA_df,
                      aes(x = length, y = perc_rRNA)) +
    geom_bar(stat="identity", position = "dodge") +
    ggpubr::theme_pubr() +
    xlim(19,60) +
    ylab("% rRNA") +
    ggtitle("rRNA %")
  
  
  ggplot(mapping_df %>% filter(length<50), aes(x = length, y = value, fill = name)) +
    geom_bar(stat="identity", position="dodge")

  
  mapping_plot <- ggplot(mapping_df2 %>% filter(length<50) %>%
           filter(name != "original_n"), aes(x = length, y = value, fill = name)) +
    geom_bar(stat="identity", position="dodge") +
    theme_classic() +
    ylab("Number of reads") +
    ggtitle("Premapping vs mapping") +
    ggeasy::easy_remove_legend() +
    scale_y_log10()
  
  premapping_plot <- ggplot(mapping_df2 %>% filter(length<50) %>%
           filter(name != "before_dedup_bam"), aes(x = length, y = value, fill = name)) +
    geom_bar(stat="identity", position="dodge") +
    theme_classic() +
    ylab("Number of reads") +
    ggtitle("Original vs premapping") +
    ggeasy::easy_remove_legend() +
    scale_y_log10()
  
  duplication_plot <- ggplot(duplicate_df %>% filter(length >= 26 & length <= 32), 
         aes(x = length, y = perc_duplicates)) +
    geom_bar(stat="identity", position="dodge") +
    ylab("% Duplicates") +
    theme_classic() +
    ggtitle("Duplication")

  useful_plot <- ggplot(useful_df, aes(x= x, y = y)) +
    geom_bar(stat="identity") +
    theme_classic() +
    ggtitle("% useful reads") +
    xlab("") +
    ylab("% useful")
  
  fig <- ((p1 + ggtitle(paste("Sample:", actual_name)) | useful_plot) +  plot_layout(widths = c(3, 1))) / (p2|p3) / 
    (length_plot | rRNA_plot) /
    (premapping_plot | mapping_plot | duplication_plot)
  fig
  
  
  ggsave(paste0("figures/", actual_name, "_results2.pdf"), plot = fig)
  
  
  duplication_summary <- duplicate_df %>%
    filter(length >= 26 & length <= 31)
  duplication_perc <- 100- 100*sum(duplication_summary$after_dedup_bam) / sum(duplication_summary$before_dedup_bam) 
  
  tx_map_summary <- 100*nrow(df3 %>% filter(rl >= 26 & rl <= 31)) / nrow(df3)
  

  
  summary_df <- useful_df %>%
    mutate(name = actual_name,
           duplication = duplication_perc,
           percent_expected_length = tx_map_summary)
  
  if(sample_name == sample_names[1]) {
    full_summary <- summary_df
  } else {
    full_summary <- bind_rows(full_summary, summary_df)
  }



# =========
# Collate results
# =========

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

ggsave("figures/summary_of_results.pdf")
