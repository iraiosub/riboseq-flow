#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-b", "--bam"), action = "store", type = "character", default=NA, help = "bam file"),
                    make_option(c("-t", "--transcript_info"), action = "store", type = "character", default=NA, help = "longest protein coding transcript details"),
                    make_option(c("-o", "--output_prefix"), action = "store", type = "character", default=NA, help = "prefix for output files"),
                    make_option(c("", "--after_premap"), action = "store", type = "character", default=NA, help = "after_premap length analysis csv file"),
                    make_option(c("", "--before_dedup"), action = "store", type = "character", default=NA, help = "before_dedup length analysis csv file"),
                    make_option(c("", "--after_dedup"), action = "store", type = "character", default=NA, help = "after_dedup length analysis csv file"),
                    make_option(c("", "--expected_length"), action = "store", type = "character", default=NA, help = "string specifying expected length range, in the format min:max"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# If a read is assigned to 2 overlapping transcripts belonging to different genes, pick the one mapped to transcript with longest CDS
keep_unique_reads <- function(bam_info) {
  
  multi.ls <- bam_info$qname[duplicated(bam_info$qname)]
  
  # If we dont separate the df, slice_max is very slow
  multi.df <- bam_info %>%
    dplyr::filter(qname %in% multi.ls)
  
  original_uniq.df <- anti_join(bam_info, multi.df, by = "qname")
  
  multi2uniq.df <- multi.df %>%
    group_by(qname) %>%
    dplyr::slice_max(cds_length, with_ties = FALSE) %>%
    ungroup()
  
  uniq.df <- rbind(original_uniq.df, multi2uniq.df)
  stopifnot(length(unique((bam_info$qname))) == nrow(uniq.df))
  
  return(uniq.df)
  
}

# A function that extracts relevant positional information from transcriptomic bam for representative protein-coding transcripts
get_info_from_bam <- function(bam, info) {

  info.df <- fread(info)

  # Add sample name, keep only pcoding transcripts mapping reads
  bam.df <- data.frame(scanBam(bam)) %>%
      mutate(rl = str_length(seq)) %>%
      dplyr::rename(transcript_id = rname) %>%
      inner_join(info.df) %>%
      dplyr::filter(strand == "+") # only keep sense reads
  
  bam.df <- keep_unique_reads(bam.df) %>%
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
                  

# Get sample name
# actual_name <- str_remove(basename(opt$bam), "transcriptome.dedup.sorted.bam")
actual_name <- opt$output_prefix
message(paste0("Analysing ", actual_name))


# Define minimum and max read length based on the expected length parameter
min_length <- str_split(opt$expected_length, ":")[[1]][1]
max_length <- str_split(opt$expected_length, ":")[[1]][2]

# =========
# Frame, distance from start and end codon analyses using transcriptomic coord
# =========

riboseq_info <- get_info_from_bam(opt$bam, opt$transcript_info)

p1 <- ggplot(riboseq_info$frame, aes(x = read_length, y = n, fill=factor(frame))) +
    geom_bar(stat="identity", position="dodge") +
    theme_classic() +
    scale_fill_discrete(name = "Frame") +
    # ggeasy::easy_add_legend_title("Frame") +
    xlim(NA, 45)

p2 <- ggplot(riboseq_info$start_dist %>% filter(rl > 18 & rl < 45), 
         aes(x = distance_from_start, y = rl, fill = n)) +
    geom_raster() +
    xlim(-50, 50) +
    scale_fill_gradient(low = "white", high="black") +
    theme_classic() +
    ylab("Read length") +
    xlab("Distance of start of read from start codon") +
    ggtitle("Reads near start codon")

p3 <- ggplot(riboseq_info$end_dist %>% filter(rl > 18 & rl < 45), 
         aes(x = distance_from_end, y = rl, fill = n)) +
    geom_raster() +
    xlim(-80, 20) +
    scale_fill_gradient(low = "white", high="black") +
    theme_classic() +
    ylab("Read length") +
    xlab("Distance of start of read from stop codon") +
    ggtitle("Reads near stop codon")


# =========
# Length analysis
# =========

original_fq <- read_csv(opt$before_dedup) %>%
  # dplyr::rename(length2 = length) %>%    # these steps were done because rGrGrG hadnt been removed in that run
  # mutate(length = length2 - 3) %>%
  # dplyr::select(-length2) %>%
  dplyr::select(length, original_n = n)

before_dedup <- read_csv(opt$before_dedup) %>%
  dplyr::select(length, before_dedup_bam = n_bam)

# Plot length distribution in the adaptor trimmed and quality and length filtered fastq
length_plot <- ggplot(original_fq, aes(x = length, y = original_n)) +
  geom_bar(stat="identity") +
  xlim(0,70) +
  # ggpubr::theme_pubr() +
  theme_classic() +
  ylab("N reads in fastq") +
  ggtitle("Read length distribution")  

# Reformat length dataframe to export for MultiQC, read length distribution of starting reads
fq_length_mqc.df <- original_fq %>%
  mutate(sample = actual_name) %>%
  dplyr::rename(number_of_reads = original_n) %>%
  dplyr::select(sample, length, number_of_reads) %>%
  mutate(length = paste0(length, "nt")) %>%
  pivot_wider(names_from = length, values_from = number_of_reads)

fwrite(fq_length_mqc.df, paste0(actual_name, "_fq_length_mqc.tsv"), sep = "\t", row.names = FALSE)

# Read length distribution of useful reads
useful_length_mqc.df <- riboseq_info$bam %>%
  mutate(sample = actual_name) %>%
  group_by(sample, rl) %>%
  summarise(number_of_reads = n()) %>%
  dplyr::rename(length = rl) %>%
  mutate(length = paste0(length, "nt")) %>%
  pivot_wider(names_from = length, values_from = number_of_reads)

fwrite(useful_length_mqc.df, paste0(actual_name, "_useful_length_mqc.tsv"), sep = "\t", row.names = FALSE)


# Regional distribution of useful reads
# Assign CDS or UTR based on relative position of the starts of useful reads to CDS start and ends
region_counts.df <- riboseq_info$bam %>%
  dplyr::select(qname, pos, rl, cds_start, cds_end) %>%
  mutate(region = case_when(pos < cds_start ~ "5UTR",
                            pos > cds_end ~ "3UTR",
                            TRUE ~ "CDS"),
         sample = actual_name) %>%
  group_by(sample, region) %>%
  summarise(region_counts = n())

# Now reformat for multiqc
region_counts_mqc.df <- region_counts.df %>%
  pivot_wider(names_from = region, values_from = region_counts)

fwrite(region_counts_mqc.df, paste0(actual_name, "_region_counts_mqc.tsv"), sep = "\t", row.names = FALSE)


# Start dist, reformat length dataframe to export for MultiQC
start_dist_mqc.df <- riboseq_info$start_dist %>%
  mutate(sample = actual_name) %>%
  dplyr::filter(distance_from_start <=50 & distance_from_start >= -50) %>%
  dplyr::filter(rl >= min_length & rl <= max_length) %>%
  group_by(sample, distance_from_start) %>%
  summarise(number_of_reads = sum(n)) %>%
  dplyr::select(sample, distance_from_start, number_of_reads) %>%
  mutate(distance_from_start = paste0(distance_from_start, "nt")) %>%
  pivot_wider(names_from = distance_from_start, values_from = number_of_reads)

fwrite(start_dist_mqc.df, paste0(actual_name, "_start_dist_mqc.tsv"), sep = "\t", row.names = FALSE)

# Frames for multiqc, summarised for the expected length range
frame_mqc.df <- riboseq_info$frame %>%
  mutate(sample = actual_name) %>%
  dplyr::filter(read_length >= min_length & read_length <= max_length) %>%
  group_by(sample, frame) %>%
  summarise(number_of_reads = sum(n)) %>%
  dplyr::select(sample, frame, number_of_reads) %>%
  mutate(frame = paste0("frame ", frame)) %>%
  pivot_wider(names_from = frame, values_from = number_of_reads)

fwrite(frame_mqc.df, paste0(actual_name, "_frame_mqc.tsv"), sep = "\t", row.names = FALSE)


# =========
# Premapping
# =========

name_colours <- c("after_premap_n" = "#f8766d", "original_n" = "#00BFC4", "before_dedup_bam" = "#7851a9")

if (!is.na(opt$after_premap)) {
# if (basename(opt$after_premap) != "optional.txt") {

  after_premap <- read_csv(opt$after_premap) %>%
  dplyr::select(length, after_premap_n = n)

  rRNA_df <- inner_join(original_fq, after_premap) %>%
  mutate(perc_rRNA = 100-100*after_premap_n/original_n)

  rRNA_plot <- ggplot(rRNA_df,
                    aes(x = length, y = perc_rRNA)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_classic() +
  # ggpubr::theme_pubr() +
  xlim(19,60) +
  ylab("% rRNA") +
  ggtitle("rRNA %")

  mapping_df <- inner_join(original_fq, after_premap) %>%
    inner_join(before_dedup) %>%
    mutate(total = 100*before_dedup_bam/original_n,
            ignoring_rRNA = 100*before_dedup_bam/after_premap_n) %>%
    dplyr::select(length, total, ignoring_rRNA) %>%
    pivot_longer(-length)

  mapping_df2 <- inner_join(original_fq, after_premap) %>%
    inner_join(before_dedup) %>%
    pivot_longer(-length)

  mapping_plot <- ggplot(mapping_df2 %>% filter(length < 50) %>%
    # drop_na(value) %>%
    dplyr::filter(name != "original_n"), aes(x = length, y = value, fill = name)) +
    geom_bar(stat="identity", position="dodge", na.rm = T) +
    scale_fill_manual(values = name_colours) +
    theme_classic() +
    ylab("Number of reads") +
    ggtitle("Premapping vs mapping") +
    theme(legend.position = "none", legend.title=element_blank()) +
    # ggeasy::easy_remove_legend() +
    scale_y_log10() +
    expand_limits(y=1)
    # annotate(geom = "rect", xmin = 26, xmax = 31, ymin = -Inf, ymax = Inf, alpha = .1, color = "black", linetype = "dashed", fill = "black")

  premapping_plot <- ggplot(mapping_df2 %>% filter(length < 50) %>%
    dplyr::filter(name != "before_dedup_bam"), aes(length, y = value, fill = name)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values = name_colours) +
    theme_classic() +
    ylab("Number of reads") +
    ggtitle("Original vs premapping") +
    # ggeasy::easy_remove_legend() +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title=element_blank(), legend.text = element_text(size=6)) +
    scale_y_log10()+
    expand_limits(y=1)
    # annotate(geom = "rect", xmin = 26, xmax = 31, ymin = -Inf, ymax = Inf, alpha = .1, color = "black", linetype = "dashed", fill = "black")

  # ggplot(mapping_df %>% filter(length<50), aes(x = length, y = value, fill = name)) +
  #   geom_bar(stat="identity", position="dodge")

} else {

  # rRNA_plot <- mapping_plot <- premapping_plot <- ggplot() + theme_void()
  rRNA_plot  <- ggplot() + theme_void() + ggtitle("rRNA %") + geom_text(aes(0,0,label='N/A'))
  mapping_plot  <- ggplot() + theme_void() + ggtitle("Premapping vs mapping") + geom_text(aes(0,0,label='N/A'))
  premapping_plot <- ggplot() + theme_void() + ggtitle("Original vs mapping") + geom_text(aes(0,0,label='N/A'))

}

# =========
# Duplication
# =========


# If UMIs were used, calculate the proprtion of duplicated reads within the expected RPF length range
# if (!is.null(opt$after_dedup)) {
# if(basename(opt$after_dedup) != "optional.txt") {
if(!is.na(opt$after_dedup)) {

  after_dedup <- read_csv(opt$after_dedup) %>%
    dplyr::select(length, after_dedup_bam = n)

  duplicate_df <- inner_join(before_dedup, after_dedup) %>%
    mutate(perc_duplicates = 100*(before_dedup_bam-after_dedup_bam)/before_dedup_bam)

  # Filter based on expected RPF length 
  duplication_summary <- duplicate_df %>%
    dplyr::filter(length >= min_length & length <= max_length)
  
  # This is the percentage reported in the ribo-seq summary
  duplication_perc <- 100- 100*sum(duplication_summary$after_dedup_bam) / sum(duplication_summary$before_dedup_bam) 

  duplication_plot <- ggplot(duplicate_df %>% filter(length >= min_length & length <= max_length), 
        aes(x = length, y = perc_duplicates)) +
  geom_bar(stat="identity", position="dodge") +
  ylab("% Duplicates") +
  theme_classic() +
  ggtitle("Duplication")

} else {

  duplication_plot <- ggplot() + theme_void() + ggtitle("Duplication") + geom_text(aes(0,0,label='N/A'))
  duplication_perc <- NA
}


# =========
# Useful reads
# =========

# Useful reads are defined as reads mapping to protein coding/number of preprocessed reads
useful_read_perc <- 100*nrow(riboseq_info$bam) / sum(original_fq$original_n)
useful_df <- data.frame(x = "", y = useful_read_perc,
                        useful_read_n = nrow(riboseq_info$bam))

# Proprtion of uniquely mapped reads to reppresentative pcoding transcripts ie useful reads that are of expected read length
tx_map_expected_length_perc <- 100*nrow(riboseq_info$bam %>% filter(rl >= min_length & rl <= max_length)) / nrow(riboseq_info$bam)
tx_map_expected_length_count <- nrow(riboseq_info$bam %>% filter(rl >= min_length & rl <= max_length))

summary_df <- useful_df %>%
  mutate(name = actual_name,
          duplication = duplication_perc,
          expected_length = opt$expected_length,
          percent_expected_length = tx_map_expected_length_perc,
          expected_length_n = tx_map_expected_length_count)


# Customise sub-title based on whether UMIs were used or not
# if(basename(opt$after_dedup) != "optional.txt") {
if(!is.na(opt$after_dedup)) {
  useful_plot_subtitle <- "UMI-deduplicated reads mapped\nuniquely to longest CDS transcripts"
} else {

  useful_plot_subtitle <- "Reads mapped\nuniquely to longest CDS transcripts"    
}

useful_plot <- ggplot(useful_df, aes(x= x, y = y)) +
  geom_bar(stat="identity") +
  theme_classic() +
  ggtitle("% Useful reads", useful_plot_subtitle) +
  xlab("") +
  ylab("Proportion (%)")

# Pull together plots
fig <- ((p1 + ggtitle(paste("Sample:", actual_name)) | useful_plot) +  plot_layout(widths = c(3, 1))) / (p2|p3) / 
  (length_plot | rRNA_plot) /
  (premapping_plot | mapping_plot | duplication_plot)

# Save results
ggsave(paste0(actual_name, ".qc_results.pdf"), plot = fig, dpi = 300, height = 18, width = 12)
fwrite(summary_df, paste0(actual_name, ".qc_results.tsv.gz"), sep = "\t", row.names = FALSE)

