#!/usr/bin/env Rscript

# Script for QC analyses. Originally written by Oscar Wilkins.
# Modified by Ira Iosub 24.01.2023

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--summary_list"), action = "store", type = "character", default=NA, help = "list of comma separated qc summary tables"),
                    make_option(c("-l", "--read_len_list"), action = "store", type = "character", default=NA, help = "list of tab separated read length distribution"),
                    make_option(c("-u", "--useful_len_list"), action = "store", type = "character", default=NA, help = "list of tab separated useful read length distribution"),
                    make_option(c("-r", "--region_counts_list"), action = "store", type = "character", default=NA, help = "list of tab separated region counts"),
                    make_option(c("", "--start_dist_list"), action = "store", type = "character", default=NA, help = "list of tab separated useful read count density around start codon"),
                    make_option(c("-m", "--mapping_counts_list"), action = "store", type = "character", default=NA, help = "list of tab separated mapping counts"),
                    make_option(c("-f", "--frame_counts_list"), action = "store", type = "character", default=NA, help = "list of tab separated frame counts"),
                    make_option(c("", "--length_filter_list"), action = "store", type = "character", default=NA, help = "list of tab length-filtered counts"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Collate summary results
# =========

# Summary files produced by riboseq_qc.R
summary_df.ls <- as.list(strsplit(opt$summary_list, ",")[[1]])

# summary_df.ls <- list(opt$input_list)
full_summary.df <- rbindlist(lapply(summary_df.ls, fread), use.names = TRUE)

# If no UMIs, the duplication column has NAs, in that case, the plot should be skipped

has_duplication <- any(!is.na(full_summary.df$duplication))

if (has_duplication) {

  dup <- ggplot(full_summary.df, aes(x = name, y = duplication)) + # , fill = name
    geom_bar(stat="identity") +
    ggtitle("Duplication %") +
    ylab("%") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
    # ggeasy::easy_rotate_x_labels()

}

exp <- ggplot(full_summary.df , aes(x = name, y = percent_expected_length)) + 
  geom_bar(stat="identity") +
  ggtitle( paste0("% Expected length: ", unique(full_summary.df$expected_length)), "Low % indicates over-digestion") +
  ylab("%") +
  theme_classic() +
  # ggeasy::easy_rotate_x_labels() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

use <- ggplot(full_summary.df, aes(x = name, y = y)) +
  geom_bar(stat="identity") +
  ggtitle("% Reads that are useful", "Reads mapped uniquely to longest CDS transcripts") +
  ylab("%") +
  theme_classic() +
  # ggeasy::easy_rotate_x_labels() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none")

fwrite(full_summary.df, "qc_summary.tsv.gz", sep = "\t", row.names = FALSE)



# =========
# MultiQC tsv preparation
# =========

pcoding_percentage_mqc.df <- full_summary.df %>%
  dplyr::select(name, y) %>%
  dplyr::rename(sample = name, pcoding_percentage = y) %>%
  dplyr::arrange(sample)
fwrite(pcoding_percentage_mqc.df, "pcoding_percentage_mqc.tsv", row.names = FALSE, sep = "\t")

expected_length_mqc.df <- full_summary.df %>%
  dplyr::select(name, percent_expected_length) %>%
  dplyr::rename(sample = name, expected_length_percentage = percent_expected_length) %>%
  dplyr::arrange(sample)

outside_expected_length.df <- full_summary.df %>%
  dplyr::filter(most_abundant_length_in_expected_range == FALSE) %>%
  dplyr::arrange(name)

if (nrow(outside_expected_length.df) > 0) {
  length_warnings <- paste0(htmltools::htmlEscape(outside_expected_length.df$name), ": ",
                           outside_expected_length.df$most_abundant_length, " nt (expected ",
                           outside_expected_length.df$expected_length, " nt)")
  writeLines(c("# description: |",
               "#   Percentage of useful reads within the expected read length range. The higher the percentage, the better.",
               paste0('#   <div class="alert alert-warning"><strong>Warning:</strong> Most abundant useful read length outside the expected range: ',
                      paste(length_warnings, collapse = "; "), ".</div>")), "expected_length_mqc.tsv")
}

fwrite(expected_length_mqc.df, "expected_length_mqc.tsv", row.names = FALSE, sep = "\t",
       append = nrow(outside_expected_length.df) > 0, col.names = TRUE)

duplication_mqc.df <- full_summary.df %>%
  dplyr::select(name, duplication) %>%
  dplyr::rename(sample = name) %>%
  dplyr::arrange(sample)
fwrite(duplication_mqc.df, "duplication_mqc.tsv", row.names = FALSE, sep = "\t")

# FASTQ Read length files produced by riboseq_qc.R (starting read length distribution)
read_length.ls <- as.list(strsplit(opt$read_len_list, ",")[[1]])
read_length.df <- rbindlist(lapply(read_length.ls , fread), use.names = TRUE, fill = TRUE)

# Strip nt from colnames to allow linegraph
read_length.df <- read_length.df %>%
  rename_with(~str_remove(., 'nt')) %>%
  dplyr::arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)
fwrite(read_length.df, "starting_length_mqc.tsv", row.names = TRUE, sep = "\t")

# Useful read length files produced by riboseq_qc.R
useful_length.ls <- as.list(strsplit(opt$useful_len_list, ",")[[1]])
useful_length.df <- rbindlist(lapply(useful_length.ls , fread), use.names = TRUE, fill = TRUE)

# Strip nt from colnames to allow linegraph
useful_length.df <- useful_length.df %>%
  rename_with(~str_remove(., 'nt')) %>%
  dplyr::arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)
fwrite(useful_length.df, "useful_length_mqc.tsv", row.names = TRUE, sep = "\t")

# Distance from start files produced by riboseq_qc.R
start_dist.ls <- as.list(strsplit(opt$start_dist_list, ",")[[1]])
start_dist.df <- rbindlist(lapply(start_dist.ls , fread), use.names = TRUE, fill = TRUE)

# Strip nt from colnames to allow linegraph
start_dist.df <- start_dist.df %>%
  rename_with(~str_remove(., 'nt')) %>%
  dplyr::arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  replace(is.na(.), 0)
fwrite(start_dist.df, "start_dist_mqc.tsv", row.names = TRUE, sep = "\t")


# Barplots data for multiqc
# Region counts
region_counts.ls <- as.list(strsplit(opt$region_counts_list, ",")[[1]])
region_counts.df <- rbindlist(lapply(region_counts.ls , fread), use.names = TRUE, fill = TRUE) %>%
  dplyr::arrange(sample) %>%
  replace(is.na(.), 0)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(region_counts.df) == nrow(region_counts.df) + 1) {

  region_counts.df <- region_counts.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}


fwrite(region_counts.df, "region_mqc.tsv", row.names = FALSE, sep = "\t")

# Mapping counts
mapping_counts.ls <- as.list(strsplit(opt$mapping_counts_list, ",")[[1]])
mapping_counts.df <- rbindlist(lapply(mapping_counts.ls , fread), use.names = TRUE) %>%
  dplyr::arrange(sample)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(mapping_counts.df) == nrow(mapping_counts.df) + 1) {

  mapping_counts.df <- mapping_counts.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}

fwrite(mapping_counts.df, "mapping_mqc.tsv", row.names = FALSE, sep = "\t")

# Frame counts
frame_counts.ls <- as.list(strsplit(opt$frame_counts_list, ",")[[1]])
frame_counts.df <- rbindlist(lapply(frame_counts.ls , fread), use.names = TRUE) %>%
  dplyr::arrange(sample)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(frame_counts.df) == nrow(frame_counts.df) + 1) {

  frame_counts.df <- frame_counts.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}

fwrite(frame_counts.df, "frame_mqc.tsv", row.names = FALSE, sep = "\t")

# Frame counts for the most abundant useful read length
most_abundant_frame_mqc.df <- full_summary.df %>%
  dplyr::select(name, most_abundant_frame_0, most_abundant_frame_1, most_abundant_frame_2) %>%
  dplyr::rename(sample = name, `frame 0` = most_abundant_frame_0, `frame 1` = most_abundant_frame_1, `frame 2` = most_abundant_frame_2) %>%
  dplyr::arrange(sample)

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(most_abundant_frame_mqc.df) == nrow(most_abundant_frame_mqc.df) + 1) {

  most_abundant_frame_mqc.df <- most_abundant_frame_mqc.df %>%
    mutate(x = 0)
}

fwrite(most_abundant_frame_mqc.df, "most_abundant_frame_mqc.tsv", row.names = FALSE, sep = "\t")

# Frame counts for the most periodic supported expected read length
most_periodic_frame_mqc.df <- full_summary.df %>%
  dplyr::filter(!is.na(most_periodic_length)) %>%
  dplyr::select(name, most_periodic_frame_0, most_periodic_frame_1, most_periodic_frame_2) %>%
  dplyr::rename(sample = name, `frame 0` = most_periodic_frame_0, `frame 1` = most_periodic_frame_1, `frame 2` = most_periodic_frame_2) %>%
  dplyr::arrange(sample)

writeLines("# plot_type: bargraph", "most_periodic_frame_mqc.tsv")
fwrite(most_periodic_frame_mqc.df, "most_periodic_frame_mqc.tsv", row.names = FALSE, sep = "\t", append = TRUE, col.names = TRUE)

# Selected read lengths, together for all samples
selected_lengths_mqc.df <- full_summary.df %>%
  dplyr::select(name, expected_length, most_abundant_length, most_periodic_length) %>%
  dplyr::rename(sample = name, `Expected length range (nt)` = expected_length,
                `Most abundant length (nt)` = most_abundant_length, `Most periodic expected length (nt)` = most_periodic_length) %>%
  dplyr::arrange(sample)

# Include table configuration here because MultiQC infers and replaces pconfig for TSV tables
writeLines(c("# plot_type: table",
             '# pconfig: {id: selected_lengths_table, title: "Selected Read Lengths", scale: false, format: "{:,.0f}"}'),
           "selected_lengths_mqc.tsv")
fwrite(selected_lengths_mqc.df, "selected_lengths_mqc.tsv", row.names = FALSE, sep = "\t", append = TRUE, col.names = TRUE)

# Dominant frame percentage for each of the three read length categories
frame_metrics.df <- bind_rows(
  frame_counts.df %>%
    left_join(full_summary.df %>% dplyr::select(sample = name, expected_length), by = "sample") %>%
    transmute(sample, category = "Expected length range", length = paste0(expected_length, " nt"),
              frame_0 = `frame 0`, frame_1 = `frame 1`, frame_2 = `frame 2`),
  full_summary.df %>%
    transmute(sample = name, category = "Most abundant length", length = paste0(most_abundant_length, " nt"),
              frame_0 = most_abundant_frame_0, frame_1 = most_abundant_frame_1, frame_2 = most_abundant_frame_2),
  full_summary.df %>%
    transmute(sample = name, category = "Most periodic expected length", length = paste0(most_periodic_length, " nt"),
              frame_0 = most_periodic_frame_0, frame_1 = most_periodic_frame_1, frame_2 = most_periodic_frame_2)
) %>%
  mutate(total = frame_0 + frame_1 + frame_2,
          weakest_frame = pmin(frame_0, frame_1, frame_2),
          dominant_percentage = if_else(total > 0, 100*pmax(frame_0, frame_1, frame_2) / total, NA_real_))

dominant_frame_mqc.df <- frame_metrics.df %>%
  dplyr::select(sample, category, dominant_percentage) %>%
  pivot_wider(names_from = category, values_from = dominant_percentage) %>%
  dplyr::arrange(sample)

writeLines("# plot_type: bargraph", "dominant_frame_mqc.tsv")
fwrite(dominant_frame_mqc.df, "dominant_frame_mqc.tsv", row.names = FALSE, sep = "\t", append = TRUE, col.names = TRUE)

# Dominant / least abundant frame ratios
frame_ratio_mqc.df <- full_summary.df %>%
  dplyr::select(name, expected_frame_ratio, most_abundant_frame_ratio, most_periodic_frame_ratio) %>%
  dplyr::rename(sample = name, `Expected length range` = expected_frame_ratio, `Most abundant length` = most_abundant_frame_ratio,
                `Most periodic expected length` = most_periodic_frame_ratio) %>%
  dplyr::arrange(sample)

# Keep unavailable ratios from changing MultiQC's inferred plot type
frame_ratio_header <- "# plot_type: bargraph"
zero_frame.df <- frame_metrics.df %>%
  dplyr::filter(total > 0, weakest_frame == 0) %>%
  dplyr::arrange(sample, category)

if (nrow(zero_frame.df) > 0) {
  frame_ratio_notes <- paste0(htmltools::htmlEscape(zero_frame.df$sample), ": ", zero_frame.df$category,
                             " (", zero_frame.df$length, "; frame counts ", zero_frame.df$frame_0, " / ",
                             zero_frame.df$frame_1, " / ", zero_frame.df$frame_2, "; n=", zero_frame.df$total, ")")
  frame_ratio_header <- c(frame_ratio_header, "# description: |",
                          "#   Ratio of the largest to the smallest 5' end frame count for each read length category. Ratios are unavailable when the smallest count is zero or no eligible length exists.",
                          paste0('#   <p class="text-muted"><small><strong>Missing ratio bars:</strong> no finite ratio because the weakest frame has zero reads.<br>',
                                 paste(frame_ratio_notes, collapse = "<br>"), "</small></p>"))
}

writeLines(frame_ratio_header, "frame_ratio_mqc.tsv")
fwrite(frame_ratio_mqc.df, "frame_ratio_mqc.tsv", row.names = FALSE, sep = "\t", append = TRUE, col.names = TRUE)

# Length filter counts
length_filter.ls <- as.list(strsplit(opt$length_filter_list, ",")[[1]])
length_filter.df <- rbindlist(lapply(length_filter.ls, fread), use.names = TRUE) %>%
  dplyr::arrange(sample) 

# If dim of df are n x n, multiqc overrides plot_type and plots a heatmap
if (ncol(length_filter.df) == nrow(length_filter.df) + 1) {

  length_filter.df <- length_filter.df %>%
    mutate(x = 0) # workaround to add dummy col so multiqc plots as barplot and not heatmap
}

fwrite(length_filter.df, "length_filter_mqc.tsv", row.names = FALSE, sep = "\t")

# Combined QC report, including the same frame metrics as MultiQC
frame_colours <- c("Expected length range" = "#7cb5ec", "Most abundant length" = "#434348", "Most periodic expected length" = "#90ed7d")
frame_labels <- c("Expected length range" = "Expected length\nrange", "Most abundant length" = "Most abundant\nlength", "Most periodic expected length" = "Most periodic\nexpected length")

dominant_frame_plot <- ggplot(frame_metrics.df, aes(x = sample, y = dominant_percentage, fill = category)) +
  geom_bar(stat = "identity", position = "dodge", na.rm = TRUE) +
  geom_text(aes(label = if_else(category == "Expected length range", "", length)),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = frame_colours, labels = frame_labels, name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1))) +
  ggtitle("Dominant frame percentage", "Labels show selected lengths; the expected range pools lengths.") +
  xlab("Sample") + ylab("%") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")

frame_ratio_plot.df <- frame_ratio_mqc.df %>%
  pivot_longer(-sample, names_to = "category", values_to = "ratio") %>%
  mutate(ratio = as.numeric(ratio))

frame_ratio_plot <- ggplot(frame_ratio_plot.df, aes(x = sample, y = ratio, fill = category)) +
  geom_bar(stat = "identity", position = "dodge", na.rm = TRUE) +
  geom_text(aes(y = 0, label = if_else(is.na(ratio), "NA", "")),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
  scale_fill_manual(values = frame_colours, labels = frame_labels, name = NULL) +
  ggtitle("Dominant / least abundant frame ratio") +
  labs(caption = "NA: zero reads in the weakest frame, or no eligible read length.") +
  xlab("Sample") + ylab("Ratio") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")

summary_plot <- use + exp
if (has_duplication) {
  summary_plot <- summary_plot + dup
}
summary_plot <- summary_plot + dominant_frame_plot + frame_ratio_plot +
  plot_layout(nrow = if (has_duplication) 5 else 4)

# Allow room for the frame legends and increase width with the number of samples
plot.width <- max(8, length(unique(full_summary.df$name)))
ggsave("qc_summary.pdf", plot = summary_plot, dpi = 300, height = if (has_duplication) 18 else 15, width = plot.width)
# ggsave("qc_summary_mqc.png", dpi = 300, height = 12, width = plot.width)
