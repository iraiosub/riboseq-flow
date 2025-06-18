#!/usr/bin/env Rscript

# Author: Oscar Wilkins

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggeasy))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(data.table))

option_list <- list(make_option(c("-i", "--input"), action = "store", type = "character", default=NA, help = "riboloco-lite output"),
                    make_option(c("-t", "--transcript_info"), action = "store", type = "character", default=NA, help = "transcript info"),
                    make_option(c("-l", "--lengths"), action = "store", type = "character", default="27:31", help = "lengths of interest"),
                    make_option(c("-c", "--min_footprints"), action = "store", type = "integer", default=10, help = "min footprints in orf"),
                    make_option(c("-p", "--min_unique_footprint_positions"), action = "store", type = "integer", default=3, help = "min unique footprint positions"),
                    make_option(c("-g", "--gene_names"), action = "store", type = "character", default=NULL,help = "Comma-separated list of gene name patterns (e.g., 'Scn,Grin1')"),
                    make_option(c("-o", "--output"), action = "store", type = "character", help = "output prefix"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load input
riboloco_lite_output <- opt$input
info <- readr::read_tsv(opt$transcript_info)

lengths_of_interest <- eval(parse(text = opt$lengths))
min_footprints_in_orf <- opt$min_footprints
min_unique_footprint_positions <- opt$min_unique_footprint_positions
prefix <- opt$output

downsample_for_heatmaps <- 20

# Genes to label
gene_pattern <- gsub(",", "|", opt$gene_names)
gene_info <- info %>%
  filter(str_detect(gene_name, gene_pattern))


#### RUN ####
df <- readr::read_csv(riboloco_lite_output) %>%
  mutate(within_orf = A_site_estimate >= orf_start & A_site_estimate <= orf_stop)

# Get ORF types counts
df_unique_orf_identifier <- df %>%
  dplyr::filter(within_orf) %>%
  ungroup() %>%
  mutate(orf_id = paste(transcript_id, orf_start, orf_stop, orf_frame, sep = "_"))


df_unique_orf_identifier_footprint_types <- df_unique_orf_identifier %>%
  group_by(orf_id, footprint_type) %>%
  summarise(footprint_type_n = n()) %>%
  group_by(orf_id) %>%
  dplyr::filter(sum(footprint_type_n) > min_footprints_in_orf)

write_csv(df_unique_orf_identifier_footprint_types,
       paste0(prefix,".footprint_types_per_orf.csv.gz"))


# Assign ORF labels based on overlap with annotated ORFs
orf_labels <- df %>%
  dplyr::select(transcript_id, orf_start, orf_stop, annotated) %>%
  distinct() %>%
  group_by(transcript_id) %>%
  mutate(annotated_start = max(ifelse(annotated == 1, orf_start, -1)),
         annotated_stop = max(ifelse(annotated == 1, orf_stop, -1))) %>%
  mutate(orf_label = case_when(annotated == 1 ~ 'Annotated CDS',
                               annotated_start == -1 ~ 'No annotated start',
                               orf_start > annotated_stop ~ 'Downstream ORF',
                               orf_stop < annotated_start ~ 'Upstream ORF',
                               orf_start < annotated_start & orf_stop > annotated_start ~ 'Upstream overlapping ORF',
                               orf_stop > annotated_stop & orf_start < annotated_stop ~ 'Downstream overlapping ORF',
                               T ~ 'Overlapping ORF'))

# Find fractions of each footprint in annotated regions
annotated_fractions <- df %>%
  dplyr::filter(annotated == 1) %>%
  dplyr::filter(within_orf) %>%
  group_by(footprint_type) %>%
  mutate(n = n()) %>%
  distinct(footprint_type, n) %>%
  mutate(length = as.numeric(word(footprint_type, 1, sep="_")),
         frame = as.numeric(word(footprint_type, 2, sep="_")),
         mismatch = str_detect(footprint_type, "MM")) %>%
  filter(length %in% lengths_of_interest) %>%
  ungroup() %>%
  mutate(frac = n/sum(n))


write_csv(annotated_fractions, paste0(prefix, ".annotated_fractions.csv.gz"))


# Produce heatmap of footprint type densities near start and stop codons
df2 <- df %>%
  group_by(transcript_id) %>%
  filter(annotated == 1) %>%
  mutate(cds_start = max(ifelse(annotated == 1, orf_start, -1000))) %>%
  mutate(cds_stop  = max(ifelse(annotated == 1, orf_stop, -1000))) %>%
  mutate(Position_relative_to_CDS_start = A_site_estimate - cds_start,
         Position_relative_to_CDS_stop = A_site_estimate - cds_stop) %>%
  ungroup() %>%
  dplyr::select(Position_relative_to_CDS_start, Position_relative_to_CDS_stop, footprint_type) %>%
  pivot_longer(cols = contains('rel')) %>%
  mutate(downsampled_pos = round(value / downsample_for_heatmaps)) %>%
  dplyr::select(-value) %>%
  group_by(name, footprint_type, downsampled_pos) %>%
  mutate(n = n()) %>%
  distinct() %>%
  ungroup() %>%
  group_by(name, downsampled_pos) %>%
  mutate(actual_frac = n/sum(n)) %>%
  inner_join(annotated_fractions %>% dplyr::select(footprint_type, expected_frac = frac)) %>%
  mutate(log2fc = log2(actual_frac/expected_frac)) %>%
  mutate(n_footprints_total = sum(n))


heatmap_data <- df2 %>% filter(abs(downsampled_pos) < 200/downsample_for_heatmaps)

if (nrow(heatmap_data) > 0) {
  p1 <- ggplot(heatmap_data, 
               aes(x = downsample_for_heatmaps*downsampled_pos, y = footprint_type, fill = log2fc)) +
    geom_tile() +
    facet_wrap(~name) +
    scale_fill_viridis_c() +
    theme_classic() +
    xlab('Relative nucleotide position') +
    ggeasy::easy_add_legend_title('Log2(enrichment)') +
    ggtitle('Enrichment of each footprint type near start and stop codons',
            'Expected values are based on annotated CDS regions') +
    ggeasy::easy_remove_legend()

  p2 <- ggplot(heatmap_data %>% distinct(name, downsampled_pos, n_footprints_total),
               aes(x = downsample_for_heatmaps*downsampled_pos, y = n_footprints_total)) +
    geom_area() +
    geom_point() +
    facet_wrap(~name) +
    xlab('Relative nucleotide position') +
    ggtitle('Total footprint density') +
    theme_classic() +
    ylab('Total_footprints')

  combined = p1 / p2 + plot_layout(heights = c(3, 1))
  ggsave(paste0(prefix, ".footprint_types_heatmap.pdf"), combined, dpi = 300, height = 25, width = 20, units = 'cm')
} else {
  message("No data available for heatmap plot. Skipping ggsave.")
}



# p3 <- ggplot(annotated_fractions, aes(x = 0, y = footprint_type, fill = frac)) +
#   geom_tile() +
#   scale_fill_continuous(low = 'white', high = 'red') +
#   theme_classic() +
#   ggeasy::easy_remove_legend()
# combined2 = p3/plot_spacer() + plot_layout(heights = c(3, 1))
# (combined | combined2) + plot_layout(widths = c(8, 1))


### Identify ORFs with high amount of unexpected footprints ####
df3 <- df %>%
  inner_join(annotated_fractions) %>%
  dplyr::filter(within_orf) %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(n_distinct(A_site_estimate) > min_unique_footprint_positions) %>%
  ungroup() %>%
  group_by(transcript_id, footprint_type, orf_start) %>%
  mutate(n_type_in_orf = n()) %>%
  distinct(n_type_in_orf, annotated) %>%
  inner_join(annotated_fractions) %>%
  ungroup() %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(sum(n_type_in_orf) > min_footprints_in_orf) %>%
  mutate(frac_in_orf = n_type_in_orf / sum(n_type_in_orf)) %>% 
  mutate(kl_div = sum(frac_in_orf * log(frac_in_orf / frac))) %>%
  mutate(n_in_orf = sum(n_type_in_orf)) %>%
  distinct(kl_div, annotated, n_in_orf) %>%
  left_join(orf_labels) %>%
  ungroup() %>%
  group_by(orf_label) %>%
  mutate(r = rank(kl_div)) %>%
  left_join(gene_info) %>%
  mutate(label = gene_name)



if (nrow(df3) > 0 && sum(!is.na(df3$kl_div)) > 0) {
    p4 <- ggplot(df3, aes(y = kl_div, x = orf_label)) +
      geom_violin(scale = "width") +
      theme_classic() +
      ggeasy::easy_rotate_x_labels(side = 'right') +
      ylab("KullbackLeibler Divergence") +
      ggtitle('"Weirdness of distribution of footprint types" score',
              'Higher score = more weird = more likely out-of-frame') +
      xlab("Type of ORF")

    ggsave(paste0(prefix,'.footprint_heatmap.pdf'), p4, height = 25, width = 20, units = 'cm')


    p5 <- ggplot(df3 %>% ungroup() %>% arrange(desc(label)), aes(x = r, y = kl_div, colour = label)) +
      geom_point() +
      theme_classic() +
      ylab("KullbackLeibler Divergence") +
      xlab('Rank') +
      facet_wrap(~orf_label, scales = 'free_x')  +
      ggrepel::geom_label_repel(aes(label = label, alpha = 0.5), box.padding = 2) +
      ggeasy::easy_remove_legend() +
      ggtitle('"Weirdness of distribution of footprint types" score',
            'Higher score = more weird = more likely out-of-frame')

    ggsave(paste0(prefix,'.kl_div_rank_plot.pdf'),
        p5, height = 25, width = 25, units = 'cm')
  } else {
    message("Skipping p4 and p5 plots: no valid data.")
  }





# combined_2 =(p4 | p5) + plot_layout(widths = c(1, 3))
# ggplot(df3 %>% ungroup() %>% arrange(desc(label)) %>% filter(orf_label == 'Upstream overlapping ORF'),
#        aes(x = r, y = kl_div, colour = label)) +
#   geom_point() +
#   theme_classic() +
#   ylab("KullbackLeibler Divergence") +
#   xlab('Rank') +
#   facet_wrap(~orf_label, scales = 'free_x')  +
#   ggrepel::geom_label_repel(aes(label = label, alpha = 0.5), box.padding = 2) +
#   ggeasy::easy_remove_legend() +
#   ggtitle('"Weirdness of distribution of footprint types" score',
#           'Higher score = more weird = more likely out-of-frame') 

write_csv(df3, paste0(prefix, ".kld_div.csv.gz"))

# Frame shift analyses
df3_0 <- df %>%
  inner_join(annotated_fractions) %>%
  dplyr::filter(within_orf) %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(n_distinct(A_site_estimate) > min_unique_footprint_positions) %>%
  ungroup() %>%
  group_by(transcript_id, footprint_type, orf_start) %>%
  mutate(n_type_in_orf = n()) %>%
  distinct(n_type_in_orf, annotated) %>%
  inner_join(annotated_fractions) %>%
  ungroup() %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(sum(n_type_in_orf) > min_footprints_in_orf) %>%
  mutate(frac_in_orf = n_type_in_orf / sum(n_type_in_orf)) %>%
  mutate(kl_div = sum(frac_in_orf * log(frac_in_orf / frac))) %>%
  mutate(n_in_orf = sum(n_type_in_orf)) %>%
  distinct(kl_div, annotated, n_in_orf) %>%
  left_join(orf_labels) %>%
  ungroup() %>%
  group_by(orf_label) %>%
  mutate(r = rank(kl_div)) %>%
  left_join(gene_info) %>%
  mutate(label = gene_name) %>%
  mutate(offset = 0)


df3_1 <- df %>%
  inner_join(annotated_fractions) %>%
  dplyr::filter(within_orf) %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(n_distinct(A_site_estimate) > min_unique_footprint_positions) %>%
  ungroup() %>%
  group_by(transcript_id, footprint_type, orf_start) %>%
  mutate(n_type_in_orf = n()) %>%
  distinct(n_type_in_orf, annotated) %>%
  inner_join(annotated_fractions) %>%
  ungroup() %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(sum(n_type_in_orf) > min_footprints_in_orf) %>%
  mutate(new_footprint_type = paste0(length, "_", (frame+1)%%3, "_", ifelse(mismatch, "MM", "m"))) %>%
  dplyr::select(transcript_id, footprint_type = new_footprint_type, orf_start, n_type_in_orf, annotated, n) %>%
  left_join(annotated_fractions, by = 'footprint_type') %>%
  mutate(frac_in_orf = n_type_in_orf / sum(n_type_in_orf)) %>%
  mutate(kl_div = sum(frac_in_orf * log(frac_in_orf / frac))) %>%
  mutate(n_in_orf = sum(n_type_in_orf)) %>%
  distinct(kl_div, annotated, n_in_orf) %>%
  left_join(orf_labels) %>%
  ungroup() %>%
  group_by(orf_label) %>%
  mutate(r = rank(kl_div)) %>%
  left_join(gene_info) %>%
  mutate(label = gene_name) %>%
  mutate(offset = 1)


df3_2 <- df %>%
  inner_join(annotated_fractions) %>%
  dplyr::filter(within_orf) %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(n_distinct(A_site_estimate) > min_unique_footprint_positions) %>%
  ungroup() %>%
  group_by(transcript_id, footprint_type, orf_start) %>%
  mutate(n_type_in_orf = n()) %>%
  distinct(n_type_in_orf, annotated) %>%
  inner_join(annotated_fractions) %>%
  ungroup() %>%
  group_by(transcript_id, orf_start) %>%
  dplyr::filter(sum(n_type_in_orf) > min_footprints_in_orf) %>%
  mutate(new_footprint_type = paste0(length, "_", (frame+2)%%3, "_", ifelse(mismatch, "MM", "m"))) %>%
  dplyr::select(transcript_id, footprint_type = new_footprint_type, orf_start, n_type_in_orf, annotated, n) %>%
  left_join(annotated_fractions, by = 'footprint_type') %>%
  mutate(frac_in_orf = n_type_in_orf / sum(n_type_in_orf)) %>% 
  mutate(kl_div = sum(frac_in_orf * log(frac_in_orf / frac))) %>%
  mutate(n_in_orf = sum(n_type_in_orf)) %>%
  distinct(kl_div, annotated, n_in_orf) %>%
  left_join(orf_labels) %>%
  ungroup() %>%
  group_by(orf_label) %>%
  mutate(r = rank(kl_div)) %>%
  left_join(gene_info) %>%
  mutate(label = gene_name) %>%
  mutate(offset = 2)



if (nrow(df3_0) > 0 && nrow(df3_1) > 0 && nrow(df3_2) > 0) {

  df4 <- bind_rows(df3_0, df3_1, df3_2) %>%
    dplyr::select(transcript_id, orf_start, kl_div, offset, orf_label) %>%
    pivot_wider(names_from = offset, values_from = kl_div) %>%
    mutate(decrease_in_kl_with_1 = `0` - `1`) %>%
    mutate(decrease_in_kl_with_2 = `0` - `2`) %>%
    left_join(gene_info) %>%
    dplyr::select(transcript_id, orf_start, orf_label, decrease_in_kl_with_1, decrease_in_kl_with_2) %>%
    pivot_longer(cols = contains('decrease')) %>%
    left_join(gene_info) %>%
    group_by(orf_label, name) %>%
    mutate(r = rank(value)) %>%
    mutate(label = ifelse(orf_label == 'Upstream overlapping ORF', gene_name, NA))

  decrease.gg <- ggplot(df4, aes(x = r, y = value, colour = label)) +
    geom_point() +
    facet_grid(cols = vars(orf_label), rows = vars(name), scales = 'free_x') +
    ggrepel::geom_label_repel(aes(label = label), alpha = 0.5, box.padding = 1)

  ggsave(paste0(prefix, '.kl_div_decrease.pdf'),
         decrease.gg, height = 20, width = 40, units = 'cm')

  write_csv(df4, paste0(prefix, ".kld_div_decrease.csv.gz"))

} else {
  message("Skipping KL-div decrease plots: one or more shifted dataframes are empty.")
}
