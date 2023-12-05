library(tidyverse)
library(Biostrings)

# Filenames:
transcript_fasta <- "~/Downloads/gencode29.longest_cds.fa"
transcript_info <- "~/Downloads/gencode29.longest_cds.tx_info.tsv.gz"
p_sites_file <- '~/Downloads/h_r1_1.psite.tsv.gz'

# Getting some useful dataframes
transcript_df <- data.frame(transcript_seq = readDNAStringSet(transcript_fasta)) %>%
  rownames_to_column('transcript_id')

info_df <- read_tsv(transcript_info) %>%
  left_join(transcript_df, by = 'transcript_id') %>%
  filter(str_length(transcript_seq) >= cds_end + cds_start)

# Read in p site data and combine with transcript sequences
p_sites_df <- read_tsv(p_sites_file) %>%
  filter(psite > cds_start & psite < cds_stop) %>%
  select(transcript_id = transcript, psite, cds_start, cds_stop, length) %>%
  inner_join(transcript_df, by = 'transcript_id')

# Extract codons
subsample <- min(c(nrow(p_sites_df), 1E6))
codon_df <- p_sites_df %>%
  sample_n(subsample) %>% # Ensure not too many reads being processed
  # Calculate RUST scores (separated by read length)
  group_by(transcript_id, length) %>% 
  mutate(total_reads_on_tx_this_read_length = n()) %>%
  mutate(expected_reads_per_codon = total_reads_on_tx_this_read_length / (cds_stop-cds_start)) %>%
  group_by(transcript_id, psite, length) %>%
  mutate(reads_this_pos_this_length = n()) %>%
  distinct(transcript_id, psite, reads_this_pos_this_length, expected_reads_per_codon, transcript_seq, length, cds_start) %>%
  ungroup() %>%
  mutate(rust_value = ifelse(reads_this_pos_this_length > expected_reads_per_codon, 1, 0)) %>%
  # Filter for dominant frame
  mutate(frame_pos = (psite - cds_start) %% 3) %>%
  ungroup() %>%
  group_by(frame_pos) %>%
  # Filter just for max frame to keep analysis simple
  mutate(n_this_frame = sum(reads_this_pos_this_length)) %>%
  ungroup() %>%
  filter(n_this_frame == max(n_this_frame))

# Width of region of interest
offset_start = -10 # codon-wise
offset_end = 10 # codon-wise


# Find nucleotide biases for each position
all_offsets <- seq(offset_start, offset_end, by = 1)*3 # nucleotide-wise
for(offset in all_offsets){
  this_df <- codon_df %>%
    group_by(length) %>%
    mutate(this_codon = str_sub(transcript_seq, psite+offset-frame_pos, psite+offset+2-frame_pos)) %>%
    group_by(this_codon, length) %>%
    mutate(n_this_codon = sum(reads_this_pos_this_length)) %>%
    ungroup() %>%
    distinct(this_codon, n_this_codon, length) %>%
    mutate(offset = offset)
  
  if(offset == all_offsets[1]){
    full_offset_df <- this_df
  } else {
    full_offset_df <- bind_rows(full_offset_df, this_df)
  }
}

# Find expected values by looking at average codon abundances in relevant coding sequences
relevant_transcripts <- unique(codon_df$transcript_id)

expected_df <- info_df %>%
  mutate(transcript_seq = str_sub(transcript_seq, cds_start, cds_end-3)) %>%
  filter(transcript_id %in% relevant_transcripts)

# Function to extract codons from a DNA sequence
extract_codons <- function(sequence) {
  codons <- gsub("(.{3})", "\\1 ", sequence, perl = TRUE)
  codons <- unlist(strsplit(codons, " "))
  codons <- codons[codons != ""]  # Remove empty elements
  return(codons)
}

# Function to count codon occurrences across all transcripts
count_codons <- function(transcripts) {
  all_codons <- unlist(lapply(transcripts, extract_codons))
  codon_table <- table(all_codons)
  return(codon_table)
}

# Generate codon usage table across all transcripts
codon_usage_table <- data.frame(count_codons(expected_df$transcript_seq)) %>%
  filter(str_length(all_codons) == 3) %>%
  filter(!all_codons %in% c('TAA', 'TGA', 'TAG')) %>%
  mutate(q = Freq / sum(Freq)) %>%
  select(this_codon = all_codons, q)

final_df <- full_offset_df %>%
  left_join(codon_usage_table) %>%
  ungroup() %>%
  group_by(length, offset) %>%
  filter(!this_codon %in% c('TAA', 'TGA', 'TAG')) %>%
  mutate(p = n_this_codon/sum(n_this_codon)) %>%
  filter(q > 0) %>%
  mutate(kld = sum(ifelse(p == 0, 0, p*log(p/q)))) %>%
  distinct(kld)

ggplot(final_df, aes(x = offset/3, y = kld, colour = factor(length))) +
  geom_line() +
  xlab('Codon position relative to P site') +
  ylab('Sequence bias') +
  ggeasy::easy_add_legend_title('Footprint length') +
  theme_classic()

