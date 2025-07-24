#!/usr/bin/env Rscript

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

option_list <- list(
  make_option(c("-g", "--gtf"), type = "character", help = "GTF annotation file"),
  make_option(c("-o", "--org"), type = "character", help = "Organism name")
)
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

make_txdb <- function(gtf, org) {
  name <- paste0(str_split(gtf, ".gtf")[[1]][1], ".sqlite")
  if (file.exists(name)) {
    TxDb <- loadDb(name)
  } else {
    TxDb <- makeTxDbFromGFF(gtf, format = "gtf", organism = org)
    saveDb(TxDb, file = name)
  }
  return(TxDb)
}

txdb <- make_txdb(opt$gtf, opt$org)

txlengths <- transcriptLengths(
  txdb,
  with.cds_len = TRUE,
  with.utr5_len = TRUE,
  with.utr3_len = TRUE
)
txlengths.dt <- data.table(txlengths, key = c("tx_name", "gene_id"))

pc <- c("protein_coding", "IG_V_gene", "TR_V_gene", "IG_C_gene", "IG_J_gene", 
        "TR_J_gene", "TR_C_gene", "IG_D_gene", "TR_D_gene")

gtf <- import.gff2(opt$gtf)
filtered_gtf <- gtf[!(gtf$tag %in% c("cds_end_NF", "mRNA_end_NF", "cds_start_NF", "mRNA_start_NF"))]

gtf.df <- as.data.frame(filtered_gtf)
gtf.dt <- data.table(gtf.df, key = c("transcript_id", "gene_id"))
gtf.dt <- gtf.dt[txlengths.dt, on = .(transcript_id = tx_name, gene_id = gene_id)]

# Assign TSL rank
gtf.dt[, transcript_support_level := as.character(transcript_support_level)]
gtf.dt[, tsl_rank := fcase(
  transcript_support_level == "1", 1,
  transcript_support_level == "2", 2,
  is.na(transcript_support_level), 3,
  transcript_support_level == "3", 4,
  transcript_support_level == "4", 5,
  transcript_support_level == "5", 6,
  default = 99
)]

# Choose Gencode vs Ensembl annotations
if ("gene_type" %in% colnames(gtf.df) & "transcript_type" %in% colnames(gtf.df)) {
  longest.pc.dt <- gtf.dt[gene_type %in% pc & transcript_type %in% pc]
} else if ("gene_biotype" %in% colnames(gtf.df) & "transcript_biotype" %in% colnames(gtf.df)) {
  longest.pc.dt <- gtf.dt[gene_biotype %in% pc & transcript_biotype %in% pc]
} else {
  stop("Your GTF cannot be used to select a representative transcript per gene.")
}

# Prioritize by: TSL rank → CDS length → structure
longest.pc.dt <- longest.pc.dt %>%
  arrange(tsl_rank, desc(cds_len), desc(tx_len), desc(nexon), desc(utr5_len), desc(utr3_len))

# Get unique longest per gene
unique.longest.pc.dt <- longest.pc.dt[!duplicated(longest.pc.dt$gene_id)]

# Rename gene if needed
if (!"gene_name" %in% colnames(unique.longest.pc.dt) & "gene" %in% colnames(unique.longest.pc.dt)) {
  unique.longest.pc.dt <- unique.longest.pc.dt %>%
    dplyr::rename(gene_name = gene)
}

# Prepare output
tx.info.dt <- unique.longest.pc.dt %>%
  dplyr::filter(cds_len > 0) %>%
  rowwise() %>%
  mutate(cds_start = utr5_len + 1, cds_end = utr5_len + cds_len) %>%
  ungroup() %>%
  dplyr::select(gene_name, gene_id, transcript_id, transcript_support_level, tsl_rank,
                cds_start, cds_len, tx_len, cds_end) %>%
  dplyr::rename(cds_length = cds_len, tx_length = tx_len)

# Save transcript info
output_prefix <- str_split(basename(opt$gtf), ".gtf")[[1]][1]
tx_info_name <- paste0(output_prefix , ".longest_cds.transcript_info.tsv")
fwrite(tx.info.dt, tx_info_name, sep = "\t")

# Save filtered GTF with selected transcripts
filtered_annotations <- subset(gtf, transcript_id %in% tx.info.dt$transcript_id)
gtf_name <- paste0(output_prefix , ".longest_cds_transcripts.gtf")
export.gff2(filtered_annotations, gtf_name)
