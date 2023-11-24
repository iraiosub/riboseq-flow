#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-g", "--gtf"), action = "store", type = "character", default=NA, help = "GTF annotation file"),
                    make_option(c("-o", "--org"), action = "store", type = "character", default=NA, help = "string specifying organism name"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

make_txdb <- function(gtf, org) {
  
  # Get name of db
  name <- paste0(str_split(gtf, ".gtf")[[1]][1], ".sqlite")
  
  if (file.exists(name)) {
    
    TxDb <- loadDb(name)
    
  } else {
    
    TxDb <- makeTxDbFromGFF(gtf, format="gtf",
                            organism = org) # chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene
    
    saveDb(TxDb, file=name)
    # TxDb <- loadDb(name)
  }
  
  # TxDb <- keepStandardChromosomes(TxDb, pruning.mode="coarse")
  return(TxDb)
}



txdb <- make_txdb(opt$gtf, opt$org)

txlengths <- transcriptLengths(txdb, with.cds_len = TRUE,
                               with.utr5_len = TRUE,
                               with.utr3_len = TRUE)

txlengths.dt <- data.table(txlengths, key = c("tx_name", "gene_id"))
pc = c("protein_coding", "IG_V_gene", "TR_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "IG_D_gene", "TR_D_gene")


gtf <- import.gff2(opt$gtf)
filtered_gtf <- gtf[!(gtf$tag %in% c("cds_end_NF", "mRNA_end_NF", "cds_start_NF", "mRNA_start_NF"))]
# filtered_gtf <- filtered_gtf[filtered_gtf$transcript_support_level %in% 1:2]


gtf.df <- as.data.frame(filtered_gtf)
gtf.dt <- data.table(gtf.df, key = c("transcript_id", "gene_id"))
gtf.dt <- gtf.dt[txlengths.dt]

if ("gene_type" %in% colnames(gtf.df) & "transcript_type" %in% colnames(gtf.df)) {
  # Gencode
  longest.pc.dt <- gtf.dt[gene_type %in% pc & transcript_type %in% pc, longest := max(cds_len), by = gene_id] # select out where both are protein coding as sometimes a processed transcript is the longest
  longest.pc.dt <- longest.pc.dt[gene_type %in% pc & transcript_type %in% pc & cds_len == longest] # selects longest

} else if ("gene_biotype" %in% colnames(gtf.df) & "transcript_biotype" %in% colnames(gtf.df)) {
  # Ensembl
  longest.pc.dt <- gtf.dt[gene_biotype %in% pc & transcript_biotype %in% pc, longest := max(cds_len), by = gene_id] # select out where both are protein coding as sometimes a processed transcript is the longest
  longest.pc.dt <- longest.pc.dt[gene_biotype %in% pc & transcript_biotype %in% pc & cds_len == longest] # selects longest
} else {

  stop("Your GTF cannot be used to select a representative transcript per gene. Either use your own transcript info file, or use Gencode or Ensembl annotations")
}

# Hierarchy: CDS > tx_len > n_exon > UTR3 > UTR3
longest.pc.dt <- longest.pc.dt %>% 
  arrange(desc(cds_len), desc(longest), desc(nexon), desc(utr5_len), desc(utr3_len))

unique.longest.pc.dt <- longest.pc.dt[!duplicated(longest.pc.dt$gene_id), ] 

tx.info.dt <- unique.longest.pc.dt %>%
  dplyr::filter(cds_len > 0) %>%
  rowwise() %>%
  mutate(cds_start = utr5_len + 1, cds_end = utr5_len + cds_len) %>%
  ungroup() %>%
  dplyr::select(gene_name, gene_id, transcript_id, cds_start, cds_len, tx_len, cds_end) %>%
  dplyr::rename(cds_length = cds_len, tx_length = tx_len)


output_name <- paste0(str_split(basename(opt$gtf), ".gtf")[[1]][1], ".longest_cds.transcript_info.tsv")

fwrite(tx.info.dt, output_name, sep = "\t")
