#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

make_hg38_txdb <- function(gtf) {
  
  # Get name of db
  name <- paste0(str_split(gtf, ".gtf")[[1]][1], ".sqlite")
  
  if (file.exists(name)) {
    
    TxDb <- loadDb(name)
    
  } else {
    
    TxDb <- makeTxDbFromGFF(opt$gtf, format="gtf",
                            organism = "Homo sapiens") # chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene
    
    saveDb(TxDb, file=name)
    # TxDb <- loadDb(name)
  }
  
  TxDb <- keepStandardChromosomes(TxDb, pruning.mode="coarse")
  return(TxDb)
}



option_list <- list(make_option(c("-g", "--gtf"), action = "store", type = "character", default=NA, help = "GTF annotation file"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

txdb <- make_hg38_txdb(opt$gtf)

txlengths <- transcriptLengths(txdb, with.cds_len = TRUE,
                               with.utr5_len = TRUE,
                               with.utr3_len = TRUE)

txlengths.dt <- data.table(txlengths, key = c("tx_name", "gene_id"))
pc = c("protein_coding", "IG_V_gene", "TR_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "IG_D_gene", "TR_D_gene")


gtf.df <- as.data.frame(import.gff2(opt$gtf))
gtf.dt <- data.table(gtf.df, key = c("transcript_id", "gene_id"))
gtf.dt <- gtf.dt[txlengths.dt]

longest.pc.dt <- gtf.dt[gene_type %in% pc & transcript_type %in% pc, longest := max(cds_len), by = gene_id] # select out where both are protein coding as sometimes a processed transcript is the longest
longest.pc.dt <- longest.pc.dt[gene_type %in% pc & transcript_type %in% pc & cds_len == longest] # selects longest

# hierarchy
longest.pc.dt <- longest.pc.dt %>% arrange(desc(cds_len), desc(longest), desc(nexon), desc(utr3_len), desc(utr5_len))

unique.longest.pc.dt <- longest.pc.dt[ !duplicated(longest.pc.dt$gene_id), ] 

tx.info.dt <- unique.longest.pc.dt %>%
  dplyr::filter(cds_len > 0) %>%
  rowwise() %>%
  mutate(cds_start = utr5_len + 1, cds_end = utr5_len + cds_len) %>%
  ungroup() %>%
  dplyr::select(gene_name, gene_id, transcript_id, cds_start, cds_len, tx_len, cds_end) %>%
  dplyr::rename(cds_length = cds_len, tx_length = tx_len)


output_name <- paste0(str_split(opt$gtf, ".gtf")[[1]][1], ".longest_cds.transcript_info.tsv.gz")

fwrite(tx.info.dt, output_name, sep = "\t")


# cds.gr <- cds(txdb, use.names = TRUE, columns = "tx_name")
# transcripts.gr <- transcripts(txdb, use.names = TRUE)
# # names(transcripts.gr) <- id2name(txdb, "tx")
# cds_tx.gr <- mapToTranscripts(x=cds.gr, transcripts=transcripts.gr)
