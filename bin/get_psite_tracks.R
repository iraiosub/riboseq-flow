#!/usr/bin/env Rscript

# get genomic coordinates of Psites identified by ribowaltz
# Author: Ira Iosub

suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# =========
# Options and paths
# =========

option_list <- list(make_option(c("-p", "--psite"), action = "store", type = "character", default=NA, help = "P-sites table"),
                    make_option(c("-g", "--gtf"), action = "store", type = "character", default=NA, help = "GTF file"),
                    make_option(c("-f", "--fai"), action = "store", type = "character", default=NA, help = "FASTA index file"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

make_txdb <- function(gtf, org) {
  
  # Get name of db
  name <- paste0(str_split(gtf, ".gtf")[[1]][1], ".sqlite")
  
  if (file.exists(name)) {
    
    TxDb <- loadDb(name)
    
  } else {
    
    TxDb <- makeTxDbFromGFF(opt$gtf, format="gtf",
                            organism = org) # chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene
    
    saveDb(TxDb, file=name)
    # TxDb <- loadDb(name)
  }
  
  TxDb <- keepStandardChromosomes(TxDb, pruning.mode="coarse")
  return(TxDb)
}



# A function that makes a GRanges object for each sample from the psite tsv from ribowaltz
table_to_granges <- function(df) {
  
  df <- df %>%
    dplyr::select(transcript, psite, sample) %>%
    mutate(strand = "*", start = psite, end = psite, score = 1, seqnames = transcript)
  
  df.gr <- GRanges(df)
  names(df.gr) <- df.gr$transcript
  
  return(df.gr)
  
}

# A function that sums scores of identical position in a GRanges object
sum_scores <- function(gr) {
  
  hits <- findOverlaps(gr, gr)
  agg.gr <- aggregate(gr, hits, score=sum(score))
  gr$score <- agg.gr$score
  gr.unique <- gr[!duplicated(gr)]
  
  return(gr.unique)
  
}

# A function that converts transcript to genomic coordinates and exports bed and bigwig files
convert_coordinates <- function(tx.gr, annotation, chr_lengths) {
  
  genomic.gr <- mapFromTranscripts(tx.gr, annotation, ignore.strand = FALSE)
  genomic.gr$score <- 1L
  
  # Sum up scores for identical positions
  genomic.gr <- sum_scores(genomic.gr)
  export.bed(genomic.gr, paste0(unique(tx.gr$sample), ".psites.bed"))
  
  # Separate plus and minus strand and export bigwigs
  plus.gr <- genomic.gr[strand(genomic.gr) == "+" ]
  minus.gr <- genomic.gr[strand(genomic.gr) == "-" ]
  
  # Set chr lengths based on fai
  # seqlengths(plus.gr) <- chr_lengths[names(seqlengths(plus.gr))]
  # seqlengths(minus.gr) <- chr_lengths[names(seqlengths(minus.gr))]
  
  # export.bw(plus.gr, paste0(unique(tx.gr$sample), ".psites.forward.bigwig"))
  # export.bw(minus.gr, paste0(unique(tx.gr$sample), ".psites.reverse.bigwig"))
  
  return(genomic.gr)
  
}


# Load psites
# psites.ls <- as.list(strsplit(opt$psite, ",")[[1]])
# psites.df.ls <- lapply(psites.ls, fread)
# psites.df <- rbindlist(psites.df.ls)

psites.df <- opt$psite

# Load annotation and subset transcripts of interest
txdb <- make_txdb(opt$gtf, org = 'Homo sapiens')

# Load fai and extract chromosome lengths
fai <- fread(opt$fai)
chr_lengths <- fai$V2
names(chr_lengths) <- fai$V1

# Get structure of the transcripts in the form of a GRangesList
# If we use only tx.gr <- transcripts(txdb, columns = c("exon_id", "tx_name")), that will extract the pre-mRNAs coordinates for each isoform, the transcript structure (i.e. what is intron and what is exon is not automatically inferred)
exons.grl <- exonsBy(txdb, "tx", use.names = T)
exons.grl <- exons.grl[names(exons.grl) %in% unique(psites.df$transcript)]

# Get transcriptomic GRanges
# psites.tx.grl <- lapply(psites.df.ls, table_to_granges)
psites.tx.gr <- table_to_granges(psites.df)

# Mapping transcript coordinates to genomic coordinates
# psites.genomic.grl <- lapply(psites.tx.grl, convert_coordinates, annotation = exons.grl, chr_lengths = chr_lengths)

psites.genomic.gr <- convert_coordinates(psites.tx.gr, annotation = exons.grl, chr_lengths = chr_lengths)