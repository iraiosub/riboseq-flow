#!/usr/bin/env Rscript

# Prepare the shared exon map and chromosome sizes for genomic P-site tracks
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-g", "--gtf"), type = "character"),
                    make_option(c("-f", "--fai"), type = "character"),
                    make_option(c("-t", "--transcript_info"), type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))

txdb <- makeTxDbFromGFF(opt$gtf, format = "gtf")
exons <- exonsBy(txdb, "tx", use.names = TRUE)
transcripts <- fread(opt$transcript_info, select = "transcript_id")$transcript_id
saveRDS(exons[names(exons) %in% transcripts], "psite_exons.rds")
fwrite(fread(opt$fai, select = 1:2), "genome.chrom.sizes", sep = "\t", col.names = FALSE)
