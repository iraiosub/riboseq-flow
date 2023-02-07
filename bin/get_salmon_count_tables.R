#!/usr/bin/env Rscript
library(tximport)
library(GenomicFeatures)
library(data.table)
library(optparse)
library(SummarizedExperiment)

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

build_table_from_se <- function(se.obj, slot) {
  
  table <- cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]])
  return(as.data.frame(table))
  
}


option_list <- list(make_option(c("", "--salmon_results"), action = "store", type = "character", default=NA, help = "folder with salmon quantifications"),
                    make_option(c("", "--gtf"), action = "store", type = "character", default=NA, help = "GTF annotation file"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create txdb from GTF
opt$gtf <- "gencode.v29.primary_assembly.annotation.gtf.gz"
opt$dir <- "../salmon_quant"
prefix <- "salmon_merged"

# Prepare annotation
gtf <- import.gff2(opt$gtf)

# Create tx2gene
txdb <- make_hg38_txdb(opt$gtf)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Load salmon quantif
files <- list.files(opt$dir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))

coldata <- data.frame(files = files, names = names(files))
rownames(coldata) = coldata[["names"]]

# Import transcripts
txi <- tximport(files, type = "salmon", txOut = TRUE)

rowdata <- tx2gene[match(rownames(txi[[1]]), as.character(tx2gene[["TXNAME"]])),]
rownames(rowdata) = rowdata[["TXNAME"]]

se = SummarizedExperiment(assays = list(counts = txi[["counts"]], abundance = txi[["abundance"]], length = txi[["length"]]),
                          colData = DataFrame(coldata),
                          rowData = rowdata)

if (!is.null(tx2gene)) {
  gi = summarizeToGene(txi, tx2gene = tx2gene)
  gi.ls = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM")
  gi.s = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="scaledTPM")
  
  growdata = unique(rowdata[,1:2])
  growdata = growdata[match(rownames(gi[[1]]), growdata[["GENEID"]]),]
  rownames(growdata) = growdata[["TXNAME"]]
  gse = SummarizedExperiment(assays = list(counts = gi[["counts"]], abundance = gi[["abundance"]], length = gi[["length"]]),
                             colData = DataFrame(coldata),
                             rowData = growdata)
  gse.ls = SummarizedExperiment(assays = list(counts = gi.ls[["counts"]], abundance = gi.ls[["abundance"]], length = gi.ls[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = growdata)
  gse.s = SummarizedExperiment(assays = list(counts = gi.s[["counts"]], abundance = gi.s[["abundance"]], length = gi.s[["length"]]),
                               colData = DataFrame(coldata),
                               rowData = growdata)
}

if(exists("gse")){
  fwrite(build_table_from_se(gse, "abundance"), paste0(prefix, ".gene_tpm.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
  fwrite(build_table_from_se(gse, "counts"), paste0(prefix, ".gene_counts.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
  fwrite(build_table_from_se(gse.ls, "abundance"), paste0(prefix, ".gene_tpm_length_scaled.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
  fwrite(build_table_from_se(gse.ls, "counts"), paste0(prefix, ".gene_counts_length_scaled.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
  fwrite(build_table_from_se(gse.s, "abundance"), paste0(prefix, ".gene_tpm_scaled.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
  fwrite(build_table_from_se(gse.s, "counts"), paste0(prefix, ".gene_counts_scaled.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
}

fwrite(build_table_from_se(se,"abundance"), paste0(prefix, ".transcript_tpm.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
fwrite(build_table_from_se(se, "counts"), paste0(prefix, ".transcript_counts.tsv.gz"), sep="\t", quote=FALSE, row.names = FALSE)
