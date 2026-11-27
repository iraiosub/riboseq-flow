#!/usr/bin/env Rscript

# Map inferred P-sites to the genome and count each read once
# Author: Ira Iosub
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-p", "--psite"), type = "character"),
                    make_option(c("-a", "--annotation"), type = "character"),
                    make_option(c("-s", "--sizes"), type = "character"),
                    make_option(c("-o", "--output_prefix"), type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))

psites <- fread(opt$psite, select = c("read_name", "transcript", "psite"),
                colClasses = list(character = c("read_name", "transcript")))
if (!"read_name" %in% names(psites) || anyNA(psites$read_name) || any(psites$read_name %in% c("", "*"))) {
  stop("Genomic P-site tracks require read names; regenerate P-sites with --genomic_psite_coverage")
}
input_records <- nrow(psites)
input_read_names <- uniqueN(psites$read_name)
exons <- readRDS(opt$annotation)
sizes <- fread(opt$sizes, col.names = c("chrom", "size"))

# Map each transcript position once, then restore read identities before counting
positions <- unique(psites[, .(transcript, psite)])
tx <- GRanges(seqnames = positions$transcript, ranges = IRanges(positions$psite, width = 1), strand = "*")
mapped <- mapFromTranscripts(tx, exons, ignore.strand = FALSE)
loci <- cbind(positions[mapped$xHits], data.table(chrom = as.character(seqnames(mapped)),
              position = start(mapped), strand = as.character(strand(mapped))))
rm(positions, tx, mapped, exons)
gc()
assignments <- merge(psites, loci, by = c("transcript", "psite"))
mapped_records <- nrow(assignments)
sites <- unique(assignments[, .(read_name, chrom, position, strand)])
duplicate_assignments <- mapped_records - nrow(sites)
rm(psites, loci, assignments)
gc()

# Transcript copies of the same read count once; conflicting P-sites are excluded
ambiguous <- sites[, .N, by = read_name][N > 1, read_name]
sites <- sites[!read_name %in% ambiguous]
if (any(!sites$chrom %in% sizes$chrom) ||
    any(sites$position < 1 | sites$position > sizes$size[match(sites$chrom, sizes$chrom)]) ||
    any(!sites$strand %in% c("+", "-"))) {
  stop("Mapped P-sites do not match the genome chromosome sizes or strands")
}
coverage <- sites[, .(count = .N), by = .(chrom, position, strand)]

stats <- data.table(sample = opt$output_prefix, input_records = input_records,
                    input_read_names = input_read_names, mapped_records = mapped_records,
                    unmapped_records = input_records - mapped_records,
                    duplicate_assignments = duplicate_assignments, ambiguous_reads = length(ambiguous),
                    forward_psites = sum(sites$strand == "+"), reverse_psites = sum(sites$strand == "-"))
fwrite(stats, paste0(opt$output_prefix, ".psites.mapping_stats.tsv"), sep = "\t")
if (length(ambiguous) > 0 || stats$unmapped_records > 0) {
  message("Excluded ", length(ambiguous), " reads with conflicting genomic P-sites; ",
          stats$unmapped_records, " input records could not be mapped. See mapping_stats.tsv.")
}
rm(sites, ambiguous)
gc()

# Keep the legacy BED output, with counts in its score column
genomic <- GRanges(seqnames = coverage$chrom, ranges = IRanges(coverage$position, width = 1),
                   strand = coverage$strand, score = coverage$count)
export.bed(genomic, paste0(opt$output_prefix, ".codon_coverage_psite.bed.gz"))

for (direction in c("forward", "reverse")) {
  strand_value <- if (direction == "forward") "+" else "-"
  track <- coverage[strand == strand_value, .(chrom, start = position - 1L, end = position, count)]
  if (nrow(track) > 0) {
    fwrite(track, paste0(opt$output_prefix, ".psites.", direction, ".unsorted.bedgraph"),
           sep = "\t", col.names = FALSE)
  }
}
