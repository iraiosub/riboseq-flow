#!/usr/bin/env Rscript

# Combine per-sample riboWaltz offsets and coverage tables.
suppressPackageStartupMessages(library(data.table))

input_dir <- commandArgs(trailingOnly = TRUE)[1]
files <- list.files(input_dir, pattern = "\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)

if (length(files) == 0) {
  stop("No samples have reads passing filters for P-site identification.")
}

for (filename in unique(basename(files))) {
  input_files <- files[basename(files) == filename]

  if (filename == "psite_offset.tsv.gz") {
    combined.df <- rbindlist(lapply(input_files, fread, colClasses = c(sample = "character")), use.names = TRUE)
    sample_order <- match(combined.df$sample, sort(unique(combined.df$sample)))
    combined.df <- combined.df[order(sample_order, combined.df$length)]
  } else {
    combined.df <- fread(input_files[1], colClasses = c(transcript = "character"))
    keys <- intersect(names(combined.df), c("transcript", "start", "end", "from_cds_start",
                                          "from_cds_stop", "region", "length_cds", "length_selection"))

    # Read one sample at a time to avoid retaining every input table in memory.
    for (input_file in input_files[-1]) {
      combined.df <- merge(combined.df, fread(input_file, colClasses = c(transcript = "character")), by = keys, all = TRUE)
    }

    setcolorder(combined.df, c(keys, sort(setdiff(names(combined.df), keys))))
    combined.df[is.na(combined.df)] <- 0
  }

  fwrite(combined.df, filename, sep = "\t", row.names = FALSE)
}
