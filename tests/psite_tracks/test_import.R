# Run from the repository root with riboWaltz 1.2.0 available:
# Rscript tests/psite_tracks/test_import.R
suppressPackageStartupMessages(library(riboWaltz))
suppressPackageStartupMessages(library(data.table))
stopifnot(packageVersion("riboWaltz") == "1.2.0")
for (expr in parse("bin/identify_psites.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("read_bam_with_names"))) eval(expr)
}

root <- tempfile("psite-import-")
dir.create(root)
annotation <- data.table(transcript = c("txP", "txM"), l_utr5 = 60L, l_cds = 210L)
for (sample in c("A", "B")) {
  lines <- c("@HD\tVN:1.6", "@SQ\tSN:txP\tLN:330", "@SQ\tSN:txM\tLN:330")
  for (transcript in annotation$transcript) {
    for (length in c(28L, 29L)) {
      center <- if (length == 28L) 47L else 46L
      for (shift in -1:1) {
        count <- if (shift == 0L) 80L else 5L
        for (i in seq_len(count)) {
          name <- paste(sample, transcript, length, shift, i, sep = "_")
          lines <- c(lines, paste(name, 0, transcript, center + shift, 255,
                                 paste0(length, "M"), "*", 0, 0,
                                 strrep("A", length), strrep("I", length), sep = "\t"))
        }
      }
      for (i in 1:20) {
        lines <- c(lines, paste(paste(sample, transcript, length, "CDS", i, sep = "_"),
                               0, transcript, 80 + 3*i, 255, paste0(length, "M"), "*", 0, 0,
                               strrep("A", length), strrep("I", length), sep = "\t"))
      }
    }
  }
  # Numeric-looking QNAMEs, reverse alignments, and the indel filter boundary
  lines <- c(lines, paste("0001", 0, "txP", 47, 255, "28M", "*", 0, 0, strrep("A", 28), strrep("I", 28), sep = "\t"),
             paste("reverse", 16, "txP", 47, 255, "28M", "*", 0, 0, strrep("A", 28), strrep("I", 28), sep = "\t"),
             paste("keep_indel", 0, "txP", 47, 255, "28M5I", "*", 0, 0, strrep("A", 33), strrep("I", 33), sep = "\t"),
             paste("drop_indel", 0, "txP", 47, 255, "28M6I", "*", 0, 0, strrep("A", 34), strrep("I", 34), sep = "\t"))
  sam <- file.path(root, paste0(sample, ".sam"))
  writeLines(lines, sam)
  Rsamtools::asBam(sam, destination = file.path(root, sample), overwrite = TRUE, indexDestination = FALSE)
}

original <- bamtolist(root, annotation, name_samples = c(A = "A", B = "B"))
retained <- setNames(lapply(c("A", "B"), function(sample) read_bam_with_names(file.path(root, paste0(sample, ".bam")), annotation)), c("A", "B"))
for (sample in names(original)) {
  stopifnot(identical(as.data.frame(original[[sample]]), as.data.frame(retained[[sample]][, !"read_name"])),
            "0001" %in% retained[[sample]]$read_name, "keep_indel" %in% retained[[sample]]$read_name,
            !any(retained[[sample]]$read_name %in% c("reverse", "drop_indel")))
}
original <- length_filter(original, length_filter_mode = "periodicity", periodicity_threshold = 50)
retained <- length_filter(retained, length_filter_mode = "periodicity", periodicity_threshold = 50)
offsets <- psite(original, flanking = 6, extremity = "auto", plot = FALSE)
new_offsets <- psite(retained, flanking = 6, extremity = "auto", plot = FALSE)
stopifnot(identical(as.data.frame(offsets), as.data.frame(new_offsets)))
original <- psite_info(original, offsets, site = NULL)
retained <- psite_info(retained, new_offsets, site = NULL)
for (sample in names(original)) {
  stopifnot(identical(as.data.frame(original[[sample]]), as.data.frame(retained[[sample]][, !"read_name"])))
}
writeLines(c("@HD\tVN:1.6", "@SQ\tSN:txP\tLN:330"), file.path(root, "empty.sam"))
Rsamtools::asBam(file.path(root, "empty.sam"), destination = file.path(root, "empty"), overwrite = TRUE, indexDestination = FALSE)
empty <- read_bam_with_names(file.path(root, "empty.bam"), annotation)
stopifnot(nrow(empty) == 0, "read_name" %in% names(empty))
unlink(root, recursive = TRUE)
cat("PASS: read names retained; import filters, offsets and P-site positions unchanged\n")
