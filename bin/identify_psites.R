#!/usr/bin/env Rscript

# Calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data
# Author: Ira Iosub
# Usage: get_p_sites.R -b $bam_folder -g $gtf -f $fasta -l $length_range

suppressPackageStartupMessages(library(riboWaltz))
# suppressPackageStartupMessages(library(optparse))

# =========
# Options and paths
# =========

# option_list <- list(make_option(c("-b", "--bam"), action = "store", type = "character", default=NA, help = "path to folder containing bam files"),
#                     make_option(c("-g", "--gtf"), action = "store", type = "character", default=NA, help = "GTF file"),
#                     make_option(c("-f", "--fasta"), action = "store", type = "character", default=NA, help = "genome FASTA file"),
#                     make_option(c("-l", "--length_range"), action = "store", type = "character", default="26:31", help = "string specifying the min and max length filter for RPFs, e.g.26:31"))
#
# opt_parser = OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)

options <- commandArgs(trailingOnly = TRUE)

bam_dir <- options[1]
gtf <- options[2]
fasta <- options[3]

# Filter out reads outside this length range for P-site indetification
length_range <- options[4]

# Filter out read lengths that below periodicity threshold for P-site indetification (good default = 50)
periodicity_thresh <- as.integer(options[5])

# P-site identification method
method <- options[6]
# Optionally exclude reads near either initiating or terminating ribosome which do not behave like elongating ribosomes. e.g. Ingolia CDS = +15th codon of the CDS to -10th codon of CDS
exclude_start <- as.numeric(options[7])
exclude_stop <- as.numeric(options[8])

longest_cds.df <- options[9]

# Usage: Rscript --vanilla identify_psites.R bam_dir gtf fasta length_range

# Create folder for plots
dir.create("ribowaltz_qc")


export_psites <- function(name, df_list) {
  
  df <- df_list[[name]]
  df$sample <- name
  
  data.table::fwrite(df, paste0(getwd(), "/", name, ".psite.tsv.gz"), sep = "\t")
  return(df)
  
}


plot_length_bins <- function(sample_name, df_list) {

  comparison_list <- list()
  comparison_list[["start_codon"]] <- df_list[[sample_name]][end5 <= cds_start & end3 >= cds_start]
  comparison_list[["whole_sample"]] <- df_list[[sample_name]]
  
  if(nrow(comparison_list[["start_codon"]]) == 0) {
    
    comparison_list <- list()
    comparison_list[["whole_sample"]] <- df_list[[sample_name]]
    
    rpf_list <- list("All" = c("whole_sample"))
    
    length_dist_split <-  rlength_distr(comparison_list,
                                        sample = rpf_list,
                                        multisamples = "average",
                                        plot_style = "split",
                                        colour = c("gray70"))
    
  } else {
    
    rpf_list <- list("Only_start" = c("start_codon"), "All" = c("whole_sample"))
    
    length_dist_split <-  rlength_distr(comparison_list,
                                        sample = rpf_list,
                                        multisamples = "average",
                                        plot_style = "split",
                                        colour = c("#699FC4", "gray70"))
  }
  
  length_dist_split.gg <- length_dist_split$plot +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
  
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".length_bins_for_psite.pdf"), length_dist_split.gg, dpi = 400, width = 10, height = 5)
  
}


plot_metaheatmap <- function(name, df_list, annotation) {
  
  ends_heatmap <- rends_heat(df_list, annotation, sample = name, 
                             cl = 100,
                             utr5l = 25, cdsl = 40, utr3l = 25)
  
  ends_heatmap.gg <- ends_heatmap$plot +
    ggplot2::ylim(20,45)
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", name, ".ends_heatmap.pdf"), ends_heatmap.gg, dpi = 400, width = 12, height = 8)
  
  return(ends_heatmap.gg)
}


save_metaprofile_psite_plot <- function(sample_name, plots_ls) {
  
  plot <- plots_ls[[sample_name]]
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", strsplit(sample_name, "plot_")[[1]][2], ".metaprofile_psite.pdf"), plot, dpi = 400, width = 12, height = 6) # save in wide format

}


plot_codon_usage <- function(sample_name, psite_info_ls) {

  psite.ls <- psite_info_ls[sample_name]

  cu_barplot <- codon_usage_psite(psite.ls, annotation = annotation.dt, sample = sample_name,
                                        fasta_genome = TRUE, 
                                        fastapath = fasta,
                                        gtfpath = gtf,
                                        frequency_normalization = FALSE) 
  
  # plot <- plots_ls[[sample_name]]
  # cu_barplot.gg <- unlist(cu_barplot[!(names(cu_barplot) %in% c("dt", "plot_comparison")) ])

  cu_barplot.gg <-cu_barplot$plot +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                  panel.grid.minor = ggplot2::element_blank(),
                  panel.grid.major = ggplot2::element_blank())

  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".codon_usage.pdf"), cu_barplot.gg, dpi = 400, width = 10, height = 7)
}

exclude_samples <- function(sample_name, df_list) {
  
  sample_list <- list()
  exclude <- c()
  sample_list[["start_codon"]] <- df_list[[sample_name]][end5 <= cds_start & end3 >= cds_start]
  
  if(nrow(sample_list[["start_codon"]]) == 0) {
    message("No reads overlapping start codon. Removing sample from analysis.")
    
    exclude <- sample_name
    
  } else {
    
    message("This sample will not be excluded")
  }
  
  return(exclude)
}


stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


# =========
# Load data and prepare annotation
# =========

# Prepare annotation: for each transcript, obtain total length, 5'UTR, CDS and 3'UTR length, respectively.
annotation.dt <- create_annotation(gtf)
data.table::fwrite(annotation.dt, 
                   paste0(getwd(), "/", strsplit(gtf, "gtf")[[1]][1],"transcript_info.tsv.gz"),
                   sep = "\t")

# Load BAM files

# We first define the "name_samples" character string vector as follow:
# Load BAM files
# bams <- list.files(bam_dir, pattern = ".transcriptome.dedup.sorted.bam")

bams <- as.list(strsplit(bam_dir, ",")[[1]])
name_of_bams <- lapply(bams, function(x) strsplit(basename(x), ".transcriptome.dedup.sorted.bam")[[1]][1])

# In case no UMIs were used
name_of_bams <- lapply(name_of_bams, function(x) strsplit(basename(x), ".Aligned.toTranscriptome.sorted.out.bam")[[1]][1])

names(name_of_bams) <- lapply(bams, function(x) strsplit(basename(x), ".bam")[[1]][1])

# Get number of samples
sample_count <- length(bams)

# Load bams
reads.ls <- bamtolist(bamfolder = getwd(), 
                      annotation = annotation.dt,
                      name_samples = unlist(name_of_bams))

# Order named list alphabetically
reads.ls <- reads.ls[order(names(reads.ls))]

# Filter to only transcript with longest CDS per gene
tx <- data.table::fread(longest_cds.df)

reads.ls <- lapply(reads.ls, function(df, tx.df) {
  df <- df[df$transcript %in% unique(tx.df$transcript_id),]
  return(df)
}, tx.df = tx)


# Get filtered reads: keep only the ones with periodicity evidence, periodicity_threshold = 50
filtered.ls <- length_filter(data = reads.ls,
                             length_filter_mode = "periodicity",
                             periodicity_threshold = periodicity_thresh)

# Additionally filter them by length
min_length <- as.integer(strsplit(length_range, ":")[[1]][1])
max_length <- as.integer(strsplit(length_range, ":")[[1]][2])

length_range <- min_length:max_length


# Remove sample if no reads left
filtered.ls <- Filter(function(x) dim(x)[1] > 0, filtered.ls)

filtered.ls <- length_filter(data = filtered.ls,
                             length_filter_mode = "custom",
                             length_range = min_length:max_length)

filtered.ls <- Filter(function(x) dim(x)[1] > 0, filtered.ls)

# Plot length bins used for P-site assignment
lapply(names(filtered.ls), plot_length_bins, df_list = filtered.ls)

# Filter out sample if no reads pass filtering, and stop analysis if no sample passes filering
exclude.ls <- lapply(names(filtered.ls), exclude_samples, df_list = filtered.ls)
filtered.ls <- filtered.ls[!names(filtered.ls) %in% exclude.ls]

if (length(filtered.ls) == 0) {

  message("No sample has reads passing filters for P-site identification. Stopping analysis.")
  stop_quietly()
}



# =========
# Identify P-sites
# =========

# Identify the exact position of the ribosome P-site within
# each read, determined by the localisation of its first nucleotide

if (method == "global_max_5end") {
  
  message("P site offset = distance between the 1st nucleotide of the TIS and the nt corresponding to the global maximum found in the read length-specific 5' end profiles")
  
  psite_offset.dt <- psite(filtered.ls, flanking = 6, extremity = "5end", 
                           plot = FALSE, plot_format = "pdf")
  
  
  # Replace ribowaltz-corrected offsets with temporary offests
  psite_offset.dt$corrected_offset_from_5 <- psite_offset.dt$offset_from_5
  psite_offset.dt$corrected_offset_from_3 <- psite_offset.dt$offset_from_3

  data.table::fwrite(psite_offset.dt, paste0(getwd(), "/psite_offset.tsv.gz"), sep = "\t")

  
  # Update reads information according to the inferred P-sites
  filtered_psite.ls <- psite_info(filtered.ls, psite_offset.dt, site = "psite",
                                  fasta_genome = TRUE, refseq_sep = " ",
                                  fastapath = fasta,
                                  gtfpath = gtf)
  
} else if (method == "ribowaltz") {
  
  message("P site offset = defined by the ribowaltz method")
  
  # Compute P-site offsets: temporary and corrected
  psite_offset.dt <- psite(filtered.ls, flanking = 6, extremity = "auto", 
                           plot = TRUE, plot_format = "pdf")

  data.table::fwrite(psite_offset.dt, paste0(getwd(), "/psite_offset.tsv.gz"), sep = "\t")
  
  # Update reads information according to the inferred P-sites
  filtered_psite.ls <- psite_info(filtered.ls, psite_offset.dt, site = "psite",
                                  fasta_genome = TRUE, refseq_sep = " ",
                                  fastapath = fasta,
                                  gtfpath = gtf)

} else {
  
  stop("Incorrect P-site offset method option. Available options are 'global_max_5end` or 'rbowaltz'")
}


# Save psite info for each sample
lapply(names(filtered_psite.ls), export_psites, df_list = filtered_psite.ls)

# codon_coverage computes the number of read footprints or P-sites mapping on each triplet of annotated coding sequences and UTRs. 
# Such data can be exploited to generate occupancy profiles at codon resolution showing the abundance of RPFs along single transcripts

# This function computes transcript-specific codon coverages, defined as the number of either read footprints or P-sites mapping on each triplet of coding sequences and UTRs
codon_coverage_rpf.dt <- codon_coverage(filtered_psite.ls, psite = FALSE, annotation = annotation.dt)
codon_coverage_psite.dt <- codon_coverage(filtered_psite.ls, psite = TRUE, annotation = annotation.dt)

data.table::fwrite(codon_coverage_rpf.dt, paste0(getwd(),"/codon_coverage_rpf.tsv.gz"), sep = "\t")
data.table::fwrite(codon_coverage_psite.dt, paste0(getwd(),"/codon_coverage_psite.tsv.gz"), sep = "\t")

# Compute the number of P-sites mapping on annotated coding sequences or whole transcripts. 
# Such data can be used as starting point for downstream quantitative analyses (e.g. differential analyses) based on ribosome protected fragments

# By default, only in-frame P-sites falling in annotated coding sequences are considered and 
# no nucleotides at the beginning or at the end of the CDSs are excluded for restricting the analysis to a portion of the original coding sequences. 
# These settings can be modifyed through the parameters in_frame, start_nts and start_nts. 
# Moreover, the parameter whole_transcript specifies if whole transcripts should be considered instead of the annotated coding sequence

# All over the CDS
cds_coverage_psite.dt <- cds_coverage(filtered_psite.ls, annotation = annotation.dt)

# Compute the number of in-frame P-sites per coding sequences exluding
# the first 15 codons and the last 10 codons.
# exclude_start <- 14*3
# exclude_stop <- 9*3

cds_coverage_psite_window.dt <- cds_coverage(filtered_psite.ls, annotation = annotation.dt, start_nts = exclude_start, stop_nts = exclude_stop)

# Export CDS coverage tables
data.table::fwrite(cds_coverage_psite.dt, paste0(getwd(), "/cds_coverage_psite.tsv.gz"), sep = "\t")
data.table::fwrite(cds_coverage_psite_window.dt, paste0(getwd(), "/cds_","plus", exclude_start, "nt_minus", exclude_stop, "nt_coverage_psite.tsv.gz"), sep = "\t")


# =========
# Diagnostic plots
# =========

# Read lengths averaged across samples
length_dist <- rlength_distr(reads.ls, sample = names(reads.ls), multisamples = "average", cl = 99)

length_dist.gg <- length_dist$plot +
  ggplot2::theme(legend.position = "none", legend.title=ggplot2::element_blank()) +
  ggplot2::scale_fill_manual(values = "grey70") +
  ggplot2::scale_color_manual(values = "grey30") +
  ggplot2::theme(plot.background = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank())

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/length_distribution.pdf"), length_dist.gg, dpi = 400)


# Metaheatmaps: the abundance of the 5' and 3' extremity of reads mapping on and around the start and the stop codon of annotated CDSs, stratified by their length.
ends_heatmap.gg.ls <- lapply(names(reads.ls), plot_metaheatmap, df_list = reads.ls, annotation = annotation.dt)

# Ribosome profiling data should highlight the CDS of transcripts as the region with the higher percentage of reads. 
# To verify this property the function region_psite computes the percentage of P-sites falling in the three annotated transcript regions (5' UTR, CDS and 3' UTR). The bar plot of the resulting values includes a bar called "RNAs" displaying the expected read distribution from a random fragmentation of RNA.
psite_region <- region_psite(filtered_psite.ls, annotation.dt)
psite_region.gg <- psite_region$plot + 
  ggplot2::theme(plot.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/psite_region.pdf"), psite_region.gg, dpi = 400, width = 10)


# A fundamental characteristic of ribosome profiling data is the trinucleotide periodicity of ribosome footprints along coding sequences. 

# Functions frame_psite_length and frame_psite show if, and to which extent, the identified P-sites results in codon periodicity on the CDS. 
# Both functions compute the percentage of P-sites falling in the three possible translation reading frames for 5’ UTRs, CDSs and 3’ UTRs with one difference: 
# frame_psite_length analyses all read lengths separately and generates a heatmap for each transcript region, while frame_psite processes all reads at once, returning three bar plots.

frames_stratified <- frame_psite_length(filtered_psite.ls, region = "all", length_range = min_length:max_length)
frames_stratified.gg <- frames_stratified$plot +
  ggplot2::scale_y_continuous(limits = c(min_length - 0.5, max_length + 0.5), breaks = seq(min(min_length + ((min_length) %% 2), max_length), max(min_length + ((min_length) %% 2), max_length), 
                    by = max(2, floor((max_length - min_length) / 7))))

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/frames_stratified.pdf"), frames_stratified.gg, dpi = 600, height = 24 , width = 18)


frames <- frame_psite(filtered_psite.ls, region = "all", length_range = min_length:max_length)
frames.gg <- frames$plot +
  ggplot2::theme(plot.background = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank())

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/frames.pdf"), frames.gg, dpi = 600, height = 24 , width = 18)

# Plots should show an enrichment of P-sites in the first frame on the coding sequence but not the UTRs, as expected for ribosome protected fragments from protein coding mRNAs.

# Trinucleotide periodicity along coding sequences is provided by the function metaprofile_psite. 
# It generates metaprofiles (the merge of single, transcript-specific profiles) based on P-sites mapping around the start and the stop codon of annotated CDSs.
metaprofile <- metaprofile_psite(filtered_psite.ls, annotation.dt, sample = names(filtered_psite.ls),
                                         utr5l = 25, cdsl = 40, utr3l = 25,
                                         plot_title = "sample.transcript")

metaprofiles.gg.ls <- metaprofile[names(metaprofile) != "dt"]
lapply(names(metaprofiles.gg.ls), save_metaprofile_psite_plot, plots_ls = metaprofiles.gg.ls)

# Codon usage
lapply(names(filtered_psite.ls), plot_codon_usage, psite_info_ls = filtered_psite.ls)
