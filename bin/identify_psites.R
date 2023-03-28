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
length_range <- options[4]

# Rscript --vanilla identify_psites.R bam_dir gtf fasta length_range

# Create folder for plots
dir.create("ribowaltz_qc")


export_psites <- function(name, df_list) {
  
  df <- df_list[[name]]
  df$sample <- name
  
  data.table::fwrite(df, paste0(getwd(), "/", name, ".psite.tsv.gz"), sep = "\t")
  return(df)
  
}

plot_metaheatmap <- function(name, df_list, annotation) {
  
  ends_heatmap <- rends_heat(df_list, annotation, sample = name, 
                             cl = 100,
                             utr5l = 25, cdsl = 40, utr3l = 25)
  
  ends_heatmap.gg <- ends_heatmap$plot +
    ggplot2::ylim(20,35)
  
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", name, ".ends_heatmap.pdf"), ends_heatmap.gg, dpi = 400)
  
  return(ends_heatmap.gg)
}


save_metaprofile_psite_plot <- function(sample_name, plots_ls) {
  
  plot <- plots_ls[[sample_name]]
  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", strsplit(sample_name, "plot_")[[1]][2], ".metaprofile_psite.pdf"), plot, dpi = 400, width = 12, height = 6) # save in wide format

}


plot_cu <- function(sample_name, psite_info_ls) {

  psite.ls <- psite_info_ls[sample_name]

  cu_barplot <- codon_usage_psite(psite.ls, annotation = annotation.dt, sample = sample_name,
                                        fasta_genome = TRUE, 
                                        fastapath = fasta,
                                        gtfpath = gtf,
                                        frequency_normalization = FALSE) 
  
  # plot <- plots_ls[[sample_name]]

  cu_barplot <- cu_barplot +
    ggplot2::theme(plot.background = ggplot2::element_blank(), 
                  panel.grid.minor = ggplot2::element_blank(),
                  panel.grid.major = ggplot2::element_blank())

  cu_barplot.gg.ls <- cu_barplot[!(names(cu_barplot) %in% c("dt", "plot_comparison")) ]


  ggplot2::ggsave(paste0(getwd(),"/ribowaltz_qc/", sample_name, ".codon_usage.pdf"), plot, dpi = 400, width = 10, height = 7)
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

name_of_bams <- lapply(bams, function(x) strsplit(x, ".transcriptome.dedup.sorted.bam")[[1]][1])
names(name_of_bams) <- lapply(bams, function(x) strsplit(x, ".bam")[[1]][1])

reads.ls <- bamtolist(bamfolder = getwd(), 
                      annotation = annotation.dt,
                      name_samples = unlist(name_of_bams))

# Get filtered reads: keep only the ones with periodicity evidence
filtered.ls <- length_filter(data = reads.ls,
                             length_filter_mode = "periodicity",
                             periodicity_threshold = 50)


# Additionally filter them by length
min_length <- as.integer(strsplit(length_range, ":")[[1]][1])
max_length <- as.integer(strsplit(length_range, ":")[[1]][2])
filtered.ls <- length_filter(data = filtered.ls,
                             length_filter_mode = "custom",
                             length_range = min_length:max_length)

# =========
# Find P-sites
# =========

# Identify P sites, produce plots for each read length
psite_offset.dt <- psite(filtered.ls, flanking = 6, extremity = "auto", 
                         plot = TRUE, plot_format = "pdf")

data.table::fwrite(psite_offset.dt, 
                   paste0(getwd(), "/psite_offset.tsv.gz"), sep = "\t")

filtered_psite.ls <- psite_info(filtered.ls, site = "psite", 
                                offset = psite_offset.dt,
                                fasta_genome = TRUE, refseq_sep = " ",
                                fastapath = fasta,
                                gtfpath = gtf)

# Save psite info for each sample
lapply(names(filtered_psite.ls), export_psites, df_list = filtered_psite.ls)

# codon_coverage computes the number of read footprints or P-sites mapping on each triplet of annotated coding sequences and UTRs. 
# Such data can be exploited to generate occupancy profiles at codon resolution showing the abundance of RPFs along single transcripts

# This function computes transcript-specific codon coverages, defined as the number of either read footprints or P-sites mapping on each triplet of coding sequences and UTRs
codon_coverage_rpf.dt <- codon_coverage(filtered_psite.ls, psite = FALSE, annotation = annotation.dt)
codon_coverage_psite.dt <- codon_coverage(filtered_psite.ls, psite = TRUE, annotation = annotation.dt)

data.table::fwrite(codon_coverage_rpf.dt, 
                   paste0(getwd(),"/codon_coverage_rpf.tsv.gz"), sep = "\t")

data.table::fwrite(codon_coverage_psite.dt, 
                    paste0(getwd(),"/codon_coverage_psite.tsv.gz"), sep = "\t")

# Compute the number of P-sites mapping on annotated coding sequences or whole transcripts. 
# Such data can be used as starting point for downstream quantitative analyses (e.g. differential analyses) based on ribosome protected fragments

# By default, only in-frame P-sites falling in annotated coding sequences are considered and 
# no nucleotides at the beginning or at the end of the CDSs are excluded for restricting the analysis to a portion of the original coding sequences. 
# These settings can be modifyed through the parameters in_frame, start_nts and start_nts. 
# Moreover, the parameter whole_transcript specifies if whole transcripts should be considered instead of the annotated coding sequence

cds_coverage_psite.dt <- cds_coverage(filtered_psite.ls, annotation = annotation.dt)

data.table::fwrite(cds_coverage_psite.dt, 
                   paste0(getwd(), "/cds_coverage_psite.tsv.gz"), sep = "\t")


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
# save rds

# Metaheatmaps: the abundance of the 5' and 3' extremity of reads mapping on and around the start and the stop codon of annotated CDSs, stratified by their length.
ends_heatmap.gg.ls <- lapply(names(reads.ls), plot_metaheatmap, df_list = reads.ls, annotation = annotation.dt)
# save rds

# Ribosome profiling data should highlight the CDS of transcripts as the region with the higher percentage of reads. 
# To verify this property the function region_psite computes the percentage of P-sites falling in the three annotated transcript regions (5' UTR, CDS and 3' UTR). The bar plot of the resulting values includes a bar called "RNAs" displaying the expected read distribution from a random fragmentation of RNA.
psite_region <- region_psite(filtered_psite.ls, annotation.dt)
psite_region.gg <- psite_region$plot + 
  ggplot2::theme(plot.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/psite_region.pdf"), psite_region.gg, dpi = 400)
# save rds

# A fundamental characteristic of ribosome profiling data is the trinucleotide periodicity of ribosome footprints along coding sequences. 

# Functions frame_psite_length and frame_psite show if, and to which extent, the identified P-sites results in codon periodicity on the CDS. 
# Both functions compute the percentage of P-sites falling in the three possible translation reading frames for 5’ UTRs, CDSs and 3’ UTRs with one difference: 
# frame_psite_length analyses all read lengths separately and generates a heatmap for each transcript region, while frame_psite processes all reads at once, returning three bar plots.

frames_stratified <- frame_psite_length(filtered_psite.ls, region = "all", cl = 100)
frames_stratified.gg <- frames_stratified$plot

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/frames_stratified.pdf"), frames_stratified.gg, dpi = 400)
# save rds

frames <- frame_psite(filtered_psite.ls, region = "all")
frames.gg <- frames$plot +
  ggplot2::theme(plot.background = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank())

ggplot2::ggsave(paste0(getwd(), "/ribowaltz_qc/frames.pdf"), frames.gg, dpi = 400)
# save rds

# Plots should show an enrichment of P-sites in the first frame on the coding sequence but not the UTRs, as expected for ribosome protected fragments from protein coding mRNAs.

# Trinucleotide periodicity along coding sequences is provided by the function metaprofile_psite. 
# It generates metaprofiles (the merge of single, transcript-specific profiles) based on P-sites mapping around the start and the stop codon of annotated CDSs.
metaprofile <- metaprofile_psite(filtered_psite.ls, annotation.dt, sample = names(filtered_psite.ls),
                                         utr5l = 25, cdsl = 40, utr3l = 25,
                                         plot_title = "sample.transcript")

metaprofiles.gg.ls <- metaprofile[names(metaprofile) != "dt"]
lapply(names(metaprofiles.gg.ls), save_metaprofile_psite_plot, plots_ls = metaprofiles.gg.ls)
# save rds

# Codon usage

# lapply(names(filtered_psite.ls), plot_cu, psite_info_ls = filtered_psite.ls)