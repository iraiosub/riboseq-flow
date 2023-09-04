#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(networkD3))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-d", "--dir"), action = "store", type = "character", default=NA, help = "results directory"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Functions
# =========

logs.dir <- opt$dir

# Parse logs
all.logs <- list.files(logs.dir, full.names = TRUE)

cutadapt.log <- all.logs[str_detect(all.logs, ".cutadapt_filter.log")]
premap.log <- all.logs[str_detect(all.logs, ".premap.log")]
map.log <- all.logs[str_detect(all.logs, ".Log.final.out")]
dedup.log <- all.logs[str_detect(all.logs, ".dedup.log")]
pcoding.log <- all.logs[str_detect(all.logs, ".qc_results.tsv.gz")]

# Extract sample name
sample_id <- str_split(basename(map.log), ".Log.final.out")[[1]][1]
message("Analysing ", sample_id)

# Extract information from logs

## LENGTH FILTER, reads that were too short are filtered out
cutadapt.log <- readLines(cutadapt.log)
total.reads <- cutadapt.log[grep("^Total reads processed:", cutadapt.log)]
too.short <- cutadapt.log[grep("^Reads that were too short:", cutadapt.log)]

# Extract min length from cutadapt command
min.length <- cutadapt.log[grep("^Command line parameters:", cutadapt.log)]
min.length <- parse_number(str_split(min.length, "--minimum-length ")[[1]][2])

# Extract the read numbers from logs
total.reads <- parse_number(total.reads)
too.short <- parse_number(too.short)
passing_length_filter.reads <- total.reads - too.short


## PREMAP, reads that pre-map to rRNA, tRNA, etc are filtered out
premap.log <- readLines(premap.log)
input_premap.reads <- parse_number(premap.log[grep("reads; of these:", premap.log)])

# Sanity check the input reads for pre-mapping is identical to output from length filter
stopifnot(input_premap.reads == passing_length_filter.reads)

not_premapped.reads <- parse_number(premap.log[grep("aligned 0 times", premap.log)])
premapped.reads <- input_premap.reads - not_premapped.reads

## MAP, only the reads that mapped uniquely to the genome are kept
map.log <- readLines(map.log)
input_map.reads <- parse_number(map.log[grep("Number of input reads", map.log)])

# Sanity check the input reads for mapping is identical to output from pre-mapping filter
stopifnot(input_map.reads == not_premapped.reads)

uniquemap.reads <- parse_number(map.log[grep("Uniquely mapped reads number", map.log)])

multimap.reads <- parse_number(map.log[grep("Number of reads mapped to too many loci", map.log)])
mismatches.reads <- parse_number(map.log[grep("Number of reads unmapped: too many mismatches", map.log)])
tooshort.reads <- parse_number(map.log[grep("Number of reads unmapped: too short", map.log)])
other.reads <- parse_number(map.log[grep("Number of reads unmapped: other", map.log)])

# Sanity check unmapped + mapped add up to total
unmapped.reads <- multimap.reads + mismatches.reads + tooshort.reads + other.reads
stopifnot(unmapped.reads + uniquemap.reads == input_map.reads)

## DUPLICATED
dedup.log <- readLines(dedup.log)

input_dedup.reads <- dedup.log[grep("INFO Reads: Input Reads:", dedup.log)]
input_dedup.reads <- parse_number(str_split(input_dedup.reads, "INFO")[[1]][2])

# Sanity check input for dedup identical to unique output from map
stopifnot(input_dedup.reads == uniquemap.reads)

dedup.reads <- dedup.log[grep("INFO Number of reads out:", dedup.log)]
dedup.reads <- parse_number(str_split(dedup.reads, "INFO")[[1]][2])
dup.reads <- input_dedup.reads - dedup.reads

## MAP to pcoding transcripts
## Map to pcoding transcripts and of expected length range
pcoding.log <- fread(pcoding.log)
pcoding.reads <- as.integer(pcoding.log$useful_read_n)
not_pcoding.reads <- dedup.reads - pcoding.reads

expected.reads <- as.integer(pcoding.log$expected_length_n)
out_expected.reads <- pcoding.reads - expected.reads

expected.length <- as.character(pcoding.log$expected_length)
print(expected.length)

# Build the data for the actual Sankey plot
read_list <- list()

# Length filter
read_list['total_reads'] <- total.reads
read_list['passing_length_filter'] <- passing_length_filter.reads
read_list['too_short'] <- tooshort.reads

# Premap
read_list['premapped'] <- premapped.reads
read_list['not_premapped'] <- not_premapped.reads

# Map
read_list['uniquely_mapped'] <- uniquemap.reads
read_list['multi_mapped'] <- multimap.reads
read_list['unmapped_too_many_mismatches'] <- mismatches.reads
read_list['unmapped_too_short'] <- tooshort.reads
read_list['unmapped_other'] <- other.reads

# Dedup
read_list['deduplicated'] <- dedup.reads
read_list['duplicated'] <- dup.reads

# Pcoding
read_list['pcoding'] <- pcoding.reads
read_list['not_pcoding'] <- not_pcoding.reads
read_list['expected'] <- expected.reads
read_list['out_expected'] <- out_expected.reads

# PLOTTING

nodes <- data.frame(
  name=c(
    'Total',
    paste0('Too short - less than ', min.length, ' nt'),
    'Passing length filter',
    'Premapped',
    'Not premapped',
    'Unmapped - multimapped',
    'Umapped - too many mismatches',
    'Unmapped - too short',
    'Unmapped - other',
    'Uniquely mapped',
    'Duplicated',
    'Deduplicated',
    'Not mapped to protein-coding transcripts',
    'Mapped to protein-coding transcripts (useful)',
    'Not in the expected length range',
    paste0('Expected length - ', expected.length, " nt")
  ),
  group=c(
    'a',
    'b',
    'b',
    'c',
    'c',
    'd',
    'd',
    'd',
    'd',
    'd',
    'e',
    'e',
    'f',
    'f',
    'g',
    'g'
  )
)

links <- as.data.frame(rbind(
  c( 0,  1, read_list[['too_short']]),
  c( 0,  2, read_list[['passing_length_filter']]),
  c( 2,  3, read_list[['premapped']]),
  c( 2,  4, read_list[['not_premapped']]),
  c( 4,  5, read_list[['multi_mapped']]),
  c( 4,  6, read_list[['unmapped_too_many_mismatches']]),
  c( 4,  7, read_list[['unmapped_too_short']]),
  c( 4,  8, read_list[['unmapped_other']]),
  c( 4,  9, read_list[['uniquely_mapped']]),
  c( 9, 10, read_list[['duplicated']]),
  c( 9, 11, read_list[['deduplicated']]),
  c( 11, 12, read_list[['not_pcoding']]),
  c( 11, 13, read_list[['pcoding']]),
  c( 13, 14, read_list[['out_expected']]),
  c( 13, 15, read_list[['expected']])
))


colnames(links) <- c('source', 'target', 'value')
links$link_source <- nodes$name[links$source + 1]


# Add group hoping for single color for links
# links$group <- as.factor(c("reads_group"))

message("Plotting Sankey diagram for ", sample_id)

# my_color <- 'd3.scaleOrdinal().domain(["a", "b", "c", "d","e", "reads_group"]).range(["#F3ECD9", "#F0F1E3", "#D7E0D8", "#C6D5D0", "#889C9B", "grey"])'
my_color <- 'd3.scaleOrdinal().domain(["a", "b", "c", "d", "e", "f", "g"]).range(["#F3ECD9", "#F0F1E3", "#D7E0D8", "#C6D5D0", "#889C9B", "#7D7A70", "#5C625C"])'


p <- sankeyNetwork(
  Links=links,
  Nodes=nodes,
  Source='source',
  Target='target',
  Value='value',
  LinkGroup='link_source',
  NodeID='name',
  NodeGroup='group',
  units='reads',
  fontSize=15,
  nodeWidth=30,
  fontFamily="sans-serif",
  width = 1600, # outer width, in pixels
  height = 1000, # outer height, in pixels
  margin = c(top = 5, right = 1, bottom = 5, left = 1),
  nodePadding = 10, # vertical separation between adjacent nodes
  iterations = 10,
  colourScale=my_color,
  sinksRight = F# effectively prevent the placement algorithm from running, so your nodes will be ordered as they were in the original data
)

saveNetwork(p, paste0(sample_id, '_sankey.html'), selfcontained = FALSE)
