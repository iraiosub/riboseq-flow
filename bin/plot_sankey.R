#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(networkD3))
suppressPackageStartupMessages(library(optparse))

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-d", "--dir"), action = "store", type = "character", default=NA, help = "results directory"),
                    make_option(c("-p", "--prefix"), action = "store", type = "character", default=NA, help = "sample and output file prefix"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# =========
# Functions
# =========



