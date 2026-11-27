#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
packages <- args[-1]
versions <- c(R = as.character(getRversion()),
              setNames(vapply(packages, function(package) as.character(packageVersion(package)), character(1)), packages))
writeLines(c(paste0('"', args[1], '":'),
             sprintf('    %s: "%s"', names(versions), versions)), "versions.yml")
