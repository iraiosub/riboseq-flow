#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOCUTTER {
tag "${sample_id}"
label 'process_medium'

conda 'bioconda::ribocutter=0.1.1'

publishDir "${params.outdir}/ribocutter"

input:
tuple val(sample_id), path(reads)

output:
path("*.csv"), emit: guides
path("*stats.csv"), emit: stats


script:

 args = " -g " + params.guide_number
 args += " -r " + params.max_reads
 args += " --min_read_length " + params.min_read_length
 args += " " + params.extra_ribocutter_args

 
"""

ribocutter -i $reads -o ${sample_id} $args --save_stats

"""
}
