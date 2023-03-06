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
path("*guides.csv"), emit: guides
path("*stats.csv"), emit: stats


script:
"""
ribocutter -i $reads -o ${sample_id} -g 50 -r 1000000 –-save_stats

"""
}
