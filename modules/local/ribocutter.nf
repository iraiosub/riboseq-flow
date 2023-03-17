#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOCUTTER {

    tag "${sample_id}"
    label 'process_single'

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
    args += " " + params.ribocutter_args


    if (is.null(params.min_length)) {

        min_read_length_arg = ""
        suffix = ""

    } else {

        min_read_length_arg = " --min_read_length " + params.min_length
        suffix = params.min_length
    }

    
    """

    ribocutter -i $reads -o ${sample_id}.min${suffix} $args $min_read_length --save_stats

    """
}
