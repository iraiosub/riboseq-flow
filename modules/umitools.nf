#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process UMITOOLS_EXTRACT {
    tag "${sample_id}"
    label "process_low"
    publishDir "${params.outdir}/umi_extract", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.umi_extract.fastq.gz"), emit: fastq
        path("*.umi_extract.log"), emit: log

    script:
    args = " --bc-pattern=" + params.umi_pattern
    args += " --extract-method=" + params.umi_extract_method

    def args = task.ext.args ?: ''
    """
    umi_tools \
        extract \
        -I $reads \
        -S ${sample_id}.umi_extract.fastq.gz \
        $args \
        -L ${sample_id}.umi_extract.log
    """
}
