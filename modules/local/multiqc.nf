#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process MULTIQC {
    // tag "${workflow.runName}"
    label 'process_medium'

    // conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"


    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    // path('fastqc/*')
    // path('premapped/*')
    // path('mapped/*')
    // path('deduplicated/*')
    path(logs)
    
    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data

    script:

    config_file = "--config ${params.multiqc_config}"

    """
    multiqc -f $config_file .
    """
}
