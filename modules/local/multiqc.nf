#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process MULTIQC {
    // tag "${workflow.runName}"
    label 'process_medium'

    // conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path(logs)
    path(multiqc_config)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "versions.yml", emit: versions

    script:

    // def config_file = params.multiqc_config ? "--config $params.multiqc_config" : ''

    """
    multiqc -f $multiqc_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
