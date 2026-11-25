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

    publishDir "${params.outdir}/multiqc", saveAs: { filename -> filename.endsWith('versions.yml') ? null : filename }, mode: 'copy', overwrite: true
    publishDir "${params.outdir}/pipeline_info", pattern: "software_versions.yml", mode: 'copy', overwrite: true

    input:
    path(logs)
    path(multiqc_config)
    path(software_versions, stageAs: 'versions/??/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "versions.yml", emit: versions
    path "software_versions.yml", emit: software_versions

    script:

    // def config_file = params.multiqc_config ? "--config $params.multiqc_config" : ''
    def pipeline_version = [workflow.manifest.version, workflow.revision, workflow.commitId?.take(7)].findAll { it }.unique().join(' / ')

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        PyYAML: \$(python -c 'import yaml; print(yaml.__version__)')
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    collect_software_versions.py ${workflow.nextflow.version} "$pipeline_version" versions.yml $software_versions
    ln -s software_versions.yml software_mqc_versions.yml
    multiqc -f $multiqc_config .
    """
}
