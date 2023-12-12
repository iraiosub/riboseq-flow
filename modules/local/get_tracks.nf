#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_COVERAGE_TRACKS {
    
    tag "${sample_id}"
    label 'process_high'

    // conda "bioconda::deeptools=3.5.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0' }"

    publishDir "${params.outdir}/coverage_tracks", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(bam), path(bai)
    
    output:
        tuple val(sample_id), path("*.forward*"), emit: forward
        tuple val(sample_id), path("*.reverse*"), emit: reverse

    script:
        
        binsize = params.bin_size
        format = params.track_format

        """

        bamCoverage \
            --bam $bam \
            --filterRNAstrand reverse \
            --numberOfProcessors ${task.cpus} \
            -bs $binsize \
            --outFileFormat $format \
            --outFileName ${sample_id}.forward.$format

        bamCoverage \
            --bam $bam \
            --filterRNAstrand forward \
            -bs $binsize \
            --outFileFormat $format \
            --numberOfProcessors ${task.cpus} \
            --outFileName ${sample_id}.reverse.$format
        """

}

