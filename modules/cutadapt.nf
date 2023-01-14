#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CUTADAPT {
    tag "${sample_id}"
    label 'process_high'

    conda 'bioconda::cutadapt=4.2'

    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: fastq
        path("*.cutadapt.log"), emit: log

    script:

    args = " -j ${task.cpus}"
    args += " -a " + params.adapter_threeprime
    args += " -g " + params.adapter_fiveprime
    args += " -n " + params.times_trimmed
    args += " -q " + params.min_quality
    args += " --minimum-length " + params.min_readlength
    args += " -o ${sample_id}.trimmed.fastq.gz"

    // if (params.adapter_threeprime && params.adapter_fiveprime)
    //     """
    //     cutadapt $args $reads > ${sample_id}.cutadapt.log
    //     """
    // else if (params.adapter_threeprime && !params.adapter_fiveprime)
    //     """
    //     cutadapt $args $reads > ${sample_id}.cutadapt.log
    //     """
    // else 
    //     error "No adapter provided"
    
    """
    cutadapt $args $reads > ${sample_id}.cutadapt.log
    """
}