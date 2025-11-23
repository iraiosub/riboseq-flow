#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process CUTADAPT {
    tag "${sample_id}"
    label 'process_high'

    // conda 'bioconda::cutadapt=4.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.2--py39hbf8eff0_0' :
        'quay.io/biocontainers/cutadapt:4.2--py39hbf8eff0_0' }"

    publishDir "${params.outdir}/preprocessed", pattern: "*.log", mode: 'copy', overwrite: true
    if (params.save_trimmed) publishDir "${params.outdir}/preprocessed", pattern: "*.fastq.gz", mode: 'copy', overwrite: true


    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_fastq // necessary for Ribocutter
        tuple val(sample_id), path("${sample_id}.cut.fastq.gz"), emit: cut_fastq // necessary for length analyses
        tuple val(sample_id), path("${sample_id}.filtered.fastq.gz"), emit: filtered_fastq // used as input for pre-mapping
        tuple val(sample_id), path("*.log"), emit: log

    script:

    // Define core args for trimming
    trim_args = " -j ${task.cpus} -q ${params.minimum_quality} -o ${sample_id}.trimmed.fastq.gz"

    // Append adapter-specific args
    adapter_args = ""
    if (params.adapter_threeprime) adapter_args += " -a ${params.adapter_threeprime}"
    if (params.adapter_fiveprime) adapter_args += " -g ${params.adapter_fiveprime}"
    if (params.adapter_threeprime && params.adapter_fiveprime && params.times_trimmed < 2) params.times_trimmed = 2
    adapter_args += " -n ${params.times_trimmed}"

    // Define args for cutting bases from the end of the reads
    cut_args = " -j ${task.cpus} -u ${params.cut_end} -o ${sample_id}.cut.fastq.gz"

    // Define args for filtering reads based on length
    filter_args = " -j ${task.cpus} --minimum-length ${params.minimum_length} -o ${sample_id}.filtered.fastq.gz"

    """
    cutadapt ${trim_args}${adapter_args} $reads > ${sample_id}.cutadapt_trim.log
    cutadapt ${cut_args} ${sample_id}.trimmed.fastq.gz > ${sample_id}.cutadapt_cut.log
    cutadapt ${filter_args} ${sample_id}.cut.fastq.gz > ${sample_id}.cutadapt_filter.log
    """

}


