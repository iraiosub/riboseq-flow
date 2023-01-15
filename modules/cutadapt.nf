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
    args += " -q " + params.min_quality
    args += " --minimum-length " + params.min_readlength
    args += " -o ${sample_id}.trimmed.fastq.gz"

    args1 = args + " -n " + params.times_trimmed
    args1 += " -a " + params.adapter_threeprime

    args2 = args + " -n " + params.times_trimmed
    args2 += " -a " + params.adapter_fiveprime

    if (params.adapter_fiveprime && params.adapter_fiveprime) {
        args3 = args + " -g " + params.adapter_fiveprime
        args3 += " -a " + params.adapter_fiveprime
        args3 += " -n " + params.times_trimmed + 1
    }


    if (params.adapter_threeprime && params.adapter_fiveprime)
        """
        cutadapt $args3 $reads > ${sample_id}.cutadapt.log
        """
    else if (params.adapter_threeprime && !params.adapter_fiveprime)
        """
        cutadapt $args1 $reads > ${sample_id}.cutadapt.log
        """
    else if (params.adapter_threeprime && !params.adapter_fiveprime)
        """
        cutadapt $args2 $reads > ${sample_id}.cutadapt.log
        """
    else 
        error "Read trimming is enabled, but adapter sequence is missing"

}