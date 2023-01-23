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
    args2 += " -g " + params.adapter_fiveprime

    if (params.adapter_fiveprime && params.adapter_fiveprime) {
        args3 = args + " -g " + params.adapter_fiveprime
        args3 += " -a " + params.adapter_threeprime

        if (params.times_trimmed < 2) {
            args3 += " -n " + params.times_trimmed + 1
        } else {
            args3 += " -n " + params.times_trimmed
        
        }
    }
    
    if (params.ts_trimming) {
        ts_args = args + " -a " + params.ts_adapter_threeprime
        ts_args +=  " -u " + trim_polyG
        ts_args += " -n " + params.times_trimmed
    }

    if (params.ts_trimming)
        """
        cutadapt $ts_args $reads > ${sample_id}.cutadapt.log
        """
    else if (params.adapter_threeprime && params.adapter_fiveprime && !params.ts_trimming)
        """
        cutadapt $args3 $reads > ${sample_id}.cutadapt.log
        """
    else if (params.adapter_threeprime && !params.adapter_fiveprime && !params.ts_trimming)
        """
        cutadapt $args1 $reads > ${sample_id}.cutadapt.log
        """
    else if (params.adapter_threeprime && !params.adapter_fiveprime && !params.ts_trimming)
        """
        cutadapt $args2 $reads > ${sample_id}.cutadapt.log
        """
    else 
        error "Read trimming is enabled, but adapter sequence is missing"


}