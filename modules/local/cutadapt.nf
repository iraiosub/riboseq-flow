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

    publishDir "${params.outdir}/preprocessed", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_fastq // necessary for ribocutter
        tuple val(sample_id), path("${sample_id}.trimmed.filtered.fastq.gz"), emit: fastq
        tuple val(sample_id), path("*.log"), emit: log

    script:


    // Define core args
    args = " -j ${task.cpus}"
    args += " -q " + params.min_quality
    args += " -o ${sample_id}.trimmed.fastq.gz"

    args_filter = " -j ${task.cpus}"
    args_filter += " --minimum-length " + params.min_readlength
    args_filter += " -o ${sample_id}.trimmed.filtered.fastq.gz"

    // args_cut = " -j ${task.cpus}"
    args_cut += " -u " + params.cut_end
    // args_cut += " -o ${sample_id}.trimmed.filtered.fastq.gz"

    // Define option-specific args
    args1 = args + " -n " + params.times_trimmed
    args1 += " -a " + params.adapter_threeprime

    args2 = args + " -n " + params.times_trimmed
    args2 += " -g " + params.adapter_fiveprime

    if (params.adapter_fiveprime && params.adapter_fiveprime) {
        args3 = args + " -g " + params.adapter_fiveprime
        args3 += " -a " + params.adapter_threeprime

        // If user forgot to adjust n to >=2 when both 5' and 3' adaptor need trimming, adjust n to 2
        if (params.times_trimmed < 2) {
            args3 += " -n " + params.times_trimmed + 1
        } else {
            args3 += " -n " + params.times_trimmed
        
        }
    }
    
    if (params.ts_trimming) {
        ts_args = args + " -a " + params.ts_adapter_threeprime
        ts_args += " -n " + params.times_trimmed

        ts_args_filter = args_filter + " -u " + params.trim_polyG
    }

    if (params.ts_trimming && !params.adapter_threeprime && !params.adapter_fiveprime)
        """
        cutadapt $ts_args $reads > ${sample_id}.cutadapt_trim.log
        cutadapt $ts_args_filter ${sample_id}.trimmed.fastq.gz > ${sample_id}.cutadapt_filter.log

        """
    else if (params.adapter_threeprime && params.adapter_fiveprime && !params.ts_trimming)
        """
        cutadapt $args3 $reads > ${sample_id}.cutadapt_trim.log
        cutadapt $args_filter $args_cut ${sample_id}.trimmed.fastq.gz > ${sample_id}.cutadapt_filter.log
        """
    else if (params.adapter_threeprime && !params.adapter_fiveprime && !params.ts_trimming)
        """
        cutadapt $args1 $reads > ${sample_id}.cutadapt_trim.log
        cutadapt $args_filter $args_cut ${sample_id}.trimmed.fastq.gz > ${sample_id}.cutadapt_filter.log
        """
    else if (params.adapter_threeprime && !params.adapter_fiveprime && !params.ts_trimming)
        """
        cutadapt $args2 $reads > ${sample_id}.cutadapt_trim.log
        cutadapt $args_filter $args_cut ${sample_id}.trimmed.fastq.gz > ${sample_id}.cutadapt_filter.log
        """
    else if (!params.ts_trimming && !params.adapter_threeprime && !params.adapter_fiveprime)
        """
        cutadapt $args $reads > ${sample_id}.cutadapt_trim.log
        cutadapt $args_filter $args_cut ${sample_id}.trimmed.fastq.gz > ${sample_id}.cutadapt_filter.log
        """

    else 
        error "Read trimming is enabled, but adapter sequence is missing"

}