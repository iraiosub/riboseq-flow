#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CUTADAPT } from '../modules/local/cutadapt.nf'
include { UMITOOLS_EXTRACT } from '../modules/local/umitools.nf'

// reads.into { reads_raw; reads_trimmed }
// Create channel for optional input
ch_optional = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, [] ] }

workflow PREPROCESS_READS {

    take:
    reads

    main:
    
    if (params.with_umi && !params.skip_umi_extract) {

        UMITOOLS_EXTRACT(reads)
        // fastq = UMITOOLS_EXTRACT.out.fastq
        // trimmed_fastq = UMITOOLS_EXTRACT.out.fastq

        if (!params.skip_trimming) {

            CUTADAPT(UMITOOLS_EXTRACT.out.fastq)
            fastq = CUTADAPT.out.fastq
            trimmed_fastq = CUTADAPT.out.trimmed_fastq
            logs = CUTADAPT.out.log
            
        }  else {

            fastq = UMITOOLS_EXTRACT.out.fastq
            trimmed_fastq = UMITOOLS_EXTRACT.out.fastq
            logs = ch_optional
        }

    } else if (!params.with_umi && !params.skip_umi_extract) {

        error "The reads have no UMI, but UMI extract is enabled! Make sure you provide the appropriate UMI options."

    } else {

        if (!params.skip_trimming) {

            CUTADAPT(reads)
            fastq = CUTADAPT.out.fastq
            trimmed_fastq = CUTADAPT.out.trimmed_fastq
            logs = CUTADAPT.out.log
            
        }  else {

            fastq = reads
            trimmed_fastq = reads
            logs = ch_optional
        }
    }

    emit:

    fastq
    trimmed_fastq
    logs

}

