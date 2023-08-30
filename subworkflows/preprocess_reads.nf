#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CUTADAPT } from '../modules/local/cutadapt.nf'
include { UMITOOLS_EXTRACT } from '../modules/local/umitools.nf'

// reads.into { reads_raw; reads_trimmed }


workflow PREPROCESS_READS {

    take:
    reads

    main:
    
    if (params.with_umi && !params.skip_umi_extract) {

        UMITOOLS_EXTRACT(reads)
        fastq = UMITOOLS_EXTRACT.out.fastq
        trimmed_fastq = UMITOOLS_EXTRACT.out.fastq

        if (!params.skip_preprocessing) {

            CUTADAPT(UMITOOLS_EXTRACT.out.fastq)
            fastq = CUTADAPT.out.fastq
            trimmed_fastq = CUTADAPT.out.trimmed_fastq
            logs = CUTADAPT.out.log
            
        }  

    } else if (!params.with_umi && !params.skip_umi_extract) {

        error "The reads have no UMI, but UMI extract is enabled! Make sure you provide the appropriate UMI options."

    } else {

        if (!params.skip_preprocessing) {

            CUTADAPT(reads)
            fastq = CUTADAPT.out.fastq
            trimmed_fastq = CUTADAPT.out.trimmed_fastq
            logs = CUTADAPT.out.log
            
        }  else {

            fastq = reads
            trimmed_fastq = reads
            logs = Channel.empty()
        }
    }

    emit:

    fastq
    trimmed_fastq
    logs

}

