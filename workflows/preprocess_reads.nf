#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CUTADAPT } from '../modules/cutadapt.nf'
include { UMITOOLS_EXTRACT } from '../modules/umitools.nf'

           
workflow PREPROCESS_READS {

    take:
    reads

    main:
    
    if (params.with_umi && !params.skip_umi_extract) {

        UMITOOLS_EXTRACT(reads)
        ch_reads = UMITOOLS_EXTRACT.out.fastq

        if (!params.skip_trimming) {

            CUTADAPT(UMITOOLS_EXTRACT.out.fastq)
            ch_reads = CUTADAPT.out.fastq
            
        }  

    } else if (!params.with_umi && !params.skip_umi_extract) {

        error "The reads have no UMI, but UMI extract is enabled! Make sure you provide the appropriate UMI options."

    } else {

        if (!params.skip_trimming) {

            CUTADAPT(reads)
            ch_reads = CUTADAPT.out.fastq
            
        }  else {

            ch_reads = reads
        }
    }

    emit:

    ch_reads

}





