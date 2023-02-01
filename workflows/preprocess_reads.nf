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
        fastq = UMITOOLS_EXTRACT.out.fastq

        if (!params.skip_trimming) {

            CUTADAPT(UMITOOLS_EXTRACT.out.fastq)
            fastq = CUTADAPT.out.fastq
            // log = CUTADAPT.out.log
            
        }  

    } else if (!params.with_umi && !params.skip_umi_extract) {

        error "The reads have no UMI, but UMI extract is enabled! Make sure you provide the appropriate UMI options."

    } else {

        if (!params.skip_trimming) {

            CUTADAPT(reads)
            fastq = CUTADAPT.out.fastq
            // log = CUTADAPT.out.log
            
        }  else {

            fastq = reads
        }
    }

    emit:

    fastq

}

