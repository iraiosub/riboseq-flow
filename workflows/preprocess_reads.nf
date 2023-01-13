#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CUTADAPT } from '../modules/cutadapt.nf'
include { UMITOOLS_EXTRACT } from '../modules/umitools.nf'

process KEEP_RAW_READS {

    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true
  
    input:
        tuple val(sample_id), path(reads)
  
    output:
        tuple val(sample_id), path("*.fastq.gz"), emit: fastq
  
    script:
    """
    zcat $reads > ${sample_id}.fastq
    gzip ${sample_id}.fastq ${sample_id}.fastq.gz
    """
}

           
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
        }  else {

            KEEP_RAW_READS(UMITOOLS_EXTRACT.out.fastq)
            ch_reads = KEEP_RAW_READS.out.fastq
        }

    } else {

        if (!params.skip_trimming) {

            CUTADAPT(reads)
            ch_reads = CUTADAPT.out.fastq
            
        }  else {

            // keep the reads as they are
            ch_reads = KEEP_RAW_READS.out.fastq
        }
    }



    emit:

    ch_reads

}





