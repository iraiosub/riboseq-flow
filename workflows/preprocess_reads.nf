#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CUTADAPT } from '../modules/cutadapt.nf'
include { UMITOOLS_EXTRACT } from '../modules/umitools.nf'

process KEEP_RAW_READS {
  
  input:
    tuple val(sample_id), path(reads)
  
  output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq
  
  script:
  """
  mv $reads ${sample_id}.fastq.gz
  # touch instead?
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





