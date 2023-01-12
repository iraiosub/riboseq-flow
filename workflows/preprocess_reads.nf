!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { CUTADAPT } from '../modules/cutadapt.nf'
include { UMITOOLS_EXTRACT } from '../modules/umitools.nf'
            
workflow GENERATE_REFERENCE_INDEX {

    take:
    reads


    main:
    
    if (params.with_umi && !params.skip_umi_extract) {

        UMITOOLS_EXTRACT(reads)

        if (!params.skip_trimming) {

            CUTADAPT(UMITOOLS_EXTRACT.out.fastq)
        }  

    } else {

        if (!params.skip_trimming) {

            CUTADAPT(reads)
        }  else {

            // keep the reads as they are
        }
    }



    emit:

    trimmed_fastq = CUTADAPT.out.fastq
    umi_extract = UMITOOLS_EXTRACT.fastq
    raw_reads = reads

}





