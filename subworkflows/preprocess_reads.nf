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
    ch_versions = Channel.empty()
    
    if (params.with_umi && !params.skip_umi_extract) {

        UMITOOLS_EXTRACT(reads)
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)

        if (!params.skip_trimming) {

            CUTADAPT(UMITOOLS_EXTRACT.out.fastq)
            ch_versions = ch_versions.mix(CUTADAPT.out.versions)

            trimmed_fastq = CUTADAPT.out.trimmed_fastq
            cut_fastq = CUTADAPT.out.cut_fastq
            fastq = CUTADAPT.out.filtered_fastq
            logs = CUTADAPT.out.log
            
        }  else {
            // Trimming disabled, so all channels are the original reads with UMIs moved to the header
            trimmed_fastq = UMITOOLS_EXTRACT.out.fastq
            cut_fastq = UMITOOLS_EXTRACT.out.fastq
            fastq = UMITOOLS_EXTRACT.out.fastq
            logs = ch_optional
        }

    } else if (!params.with_umi && !params.skip_umi_extract) {

        error "The reads have no UMI, but UMI extract is enabled! Make sure you provide the appropriate UMI options."

    } else {

        // Skip UMI-extract 
        if (!params.skip_trimming) {

            CUTADAPT(reads)
            ch_versions = ch_versions.mix(CUTADAPT.out.versions)
            trimmed_fastq = CUTADAPT.out.trimmed_fastq
            cut_fastq = CUTADAPT.out.cut_fastq
            fastq = CUTADAPT.out.filtered_fastq
            logs = CUTADAPT.out.log
            
        }  else {
            // Trimming and UMI extract disabled, so all channels are the original reads 
            trimmed_fastq = reads
            cut_fastq = reads
            fastq = reads
            logs = ch_optional
        }
    }

    emit:

    trimmed_fastq
    cut_fastq
    fastq
    logs
    versions = ch_versions

}

