#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { MAPPING_LENGTH_ANALYSIS as MAPPING_LENGTH_ANALYSIS_BEFORE_DEDUP } from '../modules/local/mapping_length.nf' addParams(length_analysis_type: params.before_dedup_length_analysis)
include { MAPPING_LENGTH_ANALYSIS as MAPPING_LENGTH_ANALYSIS_AFTER_DEDUP } from '../modules/local/mapping_length.nf' addParams(length_analysis_type: params.after_dedup_length_analysis)
include { MAPPING_LENGTH_ANALYSIS as MAPPING_LENGTH_ANALYSIS_AFTER_PREMAP } from '../modules/local/mapping_length.nf' addParams(length_analysis_type: params.after_premap_length_analysis)


// Create channels for optional inputs
ch_optional = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, file(params.optional) ] }


// Remove duplicate reads from BAM file based on UMIs

workflow MAPPING_LENGTH_ANALYSES {

    take:
    bam
    reads
    dedup_bam
    unmapped_reads
    

    main:

    MAPPING_LENGTH_ANALYSIS_BEFORE_DEDUP(bam.join(reads))

    if (params.with_umi) {

        MAPPING_LENGTH_ANALYSIS_AFTER_DEDUP(dedup_bam.join(reads)) // reads not needed here, added to kkep cardinality
        after_dedup_length_analysis = MAPPING_LENGTH_ANALYSIS_AFTER_DEDUP.out.length_analysis
        
    } else {

        after_dedup_length_analysis = ch_optional.join(reads).map { [ it[1] ] }
    }


    if (!params.skip_premap) {

         MAPPING_LENGTH_ANALYSIS_AFTER_PREMAP(bam.join(unmapped_reads))
         after_premap_length_analysis = MAPPING_LENGTH_ANALYSIS_AFTER_PREMAP.out.length_analysis

    } else {

        after_premap_length_analysis = ch_optional.join(reads).map { [ it[1] ] }
    }

   

    emit:

    before_dedup_length_analysis = MAPPING_LENGTH_ANALYSIS_BEFORE_DEDUP.out.length_analysis
    after_premap_length_analysis
    after_dedup_length_analysis
    

}
