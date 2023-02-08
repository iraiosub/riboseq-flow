#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { UMITOOLS_DEDUPLICATE as DEDUPLICATE_GENOME } from '../modules/local/umitools.nf' addParams(dedup_mode: params.dedup_genome)
include { UMITOOLS_DEDUPLICATE as DEDUPLICATE_TRANSCRIPTOME } from '../modules/local/umitools.nf' addParams(dedup_mode: params.dedup_transcriptome)


// Remove duplicate reads from BAM file based on UMIs

workflow DEDUPLICATE {

    take:
    genome_bam
    transcriptome_bam

    main:

    // Deduplication of sequences aligned to genome
    DEDUPLICATE_GENOME(
        genome_bam
        
    )

    // Deduplication of sequences aligned to transcriptome
    DEDUPLICATE_TRANSCRIPTOME(
        transcriptome_bam
    )

    emit:

    // Genome dedup
    dedup_genome_bam = DEDUPLICATE_GENOME.out.dedup_bam
    dedup_genome_bed = DEDUPLICATE_GENOME.out.dedup_bed

    // Transcriptome dedup
    dedup_transcriptome_bam = DEDUPLICATE_TRANSCRIPTOME.out.dedup_bam
    dedup_transcriptome_bed = DEDUPLICATE_TRANSCRIPTOME.out.dedup_bed

}


// 
//             DEDUP_UMI_UMITOOLS_TRANSCRIPTOME (
//                 ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
//                 params.umitools_dedup_stats
//             )
