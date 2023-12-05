#!/usr/bin/env nextflow

//
// Prepare all genome files for running the riboseq analysis pipeline
//


// Specify DSL2
nextflow.enable.dsl=2

include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_SMALLRNA_FASTA } from '../modules/nf-core/gunzip/main'
include { GENERATE_REFERENCE_INDEX } from '../subworkflows/generate_index.nf'
include { GET_TRANSCRIPT_INFO } from '../modules/local/transcript_info.nf'
include { GET_TRANSCRIPT_FASTA } from '../modules/local/transcript_info.nf'


workflow PREPARE_RIBOSEQ_REFERENCE {

    take:
    genome_fasta
    genome_gtf
    smallrna_fasta

    main:

    // Prepare annotation: unzip annotation and genome files if necessary
    ch_genome_fasta = Channel.empty()
    if (genome_fasta.toString().endsWith('.gz')) {
        ch_genome_fasta = GUNZIP_FASTA ( [ [:], genome_fasta ] ).gunzip
    } else {
        ch_genome_fasta = Channel.from( [ [ [:], genome_fasta ] ] )
    }


    ch_genome_gtf = Channel.empty()
    if (genome_gtf.toString().endsWith('.gz')) {
        ch_genome_gtf = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip
    } else {

        ch_genome_gtf = Channel.from( [ [ [:], genome_gtf ] ] )
    }

    ch_smallrna_fasta = Channel.empty()
    if (!params.skip_premap && smallrna_fasta.toString().endsWith('.gz')) {
        ch_smallrna_fasta = GUNZIP_SMALLRNA_FASTA ( [ [:], smallrna_fasta ] ).gunzip
    } else {
        // ch_smallrna_fasta = file(params.smallrna_fasta)
        ch_smallrna_fasta = Channel.from( [ [ [:], smallrna_fasta ] ] )
    }

    // Prepare annotation: create index for alignment
    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta, ch_genome_fasta, ch_genome_gtf)

    if (!params.transcript_info) {
        GET_TRANSCRIPT_INFO(ch_genome_gtf.map{ it[1] })
        ch_transcript_info = GET_TRANSCRIPT_INFO.out.transcript_info
        ch_transcript_info_gtf = GET_TRANSCRIPT_INFO.out.transcripts_gtf
        
        GET_TRANSCRIPT_FASTA(ch_genome_fasta.map{ it[1] }, ch_transcript_info_gtf)
        ch_transcript_info_fa = GET_TRANSCRIPT_FASTA.out.transcripts_fa

    } else {

        ch_transcript_info = Channel.fromPath(params.transcript_info, checkIfExists: true)
        ch_transcript_info_fa = Channel.fromPath(params.transcript_fasta, checkIfExists: true)
    }

    emit:

    smallrna_bowtie2_index = GENERATE_REFERENCE_INDEX.out.smallrna_bowtie2_index
    genome_star_index = GENERATE_REFERENCE_INDEX.out.genome_star_index
    genome_gtf = ch_genome_gtf
    genome_fasta = ch_genome_fasta
    smallrna_fasta = ch_smallrna_fasta
    transcript_info = ch_transcript_info
    // transcript_info_gtf = ch_transcript_info_gtf
    transcript_info_fa = ch_transcript_info_fa

}