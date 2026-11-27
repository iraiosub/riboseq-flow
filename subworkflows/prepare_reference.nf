#!/usr/bin/env nextflow

//
// Prepare all genome files for running the riboseq analysis pipeline
//


// Specify DSL2
nextflow.enable.dsl=2

include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_CONTAMINANTS_FASTA } from '../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { GENERATE_REFERENCE_INDEX } from '../subworkflows/generate_index.nf'
include { GET_TRANSCRIPT_INFO } from '../modules/local/transcript_info.nf'
include { GET_TRANSCRIPT_FASTA } from '../modules/local/transcript_info.nf'
include { GET_TRANSCRIPT_INFO as GET_PSITE_TRANSCRIPT_INFO } from '../modules/local/transcript_info.nf'
include { GET_TRANSCRIPT_FASTA as GET_PSITE_TRANSCRIPT_FASTA } from '../modules/local/transcript_info.nf'


workflow PREPARE_RIBOSEQ_REFERENCE {

    take:
    genome_fasta
    genome_gtf
    contaminants_fasta

    main:
    ch_versions = Channel.empty()

    // Prepare annotation: unzip annotation and genome files if necessary
    ch_genome_fasta = Channel.empty()
    if (genome_fasta.toString().endsWith('.gz')) {
        ch_genome_fasta = GUNZIP_FASTA ( [ [:], genome_fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_genome_fasta = Channel.from( [ [ [:], genome_fasta ] ] )
    }


    ch_genome_gtf = Channel.empty()
    if (genome_gtf.toString().endsWith('.gz')) {
        ch_genome_gtf = GUNZIP_GTF ( [ [:], genome_gtf ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

        // println "DEBUG: genome_gtf type: ${genome_gtf.getClass()}"
        // println "DEBUG: params.gtf type: ${params.gtf.getClass()}"

    } else {

        ch_genome_gtf = Channel.from( [ [ [:], genome_gtf ] ] )

    }

    // println "ch_genome_gtf class: ${ch_genome_gtf.getClass()}"

    ch_contaminants_fasta = Channel.empty()
    if (!params.skip_premap) {
        if (contaminants_fasta.toString().endsWith('.gz')) {
            ch_contaminants_fasta = GUNZIP_CONTAMINANTS_FASTA ( [ [:], contaminants_fasta ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_CONTAMINANTS_FASTA.out.versions)
        } else {
            // ch_contaminants_fasta = file(params.contaminants_fasta)
            ch_contaminants_fasta = Channel.from( [ [ [:], contaminants_fasta ] ] )
        }
    }

    // Prepare fasta index
    SAMTOOLS_FAIDX(ch_genome_fasta)
    ch_genome_fai = SAMTOOLS_FAIDX.out.fai.map{ it[1] }

    // Prepare annotation: create index for alignment
    GENERATE_REFERENCE_INDEX(ch_contaminants_fasta, ch_genome_fasta, ch_genome_gtf)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions, GENERATE_REFERENCE_INDEX.out.versions)

    if (!params.transcript_info) {
        GET_TRANSCRIPT_INFO(ch_genome_gtf.map{ it[1] }, 'representative')
        ch_transcript_info = GET_TRANSCRIPT_INFO.out.transcript_info
        ch_transcript_info_gtf = GET_TRANSCRIPT_INFO.out.transcripts_gtf

        GET_TRANSCRIPT_FASTA(ch_genome_fasta.map{ it[1] }, ch_genome_fai, ch_transcript_info_gtf)
        ch_versions = ch_versions.mix(GET_TRANSCRIPT_INFO.out.versions, GET_TRANSCRIPT_FASTA.out.versions)
        ch_transcript_info_fa = GET_TRANSCRIPT_FASTA.out.transcripts_fa

    } else {

        ch_transcript_info = Channel.fromPath(params.transcript_info, checkIfExists: true)
        ch_transcript_info_fa = !params.skip_psite && params.psite_transcripts == 'all_cds' ? Channel.empty() : Channel.fromPath(params.transcript_fasta, checkIfExists: true)
    }

    ch_psite_transcript_info = ch_transcript_info
    ch_psite_transcript_fa = ch_transcript_info_fa
    if (!params.skip_psite && params.psite_transcripts == 'all_cds') {
        GET_PSITE_TRANSCRIPT_INFO(ch_genome_gtf.map{ it[1] }, 'all_cds')
        GET_PSITE_TRANSCRIPT_FASTA(ch_genome_fasta.map{ it[1] }, ch_genome_fai, GET_PSITE_TRANSCRIPT_INFO.out.transcripts_gtf)
        ch_psite_transcript_info = GET_PSITE_TRANSCRIPT_INFO.out.transcript_info
        ch_psite_transcript_fa = GET_PSITE_TRANSCRIPT_FASTA.out.transcripts_fa
        ch_versions = ch_versions.mix(GET_PSITE_TRANSCRIPT_INFO.out.versions, GET_PSITE_TRANSCRIPT_FASTA.out.versions)
    }

    emit:

    contaminants_bowtie2_index = GENERATE_REFERENCE_INDEX.out.contaminants_bowtie2_index
    genome_star_index = GENERATE_REFERENCE_INDEX.out.genome_star_index
    genome_gtf = ch_genome_gtf
    genome_fasta = ch_genome_fasta
    genome_fai = ch_genome_fai
    contaminants_fasta = ch_contaminants_fasta
    transcript_info = ch_transcript_info
    // transcript_info_gtf = ch_transcript_info_gtf
    transcript_info_fa = ch_transcript_info_fa
    psite_transcript_info = ch_psite_transcript_info
    psite_transcript_fa = ch_psite_transcript_fa
    versions = ch_versions

}
