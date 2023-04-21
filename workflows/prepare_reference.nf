#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { GENERATE_REFERENCE_INDEX } from '../workflows/generate_index.nf'
include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_GTF } from '../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_SMALLRNA_FASTA } from '../modules/nf-core/gunzip/main.nf'
include { GET_TRANSCRIPT_INFO } from '../modules/local/reference.nf'


workflow {

    take:
    genome_fasta
    genome_gtf
    smallrna_fasta


    // Prepare annotation: unzip annotation and genome files if necessary
    if (genome_fasta.endsWith('.gz')) {
        ch_genome_fasta = GUNZIP_FASTA ( [ [:], genome_fasta ] ).gunzip.map { it[1] }
    } else {
        ch_genome_fasta = genome_fasta
    }


    if (genome_gtf.endsWith('.gz')) {
        ch_genome_gtf = GUNZIP_GTF ( [ [:], genome_gtf ] ).gunzip.map { it[1] }
    } else {

        ch_genome_gtf = genome_gtf
    }

    if (!params.skip_premap && smallrna_fasta.endsWith('.gz')) {
        ch_smallrna_fasta = GUNZIP_SMALLRNA_FASTA ( [ [:], smallrna_fasta ] ).gunzip.map { it[1] }
    } else {
        // ch_smallrna_fasta = file(params.smallrna_fasta)
        ch_smallrna_fasta = smallrna_fasta
    }

    // Prepare annotation: create index for alignment
    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta, ch_genome_fasta, ch_genome_gtf)

    if (!params.skip_qc) {
        GET_TRANSCRIPT_INFO(ch_genome_gtf)
    }



}