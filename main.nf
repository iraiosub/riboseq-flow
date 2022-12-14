#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

ch_input = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, file(row.fastq, checkIfExists: true) ] }
            // .view()

ch_genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.fai, checkIfExists: true)
ch_genome_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
ch_smallrna_fasta = Channel.fromPath(params.smallrna_genome, checkIfExists: true)
ch_star_index = Channel.fromPath(params.star_index, checkIfExists: true)


include { GENERATE_REFERENCE_INDEX } from './workflows/generate_reference.nf'
include { PREMAP } from './modules/premap.nf'
include { MAP } from './modules/map.nf'
include { DEDUPLICATE } from './workflows/dedup.nf'

workflow {

    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta)

    PREMAP(ch_input, GENERATE_REFERENCE_INDEX.out.smallrna_bowtie2_index)

    MAP(PREMAP.out.unmapped, ch_star_index)

    DEDUPLICATE(MAP.out.aligned_genome, MAP.out.aligned_transcriptome)

}


