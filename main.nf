#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

/* 
PREPARE INPUTS
*/
ch_input = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, file(row.fastq, checkIfExists: true) ] }
            // .view()


// Genome variables
params.fasta = params.genomes[ params.org ].fasta
params.fai = params.genomes[ params.org ].fai
params.gtf = params.genomes[ params.org ].gtf
params.star_genome = params.genomes[ params.org ].star_genome
params.smallrna_genome = params.genomes[ params.org ].smallrna_genome

// Create channels for static files
ch_genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.fai, checkIfExists: true)
ch_genome_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
ch_smallrna_fasta = Channel.fromPath(params.smallrna_genome, checkIfExists: true)
// ch_star_index = Channel.fromPath(params.star_index, checkIfExists: true)


include { GENERATE_REFERENCE_INDEX } from './workflows/generate_reference.nf'
include { PREMAP } from './modules/premap.nf'
include { MAP } from './modules/map.nf'

workflow {

    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta, ch_genome_fasta, ch_genome_gtf)

    if (!params.skip_premap) {
        PREMAP(ch_input, GENERATE_REFERENCE_INDEX.out.smallrna_bowtie2_index)
        MAP(PREMAP.out.unmapped, GENERATE_REFERENCE_INDEX.out.genome_star_index)
    } else {
        MAP(ch_input, GENERATE_REFERENCE_INDEX.out.genome_star_index)
    }

     DEDUPLICATE(MAP.out.aligned_genome, MAP.out.aligned_transcriptome)

}