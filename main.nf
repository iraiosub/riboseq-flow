#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

/* 
PREPARE INPUTS
*/

// Check mandatory parameters
if (!params.input) {
    exit 1, 'Input samplesheet not specified!'
} else {
    ch_input = Channel
            .fromPath( params.input )
            .splitCsv(header:true)
            .map { row -> [ row.sample, file(row.fastq, checkIfExists: true) ] }
            // .view()
}


// Genome variables

if(params.org) {

    params.fasta = params.genomes[ params.org ].fasta
    params.fai = params.genomes[ params.org ].fai
    params.gtf = params.genomes[ params.org ].gtf
    params.star_index = params.genomes[ params.org ].star_index
    params.smallrna_fasta = params.genomes[ params.org ].smallrna_fasta

}  else {

    if(!params.fasta ) { exit 1, "--fasta is not specified." } 
    if(!params.gtf ) { exit 1, "--gtf is not specified." } 
    if(!params.smallrna_fasta && params.premap ) {exit 1, "--smallrna_fasta is not specified." }

}


// Create channels for static files
ch_genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.fai, checkIfExists: true)
ch_genome_gtf = Channel.fromPath(params.gtf, checkIfExists: true)

if (params.premap) {
    ch_smallrna_fasta = Channel.fromPath(params.smallrna_fasta, checkIfExists: true)
    if (ch_smallrna_fasta.isEmpty()) {exit 1, "File provided with --smallrna_fasta is empty: ${ch_smallrna_fasta.getName()}!"}
}
// ch_star_index = Channel.fromPath(params.star_index, checkIfExists: true)


include { GENERATE_REFERENCE_INDEX } from './workflows/generate_reference.nf'
include { PREPROCESS_READS } from './workflows/preprocess_reads.nf'
include { PREMAP } from './modules/premap.nf'
include { MAP } from './modules/map.nf'
include { DEDUPLICATE } from './workflows/dedup.nf'

workflow {

    // Prepare annotation
    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta, ch_genome_fasta, ch_genome_gtf)

    // Extract UMIs and/or trim adapters
    PREPROCESS_READS(ch_input)

    if (!params.skip_premap) {

        // Premap to the small RNA genome
        PREMAP(PREPROCESS_READS.out.fastq, GENERATE_REFERENCE_INDEX.out.smallrna_bowtie2_index.collect())
        MAP(PREMAP.out.unmapped, GENERATE_REFERENCE_INDEX.out.genome_star_index.collect())
    } else {
        MAP(PREPROCESS_READS.out.fastq, GENERATE_REFERENCE_INDEX.out.genome_star_index.collect())
    }
    
    if (params.with_umi) {
        // Remove duplicate reads from BAM file based on UMIs
        DEDUPLICATE(MAP.out.aligned_genome, MAP.out.aligned_transcriptome)
    }

    // Count reads from BAM alignments
}


workflow.onComplete {

    if (workflow.success) {
        log.info "-\033[0;34m[riboseq]\033[0;32m Pipeline completed successfully\033[0m-\n"
    } else {
        log.info "-\033[0;34m[riboseq]\033[1;91m Pipeline completed with errors\033[0m\n"
    }

}

