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


// Check if genome exists in the config file
if (params.org  && !params.genomes.containsKey(params.org)) {
    exit 1, "The provided genome '${params.org}' is not available. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Genome variables

if(params.org) {

    params.fasta = params.genomes[ params.org ].fasta
    params.fai = params.genomes[ params.org ].fai
    params.gtf = params.genomes[ params.org ].gtf
    params.star_index = params.genomes[ params.org ].star_index
    params.smallrna_fasta = params.genomes[ params.org ].smallrna_fasta

}  else {

    if(!params.fasta ) { exit 1, '--fasta is not specified.' } 
    if(!params.gtf ) { exit 1, '--gtf is not specified.' } 
    if(!params.smallrna_fasta && !params.skip_premap ) {exit 1, '--smallrna_fasta is not specified.' }

}


// Create channels for static files
ch_genome_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
ch_genome_fai = Channel.fromPath(params.fai, checkIfExists: true)
ch_genome_gtf = Channel.fromPath(params.gtf, checkIfExists: true)

if (!params.skip_premap) {
    ch_smallrna_fasta = Channel.fromPath(params.smallrna_fasta, checkIfExists: true)
    // if (ch_smallrna_fasta.isEmpty()) {exit 1, "File provided with --smallrna_fasta is empty: ${ch_smallrna_fasta.getName()}!"}
} else {

    // Create empty channel so GENERATE_REFERENCE_INDEX doesn't break
    ch_smallrna_fasta = Channel.value()
}


include { GENERATE_REFERENCE_INDEX } from './workflows/generate_reference.nf'
include { PREPROCESS_READS } from './workflows/preprocess_reads.nf'
include { FASTQC } from './modules/fastqc.nf'
include { PREMAP } from './modules/premap.nf'
include { MAP } from './modules/map.nf'
include { DEDUPLICATE } from './workflows/dedup.nf'
include { MAPPING_LENGTH_ANALYSES } from './workflows/mapping_length_analyses.nf'

workflow {

    // Prepare annotation
    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta, ch_genome_fasta, ch_genome_gtf)

    // Extract UMIs and/or trim adapters
    PREPROCESS_READS(ch_input)
    FASTQC(PREPROCESS_READS.out.fastq)

    if (!params.skip_premap) {

        // Premap to the small RNA genome
        PREMAP(PREPROCESS_READS.out.fastq, GENERATE_REFERENCE_INDEX.out.smallrna_bowtie2_index.collect())
        MAP(PREMAP.out.unmapped, GENERATE_REFERENCE_INDEX.out.genome_star_index.collect())

    } else {
        MAP(PREPROCESS_READS.out.fastq, GENERATE_REFERENCE_INDEX.out.genome_star_index.collect())
    }
    
    if (params.with_umi) {
        // Remove duplicate reads from BAM file based on UMIs
        DEDUPLICATE(MAP.out.genome_bam, MAP.out.transcriptome_bam)
    }


    if (!params.skip_qc) {
    // Mapping length analysis
    MAPPING_LENGTH_ANALYSES(MAP.out.genome_bam, PREPROCESS_READS.out.fastq, DEDUPLICATE.out.dedup_genome_bam, PREMAP.out.unmapped)
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

