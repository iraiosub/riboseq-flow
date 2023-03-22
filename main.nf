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
    ch_smallrna_fasta = Channel.empty()
}


include { GUNZIP as GUNZIP_FASTA } from './modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_GTF } from './modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_SMALLRNA_FASTA } from './modules/nf-core/gunzip/main.nf'
include { GENERATE_REFERENCE_INDEX } from './workflows/generate_reference.nf'
include { GET_TRANSCRIPT_INFO } from './modules/local/reference.nf'
include { PREPROCESS_READS } from './workflows/preprocess_reads.nf'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { PREMAP } from './modules/local/premap.nf'
include { MAP } from './modules/local/map.nf'
include { DEDUPLICATE } from './workflows/dedup.nf'
include { MAPPING_LENGTH_ANALYSES } from './workflows/mapping_length_analyses.nf'
include { RIBOSEQ_QC } from './modules/local/riboseq_qc.nf'
include { SUMMARISE_RIBOSEQ_QC } from './modules/local/riboseq_qc.nf'
include { GENE_COUNTS_FEATURECOUNTS } from './modules/local/featurecounts.nf'
include { MERGE_FEATURECOUNTS } from './modules/local/featurecounts.nf'
include { MULTIQC } from './modules/local/multiqc.nf'
include { RUN_RIBOCUTTER } from './workflows/ribocutter_analysis.nf'

workflow {

    // Prepare annotation: unzip annotation and genome files if necessary
    if (params.fasta.endsWith('.gz')) {
        ch_genome_fasta = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
    } else {
        ch_genome_fasta = ch_genome_fasta
    }


    if (params.gtf.endsWith('.gz')) {
        ch_genome_gtf = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map { it[1] }
    } else {

        ch_genome_gtf = ch_genome_gtf
    }

    if (!params.skip_premap && params.smallrna_fasta.endsWith('.gz')) {
        ch_smallrna_fasta = GUNZIP_SMALLRNA_FASTA ( [ [:], params.smallrna_fasta ] ).gunzip.map { it[1] }
    } else {
        // ch_smallrna_fasta = file(params.smallrna_fasta)
        ch_smallrna_fasta = ch_smallrna_fasta
    }

    // Prepare annotation: create index for alignment
    GENERATE_REFERENCE_INDEX(ch_smallrna_fasta, ch_genome_fasta, ch_genome_gtf)

    if (!params.skip_qc) {
        GET_TRANSCRIPT_INFO(ch_genome_gtf)
    }
    
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

        if (params.with_umi && !params.skip_premap) {

        RIBOSEQ_QC(DEDUPLICATE.out.dedup_transcriptome_bam.join(MAPPING_LENGTH_ANALYSES.out.before_dedup_length_analysis).join(MAPPING_LENGTH_ANALYSES.out.after_premap_length_analysis).join(MAPPING_LENGTH_ANALYSES.out.after_dedup_length_analysis), GET_TRANSCRIPT_INFO.out.transcript_info.collect())


        ch_merge_qc = RIBOSEQ_QC.out.qc
            .map { [ it[1] ] }
            .collect()

        SUMMARISE_RIBOSEQ_QC(ch_merge_qc)

        }
    }

    // Get gene-level counts from BAM alignments using featureCounts
    if (params.with_umi) {
        GENE_COUNTS_FEATURECOUNTS(DEDUPLICATE.out.dedup_genome_bam, ch_genome_gtf.collect())
    
    } else {
        GENE_COUNTS_FEATURECOUNTS(MAP.out.genome_bam, ch_genome_gtf.collect())

    }

    MERGE_FEATURECOUNTS(GENE_COUNTS_FEATURECOUNTS.out.counts.map { [ it[1] ] }.collect())

    // ch_logs = FASTQC.out.html.map { [ it[1] ] }.collect().mix(PREMAP.out.log.collect(), MAP.out.log.collect()).collect()
    ch_logs = FASTQC.out.html.join(FASTQC.out.zip).map { [ it[1], it[2] ] }.collect().mix(PREMAP.out.log.collect(), MAP.out.log.collect()).collect()
    MULTIQC(ch_logs)

    // Run ribocutter on trimmed but not length filtered (for ts_trimming, on trimmed but not rGrGrG-cut or length filtered)
    if (!params.skip_ribocutter) {
        RUN_RIBOCUTTER(PREPROCESS_READS.out.trimmed_fastq)
    }
    
}


workflow.onComplete {

    if (workflow.success) {
        log.info "-\033[0;34m[riboseq]\033[0;32m Pipeline completed successfully\033[0m-\n"
    } else {
        log.info "-\033[0;34m[riboseq]\033[1;91m Pipeline completed with errors\033[0m\n"
    }

}

