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
    // params.fai = params.genomes[ params.org ].fai
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
// ch_genome_fai = Channel.fromPath(params.fai, checkIfExists: true)
ch_genome_gtf = Channel.fromPath(params.gtf, checkIfExists: true)

if (!params.skip_premap) {
    ch_smallrna_fasta = Channel.fromPath(params.smallrna_fasta, checkIfExists: true)
    // if (ch_smallrna_fasta.isEmpty()) {exit 1, "File provided with --smallrna_fasta is empty: ${ch_smallrna_fasta.getName()}!"}
} else {

    // Create empty channel so GENERATE_REFERENCE_INDEX doesn't break
    ch_smallrna_fasta = Channel.empty()
}


include { PREPARE_RIBOSEQ_REFERENCE } from './workflows/prepare_reference.nf'
include { PREPROCESS_READS } from './workflows/preprocess_reads.nf'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { PREMAP } from './modules/local/premap.nf'
include { MAP } from './modules/local/map.nf'
include { DEDUPLICATE } from './workflows/dedup.nf'
include { MAPPING_LENGTH_ANALYSES } from './workflows/mapping_length_analyses.nf'
include { RIBOSEQ_QC } from './modules/local/riboseq_qc.nf'
include { SUMMARISE_RIBOSEQ_QC } from './modules/local/riboseq_qc.nf'
include { GET_GENE_LEVEL_COUNTS } from './workflows/gene_level_counts.nf'
include { IDENTIFY_PSITES } from './modules/local/ribowaltz.nf'
include { PCA } from './modules/local/featurecounts.nf'
include { MULTIQC } from './modules/local/multiqc.nf'
include { RUN_RIBOCUTTER } from './workflows/ribocutter_analysis.nf'

workflow {

    // Prepare annotation
    PREPARE_RIBOSEQ_REFERENCE(ch_genome_fasta, ch_genome_gtf, ch_smallrna_fasta)
    
    // Extract UMIs and/or trim adapters
    PREPROCESS_READS(ch_input)
    FASTQC(PREPROCESS_READS.out.fastq)

    if (!params.skip_premap) {

        // Premap to the small RNA genome
        PREMAP(PREPROCESS_READS.out.fastq, PREPARE_RIBOSEQ_REFERENCE.out.smallrna_bowtie2_index.collect())
        MAP(PREMAP.out.unmapped, PREPARE_RIBOSEQ_REFERENCE.out.genome_star_index.collect())

    } else {
        MAP(PREPROCESS_READS.out.fastq, PREPARE_RIBOSEQ_REFERENCE.out.genome_star_index.collect())
    }
    
    if (params.with_umi) {
        // Remove duplicate reads from BAM file based on UMIs
        DEDUPLICATE(MAP.out.genome_bam, MAP.out.transcriptome_bam)
    }


    if (!params.skip_qc) {
        // Mapping length analysis
        MAPPING_LENGTH_ANALYSES(MAP.out.genome_bam, PREPROCESS_READS.out.fastq, DEDUPLICATE.out.dedup_genome_bam, PREMAP.out.unmapped)

        if (params.with_umi && !params.skip_premap) {

        RIBOSEQ_QC(DEDUPLICATE.out.dedup_transcriptome_bam.join(MAPPING_LENGTH_ANALYSES.out.before_dedup_length_analysis).join(MAPPING_LENGTH_ANALYSES.out.after_premap_length_analysis).join(MAPPING_LENGTH_ANALYSES.out.after_dedup_length_analysis), PREPARE_RIBOSEQ_REFERENCE.out.transcript_info.collect())

        ch_merge_qc = RIBOSEQ_QC.out.qc
            .map { [ it[1] ] }
            .collect()

        SUMMARISE_RIBOSEQ_QC(ch_merge_qc)

        }
    }

    // Get gene-level counts from BAM alignments using featureCounts
    if (params.with_umi) {
        GET_GENE_LEVEL_COUNTS(DEDUPLICATE.out.dedup_genome_bam, ch_genome_gtf.collect())
    
    } else {
        GET_GENE_LEVEL_COUNTS(MAP.out.genome_bam, ch_genome_gtf.collect())

    }


    if (!params.skip_psite) {

         if (params.with_umi) {
            IDENTIFY_PSITES(DEDUPLICATE.out.dedup_transcriptome_bam.map { [ it[1] ] }.collect(), ch_genome_gtf.collect(), ch_genome_fasta.collect())
        
        } else {
            
            IDENTIFY_PSITES(MAP.out.transcriptome_bam.map { [ it[1] ] }.collect(), ch_genome_gtf.collect(), ch_genome_fasta.collect())

        }

    }
   

    if (!params.skip_psite) {
        PCA(GET_GENE_LEVEL_COUNTS.out.merged_counts_table, IDENTIFY_PSITES.out.cds_coverage, IDENTIFY_PSITES.out.cds_window_coverage, GET_TRANSCRIPT_INFO.out.transcript_info)
    
    } else {

        PCA(GET_GENE_LEVEL_COUNTS.out.merged_counts_table, Channel.empty(), Channel.empty(), GET_TRANSCRIPT_INFO.out.transcript_info)
    }

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

