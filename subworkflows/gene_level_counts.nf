#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include { GENE_COUNTS_FEATURECOUNTS } from '../modules/local/featurecounts.nf'
include { MERGE_FEATURECOUNTS } from '../modules/local/featurecounts.nf'


workflow GET_GENE_LEVEL_COUNTS {

    take:
    bam
    genome_gtf


    main:

    // Get gene-level counts from BAM alignments using featureCounts
    GENE_COUNTS_FEATURECOUNTS(bam, genome_gtf.first())
    MERGE_FEATURECOUNTS(GENE_COUNTS_FEATURECOUNTS.out.counts.map { [ it[1] ] }.collect())

    emit:

    merged_counts_table = MERGE_FEATURECOUNTS.out.counts


}

