#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process GENE_COUNTS_FEATURECOUNTS {
 
    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::subread=2.0.3'

    publishDir "${params.outdir}/feature_counts", pattern: "*.featureCounts*", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(sample_id), path("*.featureCounts.txt")        , emit: counts
    tuple val(sample_id), path("*.featureCounts.txt.summary"), emit: summary

    script:
    
    """
    featureCounts -a $gtf --largestOverlap -s 1 -T ${task.cpus} -o ${sample_id}.featureCounts.txt $bam
    
    """

}



