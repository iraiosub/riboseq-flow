#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process GENE_LEVEL_COUNTS {
 
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



process GENETYPE_COUNTS {
 
    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::subread=2.0.3'

    publishDir "${params.outdir}/feature_counts", pattern: "*gene_type.featureCounts*", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(sample_id), path("*.gene_type.featureCounts.txt")        , emit: counts
    tuple val(sample_id), path("*.gene_type.featureCounts.txt.summary"), emit: summary

    script:
    
    """
    featureCounts -f -g "gene_type" -s 1 -a $gtf --largestOverlap -T ${task.cpus} -o ${sample_id}.gene_type.featureCounts.txt $bam
    
    """

}