#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process GENE_LEVEL_COUNTS {
 
    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::subread=2.0.3'

    publishDir "${params.outdir}/feature_counts", pattern: "*..rfp_counts", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    path("*.featureCounts.txt"), emit: counts

    script:
    
    """
    featureCounts -a $gtf -s 1 -T ${task.cpus} -o ${sample_id}.featureCounts.txt $bam
    
    """

}



process GENETYPE_COUNTS {
 
    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::subread=2.0.3'

    publishDir "${params.outdir}/feature_counts", pattern: "*..rfp_counts", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(meta), path("*.biotype.featureCounts.txt")        , emit: counts
    tuple val(meta), path("*.featureCounts.txt.summary"), emit: summary

    script:
    
    """
    featureCounts -f -g "gene_type" -s 1 -a $gtf  --largestOverlap -T ${task.cpus} -o ${sample_id}.biotype.featureCounts.txt $bam
    
    """

}