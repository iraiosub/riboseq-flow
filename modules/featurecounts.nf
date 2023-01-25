#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process GENE_LEVEL_COUNTS {
 
    tag "${sample_id}"
    label 'process_high'

    conda 'bioconda::subread=2.0.3'

    publishDir "${params.outdir}/feature_counts", pattern: "*..rfp_counts", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    path("*.rfp_counts"), emit: counts

    script:
    
    """
    featureCounts -i $bam -a $gtf -T ${task.cpus} -o ${sample_id}.rfp_counts
    
    """

}