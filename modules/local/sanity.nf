!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process FILTER_OUT_RRNA {
 
    tag "${sample_id}"
    label 'process_medium'

    // conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'quay.io/biocontainers/subread:2.0.1--hed695b0_0' }"

    publishDir "${params.outdir}/featurecounts", pattern: "*.featureCounts*", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(sample_id), path("*.featureCounts.txt"), emit: counts
    tuple val(sample_id), path("*.featureCounts.txt.summary"), emit: summary

    script:


    if (params.strandedness == "forward") {
        
        strandedness = 1
    } else if (params.strandedness == "reverse") {
        
        strandedness = 2 
    } else if (params.strandedness == "unstranded") {

        strandedness = 0
    } else {

        error "Library strandedness must be a valid argument. Options are 'forward', 'reverse', 'unstranded'. "
    }
    
    
        """
        featureCounts -a $gtf -s $strandedness -T ${task.cpus} -o ${sample_id}.featureCounts.txt $bam

        """

}
