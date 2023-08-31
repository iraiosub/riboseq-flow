#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAPPING_LENGTH_ANALYSIS {
 
    tag "${sample_id}"
    label 'process_medium'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env.yml'
    container 'iraiosub/nf-riboseq-qc:latest'

    publishDir "${params.outdir}/mapping_length_analysis", pattern: "*.csv", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(fastq) // before deduplication; need to use join in the workflow to make sure samples arevmatched by ID! .bai not required, but necessary to keep cardinality

    output:
    tuple val(sample_id), path("*.csv"), emit: length_analysis

    script:
    
    analysis_type = params.length_analysis_type

    if (analysis_type == 'after_dedup') 
        """
        mapping_length_analysis -b $bam -o ${sample_id}.${analysis_type}.csv
        """
    else
        """
        mapping_length_analysis -f $fastq -b $bam -o ${sample_id}.${analysis_type}.csv
        """
}