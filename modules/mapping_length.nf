#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAPPING_LENGTH_ANALYSIS {
 
    tag "${sample_id}"
    label 'process_high'

    conda '/camp/lab/luscomben/home/users/iosubi/projects/riboseq_nf/riboseq/env.yml'

    publishDir "${params.outdir}/mapping_length_analysis", pattern: "*.Aligned.sortedByCoord.out.ba*", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(fastq) // before deduplication; need to use join in the workflow to make sure samples arevmatched by ID!

    output:
    tuple val(sample_id), path("*.csv"), emit: length_analysis

    script:
    
    analysis_type = params.length_analysis_type

    if (analysis_type == 'after_dedup') 
        """
        python3 mapping_length_analysis.py -b $bam -o ${sample_id}.${analysis_type}.csv
        """
    else
        """
        mapping_length_analysis -f $fastq -b $bam -o ${sample_id}.${analysis_type}.csv
        """
}