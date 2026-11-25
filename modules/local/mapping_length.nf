#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAPPING_LENGTH_ANALYSIS {
 
    tag "${sample_id}"
    label 'process_single'

    container 'iraiosub/nf-riboseq-qc@sha256:719e18799ff01b3071cb2187fcae78efaafe53faa604ca4a0bf224663f0cefba'

    publishDir "${params.outdir}/riboseq_qc/mapping_length_analysis", pattern: "*.csv", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(fastq) // before deduplication; need to use join in the workflow to make sure samples arevmatched by ID! .bai not required, but necessary to keep cardinality

    output:
    tuple val(sample_id), path("*.csv"), emit: length_analysis
    path("versions.yml"), emit: versions

    script:
    
    def analysis_type = params.length_analysis_type
    def version_command = """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        mapping-length-analysis: \$(python -c 'from pkg_resources import get_distribution; print(get_distribution("mapping-length-analysis").version)')
    END_VERSIONS
    """

    if (analysis_type == 'after_dedup') 
        """
        mapping_length_analysis -b $bam -o ${sample_id}.${analysis_type}.csv
        $version_command
        """
    else
        """
        mapping_length_analysis -f $fastq -b $bam -o ${sample_id}.${analysis_type}.csv
        $version_command
        """
}