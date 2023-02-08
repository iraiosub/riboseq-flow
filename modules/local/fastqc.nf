#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2


process FASTQC {
    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::fastqc=0.11.9'

    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
        tuple val(sample_id), path(reads)

    output:
        path("*fastqc.{zip,html}"), emit: fastqc
    
    script:
    read_ext = reads.getName().split('\\.', 2)[1]
    read_name = reads.getName().split('\\.', 2)[0]
    new_reads = "${sample_id}_fastqc.${read_ext}"
    new_reads_simple = "${sample_id}_fastqc"
    """
    cp ${reads} ${new_reads}
    fastqc --quiet --threads ${task.cpus} ${new_reads}
    mv ${new_reads_simple}*.html ${sample_id}_fastqc.html
    mv ${new_reads_simple}*.zip ${sample_id}_fastqc.zip
    """
}