#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process PREMAP {

    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::bowtie2=2.5.0 bioconda::samtools=1.16.1'

    publishDir "${params.outdir}/premapped", pattern: "*.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/premapped", pattern: "*.bam.seqs.gz", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/premapped", pattern: "*.unmapped.fastq.gz", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/premapped", pattern: "*.premap.log", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)
    path(smallrna_index)
    

    output:
    tuple val(sample_id), path("*.bam"), emit: premapped_bam
    tuple val(sample_id), path("*.seqs.gz"), emit: seqs
    tuple val(sample_id), path("*unmapped.fastq.gz"), emit: unmapped
    path("*.premap.log"), emit: log

    script:

    args = params.bowtie2_args

    
    """
    bowtie2 $args \
        -U $reads \
        -p ${task.cpus} \
        -x ${smallrna_index[0].simpleName}  \
        --un-gz ${sample_id}.unmapped.fastq.gz \
        -S ${sample_id}.sam \
        2> ${sample_id}.premap.log

    samtools sort -o ${sample_id}.bam -O bam -@ ${task.cpus} ${sample_id}.sam
    samtools index ${sample_id}.bam
    samtools view ${sample_id}.bam | cut -f3-5,10 | gzip > ${sample_id}.bam.seqs.gz
    rm ${sample_id}.sam

    """
}