#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process UMITOOLS_EXTRACT {
    tag "${sample_id}"
    label "process_low"

    conda 'bioconda::umi_tools=1.1.2'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.umi_extract.fastq.gz"), emit: fastq
        path("*.umi_extract.log"), emit: log

    script:
    args = " --bc-pattern=" + params.umi_pattern
    args += " --extract-method=" + params.umi_extract_method

    """
    umi_tools \
        extract \
        -I $reads \
        -S ${sample_id}.umi_extract.fastq.gz \
        $args \
        -L ${sample_id}.umi_extract.log
    """
}


process UMITOOLS_DEDUPLICATE {

    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::umi_tools=1.1.2 conda bioconda::samtools=1.16.1 bioconda::bedtools=2.30.0'

    publishDir "${params.outdir}/deduplicated", pattern: "*.dedup.sorted.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated", pattern: "*.dedup.sorted.bai", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated", pattern: "*.dedup.bed.gz", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated", pattern: "*.log", mode: 'copy', overwrite: true


    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("*.dedup.sorted.bam"), path("*dedup.sorted.bai"), emit: dedup_bam
    tuple val(sample_id), path("*.dedup.bed.gz"), emit: dedup_bed
    tuple val(sample_id), path("*.log"), emit: log


    script:

    suffix = params.dedup_mode

    """ 
    umi_tools dedup --umi-separator ${params.umi_separator} -I $bam -S ${sample_id}.${suffix}.dedup.unsorted.bam --log ${sample_id}.${suffix}.dedup.log
    samtools sort -@ ${task.cpus} ${sample_id}.${suffix}.dedup.unsorted.bam > ${sample_id}.${suffix}.dedup.sorted.bam
    samtools index ${sample_id}.${suffix}.dedup.sorted.bam > ${sample_id}.${suffix}.dedup.sorted.bai

    bedtools bamtobed -i ${sample_id}.${suffix}.dedup.sorted.bam | bedtools sort > ${sample_id}.${suffix}.dedup.bed
    gzip ${sample_id}.${suffix}.dedup.bed
    """


}


//  # samtools view -q 20 -h $aligned_genome > ${sample_id}.um.bam  # -q 20 is probably unnecessary as we don't allow multimapping reads.