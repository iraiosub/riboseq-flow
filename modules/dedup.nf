#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DEDUPLICATE_GENOME {

    tag "${sample_id}"
    label "mid_memory"

    conda '/camp/home/rebselj/.conda/envs/riboseq_env'

    publishDir "${params.outdir}/deduplicated_genome", pattern: "*.sorted.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_genome", pattern: "*.bai", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_genome", pattern: "*.bed.gz", mode: 'copy', overwrite: true


    input:
    tuple val(sample_id), path(aligned_genome)

    output:
    tuple val(sample_id), path("*.sorted.bam"), emit: dedup_genome_bam
    tuple val(sample_id), path("*.bai"), emit: dedup_genome_bai
    tuple val(sample_id), path("*.bed.gz"), emit: dedup_genome_bed


    script:
    """
    samtools view -q 20 -h $aligned_genome > ${sample_id}.um.bam
    umi_tools dedup --umi-separator 'rbc:' -I ${sample_id}.um.bam -S ${sample_id}.unsorted.bam
    samtools sort -@ ${task.cpus} ${sample_id}.unsorted.bam > ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam > ${sample_id}.bai
    rm ${sample_id}.unsorted.bam        
    bam2bed < ${sample_id}.sorted.bam | cut -f1-3,6 -d$'\t' > ${sample_id}.bed
    gzip ${sample_id}.bed
    """
//samtools: #-q 20 is probably unnecessary as we don't allow multimapping reads.

}

process DEDUPLICATE_TRANSCRIPTOME {

    tag "${sample_id}"
    label "mid_memory"

    conda '/camp/home/rebselj/.conda/envs/riboseq_env'

    publishDir "${params.outdir}/deduplicated_transcriptome", pattern: "*.sorted.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_transcriptome", pattern: "*.bai", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_transcriptome", pattern: "*.bed.gz", mode: 'copy', overwrite: true


    input:
    tuple val(sample_id), path(aligned_transcriptome)

    output:
    tuple val(sample_id), path("*.sorted.bam"), emit: dedup_transcriptome_bam
    tuple val(sample_id), path("*.bai"), emit: dedup_transcriptome_bai
    tuple val(sample_id), path("*.bed.gz"), emit: dedup_transcriptome_bed

    script:
    """
    samtools sort -@ ${task.cpus} $aligned_transcriptome > ${sample_id}.temp.bam
    samtools index ${sample_id}.temp.bam > ${sample_id}.temp.bai
    echo umitools
    umi_tools dedup --umi-separator 'rbc:' -I ${sample_id}.temp.bam -S ${sample_id}.unsorted.bam
    samtools sort -@ ${task.cpus} ${sample_id}.unsorted.bam > ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam > ${sample_id}.bai
    rm ${sample_id}.temp.bam
    rm ${sample_id}.temp.bai
    rm ${sample_id}.unsorted.bam
    bam2bed < ${sample_id}.sorted.bam | cut -f1-3,6 -d$'\t' > ${sample_id}.bed
    gzip ${sample_id}.bed
    """

}

workflow DEDUPLICATION {

    take:
    aligned_genome
    aligned_transcriptome

    main:

    //Deduplication of sequences aligned to genome
    DEDUPLICATE_GENOME(
        aligned_genome
    )

    //Deduplication of sequences aligned to transcriptome
    DEDUPLICATE_TRANSCRIPTOME(
        aligned_transcriptome
    )

    emit:
    //Genome dedup
    dedup_genome_bam = DEDUPLICATE_GENOME.out.bam
    dedup_genome_bai = DEDUPLICATE_GENOME.out.bai
    dedup_genome_bed = DEDUPLICATE_GENOME.out.bed
    //Transcriptome dedup
    dedup_transcriptome_bam = DEDUPLICATE_TRANSCRIPTOME.out.bam
    dedup_transcriptome_bai = DEDUPLICATE_TRANSCRIPTOME.out.bai
    dedup_transcriptome_bed = DEDUPLICATE_TRANSCRIPTOME.out.bed


}