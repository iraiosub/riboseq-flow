#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DEDUPLICATE_GENOME {

    tag "${sample_id}"
    label 'process_medium'

    // conda '/camp/home/iosubi/miniconda3/envs/riboseq_nf_env'
    conda 'bioconda::umi_tools=1.1.2 conda bioconda::samtools=1.16.1 bioconda::bedtools=2.30.0'

    publishDir "${params.outdir}/deduplicated_genome", pattern: "*.dedup.sorted.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_genome", pattern: "*.dedup.bai", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_genome", pattern: "*.dedup.bed.gz", mode: 'copy', overwrite: true


    input:
    tuple val(sample_id), path(aligned_genome)

    output:
    tuple val(sample_id), path("*.dedup.sorted.bam"), path("*dedup.sorted.bai"), emit: genome_bam
    tuple val(sample_id), path("*.dedup.bed.gz"), emit: genome_bed


    script:

    """
    samtools view -q 20 -h $aligned_genome > ${sample_id}.um.bam  # -q 20 is probably unnecessary as we don't allow multimapping reads.
    umi_tools dedup --umi-separator ${params.umi_separator} -I ${sample_id}.um.bam -S ${sample_id}.unsorted.bam
    samtools sort -@ ${task.cpus} ${sample_id}.unsorted.bam > ${sample_id}.dedup.sorted.bam
    samtools index ${sample_id}.dedup.sorted.bam > ${sample_id}.dedup.sorted.bai
    rm ${sample_id}.unsorted.bam 
    rm ${sample_id}.um.bam

    # bam2bed < ${sample_id}.dedup.sorted.bam | cut -f1-3,6 > ${sample_id}.dedup.bedops.bed
    bedtools bamtobed -i ${sample_id}.dedup.sorted.bam | bedtools sort > ${sample_id}.dedup.bed
    gzip ${sample_id}.dedup.bed
    """


}

process DEDUPLICATE_TRANSCRIPTOME {

    tag "${sample_id}"
    label 'process_medium'

    conda 'bioconda::umi_tools=1.1.2 conda bioconda::samtools=1.16.1 bioconda::bedtools=2.30.0'

    publishDir "${params.outdir}/deduplicated_transcriptome", pattern: "*.tx.dedup.sorted.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_transcriptome", pattern: "*.tx.dedup.sorted.bai", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/deduplicated_transcriptome", pattern: "*.tx.dedup.bed.gz", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(aligned_transcriptome)

    output:
    tuple val(sample_id), path("*.tx.dedup.sorted.bam"), path("*.tx.dedup.sorted.bai"), emit: transcriptome_bam
    tuple val(sample_id), path("*.tx.dedup.bed.gz"), emit: transcriptome_bed

    script:
    """
    samtools sort -@ ${task.cpus} $aligned_transcriptome > ${sample_id}.temp.bam
    samtools index ${sample_id}.temp.bam > ${sample_id}.temp.bai
    
    umi_tools dedup --umi-separator 'rbc:' -I ${sample_id}.temp.bam -S ${sample_id}.unsorted.bam
    samtools sort -@ ${task.cpus} ${sample_id}.unsorted.bam > ${sample_id}.tx.dedup.sorted.bam
    samtools index ${sample_id}.tx.dedup.sorted.bam > ${sample_id}.tx.dedup.sorted.bai
    rm ${sample_id}.temp.bam
    rm ${sample_id}.temp.bai
    rm ${sample_id}.unsorted.bam

    # bam2bed < ${sample_id}.tx.dedup.sorted.bam | cut -f1-3,6 > ${sample_id}.tx.dedup.bedops.bed
    bedtools bamtobed -i ${sample_id}.tx.dedup.sorted.bam | bedtools sort > ${sample_id}.tx.dedup.bed
    gzip ${sample_id}.tx.dedup.bed
    """

}

// Remove duplicate reads from BAM file based on UMIs

workflow DEDUPLICATE {

    take:
    aligned_genome
    aligned_transcriptome

    main:

    // Deduplication of sequences aligned to genome
    DEDUPLICATE_GENOME(
        aligned_genome
    )

    // Deduplication of sequences aligned to transcriptome
    DEDUPLICATE_TRANSCRIPTOME(
        aligned_transcriptome
    )

    emit:

    // Genome dedup
    dedup_genome_bam = DEDUPLICATE_GENOME.out.genome_bam
    dedup_genome_bed = DEDUPLICATE_GENOME.out.genome_bed

    // Transcriptome dedup
    dedup_transcriptome_bam = DEDUPLICATE_TRANSCRIPTOME.out.transcriptome_bam
    dedup_transcriptome_bed = DEDUPLICATE_TRANSCRIPTOME.out.transcriptome_bed

}


// // Co-ordinate sort, index and run stats on transcriptome BAM
//             BAM_SORT_SAMTOOLS (
//                 ch_transcriptome_bam
//             )

// 
//             DEDUP_UMI_UMITOOLS_TRANSCRIPTOME (
//                 ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
//                 params.umitools_dedup_stats
//             )
