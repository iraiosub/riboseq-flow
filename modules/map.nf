#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAP {
 
    tag "${sample_id}"
    label 'process_high'

    conda '/camp/home/iosubi/miniconda3/envs/riboseq_nf_env'

    publishDir "${params.outdir}/mapped", pattern: "*.Aligned.sortedByCoord.out.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/mapped", pattern: "*.Aligned.toTranscriptome.out.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/mapped", pattern: "*.Log.final.out", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(unmapped)
    path(genome_star_index)

    output:
    tuple val(sample_id), path("*.Aligned.sortedByCoord.out.bam"), emit: aligned_genome
    tuple val(sample_id), path("*.Aligned.toTranscriptome.out.bam"), emit: aligned_transcriptome
    path("*.Log.final.out"), emit: log

    script:
    
    // map_params = params.star_args
    
    """
    STAR --runThreadN ${task.cpus} \
    --genomeDir $genome_star_index \
    --readFilesIn $unmapped \
    --readFilesCommand zcat \
    --outFileNamePrefix ${sample_id}. \
    --outSAMtype BAM SortedByCoordinate \
    --seedSearchStartLmax 15 \
    --outReadsUnmapped Fastx \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNoverReadLmax 0.08 \
    --alignEndsType EndToEnd \
    --quantMode TranscriptomeSAM \
    --outSAMattributes Standard \
    --quantTranscriptomeBan Singleend
    
    """

}