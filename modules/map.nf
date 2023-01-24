#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAP {
 
    tag "${sample_id}"
    label 'process_high'

    conda 'bioconda::star=2.7.10a bioconda::samtools=1.16.1'

    publishDir "${params.outdir}/mapped", pattern: "*.Aligned.sortedByCoord.out.ba*", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/mapped", pattern: "*.Aligned.toTranscriptome.sorted.out.ba*", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/mapped", pattern: "*.Log.final.out", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(unmapped)
    path(genome_star_index)

    output:
    tuple val(sample_id), path("*.Aligned.sortedByCoord.out.bam"), path("*.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bam
    tuple val(sample_id), path("*.Aligned.toTranscriptome.sorted.out.bam"), path("*.Aligned.toTranscriptome.sorted.out.bam.bai"), emit: transcriptome_bam
    path("*.Log.final.out"), emit: log

    script:
    
    args = params.star_args
    
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
    --quantTranscriptomeBan Singleend \
    $args

    samtools index -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam
    samtools sort -@ ${task.cpus} ${sample_id}.Aligned.toTranscriptome.out.bam > ${sample_id}.Aligned.toTranscriptome.sorted.out.bam
    samtools index -@ ${task.cpus} ${sample_id}.Aligned.toTranscriptome.sorted.out.bam
    
    """

}