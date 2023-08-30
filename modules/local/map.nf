#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAP {
 
    tag "${sample_id}"
    label 'process_high'

    // conda 'bioconda::star=2.7.10a bioconda::samtools=1.16.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    publishDir "${params.outdir}/mapped", pattern: "*.Aligned.sortedByCoord.out.ba*", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/mapped", pattern: "*.Aligned.toTranscriptome.sorted.out.ba*", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/mapped", pattern: "*.Log.final.out", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(unmapped)
    path(genome_star_index)

    output:
    tuple val(sample_id), path("*.Aligned.sortedByCoord.out.bam"), path("*.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bam
    tuple val(sample_id), path("*.Aligned.toTranscriptome.sorted.out.bam"), path("*.Aligned.toTranscriptome.sorted.out.bam.bai"), emit: transcriptome_bam
    tuple val(sample_id), path("*.Log.final.out"), emit: log

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

//  if a read maps to the exon that is common to 15 transcripts, it will be reported as genomic unique mapper in the Aligned.out.bam and as 15-mapper in the Aligned.toTranscriptome.bam