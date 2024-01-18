#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAP {
 
    tag "${sample_id}"
    label 'process_high'

    // conda 'bioconda::star=2.7.10a bioconda::samtools=1.18'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

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
    --outFileNamePrefix ${sample_id}. \
    $args

    samtools index -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam

    samtools sort -@ ${task.cpus} ${sample_id}.Aligned.toTranscriptome.out.bam > ${sample_id}.Aligned.toTranscriptome.sorted.out.bam
    samtools index -@ ${task.cpus} ${sample_id}.Aligned.toTranscriptome.sorted.out.bam

    """

}

//  if a read maps to the exon that is common to 15 transcripts, it will be reported as genomic unique mapper in the Aligned.out.bam and as 15-mapper in the Aligned.toTranscriptome.bam

// Steps needed with STAR 2.7.9a bug
// samtools view -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam | cut -f1 | sort | uniq > ${sample_id}.qnames_list.txt
// samtools view -@ ${task.cpus} -N ${sample_id}.qnames_list.txt -o ${sample_id}.Aligned.toTranscriptome.filtered.out.bam ${sample_id}.Aligned.toTranscriptome.out.bam