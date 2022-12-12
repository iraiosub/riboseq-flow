#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process MAP {
 
    tag "${sample_id}"
    label "high_memory"


    //container 'quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0'
    // add conda!
    conda '/camp/home/rebselj/.conda/envs/riboseq_env'

    publishDir "${params.outdir}/map", pattern: "*.Aligned.sortedByCoord.out.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/map", pattern: "*.Aligned.toTranscriptome.out.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/map", pattern: "*.Log.final.out", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(unmapped)
    path(genome_star_index)

    output:
    tuple val(sample_id), path("*.Aligned.sortedByCoord.out.bam"), emit: aligned_genome
    tuple val(sample_id), path("*.Aligned.toTranscriptome.out.bam"), emit: aligned_transcriptome
    path("*.Log.final.out"), emit: log

    script:
    
    map_params = params.star_args
    """
    STAR --runThreadN ${task.cpus} \
    --genomeDir $genome_star_index \
    --readFilesIn $unmaped \
    --readFilesCommand zcat \
    $map_params
    """
    
    //--outFileNamePrefix {params.outprefix} \
    // mv {params.log} {output.logfolder}

}