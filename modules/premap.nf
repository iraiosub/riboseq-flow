#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process PREMAP {

    tag "${sample_id}"
    label "mid_memory"


    //container 'quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0'
    // add conda!
    conda '/camp/home/rebselj/.conda/envs/riboseq_env'

    publishDir "${params.outdir}/premap", pattern: "*.bam", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/premap", pattern: "*.bam.seqs.gz", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/premap", pattern: "*.unmapped.fastq.gz", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/premap", pattern: "*.premap.log", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)
    path(smallrna_index)
    

    output:
    tuple val(sample_id), path("*.bam"), emit: premapped_bam
    tuple val(sample_id), path("*.seqs.gz"), emit: seqs
    tuple val(sample_id), path("*unmapped.fastq.gz"), emit: unmapped
    path("*.premap.log"), emit: log

    script:

    map_params = params.bowtie2_args

    //  cmd = "STAR $args && samtools index -@ ${task.cpus} ${sample_id}.Aligned.sortedByCoord.out.bam"
    
    """
    bowtie2 $map_params \
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

//32g,6h