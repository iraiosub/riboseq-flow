#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_TRANSCRIPT_INFO {
    
    tag "$gtf"
    label 'process_single'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env2.yml'
    container 'iraiosub/nf-riboseq:latest'

    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
        path(gtf)

    output:
        path("*.longest_cds.transcript_info.tsv"), emit: transcript_info
        path("*.longest_cds_transcripts.gtf"), emit: transcripts_gtf

    script:
    """
    get_transcript_info.R -g $gtf -o 'Homo sapiens'
    
    """
}


process GET_TRANSCRIPT_FASTA {
    
    tag "$gtf"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
        path(gtf)
        path(fasta)

    output:
        path("*.longest_cds_transcripts.fa"), emit: transcripts_fa

    script:

    def prefix = "${gtf.baseName}"
    
    """
    gffread -w ${prefix}.longest_cds_transcripts.fa -g $fasta $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    
    """
    
}