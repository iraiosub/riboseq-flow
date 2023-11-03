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

    script:
    """
    get_transcript_info.R -g $gtf -o 'Homo sapiens'
    
    """
}
