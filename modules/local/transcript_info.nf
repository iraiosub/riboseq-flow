#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_TRANSCRIPT_INFO {
    
    tag "$gtf"
    label 'process_single'

    container 'iraiosub/nf-riboseq:latest'

    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
        path(gtf)

    output:
        path("*.longest_cds.transcript_info.tsv"), emit: transcript_info
        path("*.longest_cds_transcripts.gtf"), emit: transcripts_gtf

    script:

    def organism = params.org_name ?: 'Homo sapiens'
    
    """
    get_transcript_info.R -g $gtf -o "$organism"
    
    """
}


process GET_TRANSCRIPT_FASTA {
    
    tag "$gtf"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'quay.io/biocontainers/gffread:0.12.1--h8b12597_0' }"

    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true

    input:
        path(fasta)
        path(fai)
        path(gtf)
        

    output:
        path("*.longest_cds_transcripts.fa"), emit: transcripts_fa

    script:

    def prefix = "${gtf.baseName}"
    
    """
    gffread -w ${prefix}_gffread.fa -g $fasta $gtf

    sed -e "s/^>\\([^ ]*\\) .*/>\\1/" ${prefix}_gffread.fa > ${prefix}.fa


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    
    """
    
}
