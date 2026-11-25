#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GET_TRANSCRIPT_INFO {
    
    tag "$gtf"
    label 'process_single'

    container 'iraiosub/nf-riboseq@sha256:0c0c33861d533653fb8fe6091f247e302d68c1e9b55ee5f5faf46979a9750a4b'

    publishDir "${params.outdir}/annotation", pattern: "*.{tsv,gtf}", mode: 'copy', overwrite: true

    input:
        path(gtf)
        val(transcripts)

    output:
        path("*.transcript_info.tsv"), emit: transcript_info
        path("*_transcripts.gtf"), emit: transcripts_gtf
        path("versions.yml"), emit: versions

    script:

    def organism = params.org_name ?: 'Homo sapiens'
    def transcript_arg = transcripts == 'all_cds' ? '--transcripts all_cds' : ''
    
    """
    get_transcript_info.R -g $gtf -o "$organism" $transcript_arg

    software_versions.R "${task.process}" optparse GenomicFeatures rtracklayer data.table tidyverse dplyr
    
    """
}


process GET_TRANSCRIPT_FASTA {
    
    tag "$gtf"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'quay.io/biocontainers/gffread:0.12.1--h8b12597_0' }"

    publishDir "${params.outdir}/annotation", pattern: "*.fa", mode: 'copy', overwrite: true

    input:
        path(fasta)
        path(fai)
        path(gtf)
        

    output:
        path("*_transcripts.fa"), emit: transcripts_fa
        path("versions.yml"), emit: versions

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
