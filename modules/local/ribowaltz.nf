#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2
 
process IDENTIFY_PSITES {
 
    tag "${sample_id}"
    label 'process_medium'

    // conda "bioconda::ribowaltz=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:1.2.0--r42hdfd78af_1' :
        'quay.io/biocontainers/ribowaltz:1.2.0--r42hdfd78af_1' }"

    publishDir "${params.outdir}/ribowaltz", mode: 'copy', overwrite: true
    
    input:
    path(bam_folder)
    path(gtf)
    path(fasta)

    output:
    path("psite_offset.tsv.gz"), emit: psite_offset
    // path("offset_plot"), emit: offset_plot
    path("*.psite.tsv.gz"), emit: psites
    path("codon_coverage_rpf.tsv.gz"), emit: codon_coverage_rpf
    path("codon_coverage_psite.tsv.gz"), emit: codon_coverage_psite
    path("cds_coverage_psite.tsv.gz"), emit: cds_coverage
    path("ribowaltz_qc/*.pdf"), emit: ribowaltz_qc

    script:

    length_range = params.length_range

    //     identify_p_sites.R -b $bam_folder -g $gtf -f $fasta -l $length_range

        """

        Rscript --vanilla identify_psites.R $bam_folder $gtf $fasta $length_range

        """

}
