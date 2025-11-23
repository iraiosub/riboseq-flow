#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process IDENTIFY_PSITES {

    label 'process_high'

    // conda "bioconda::ribowaltz=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:1.2.0--r42hdfd78af_1' :
        'quay.io/biocontainers/ribowaltz:1.2.0--r42hdfd78af_1' }"

    publishDir "${params.outdir}/psites", mode: 'copy', overwrite: true

    input:
    path(bam_list)
    path(gtf)
    path(fasta)
    path(transcript_info)

    output:
    path("psite_offset.tsv.gz"), emit: psite_offset, optional: true
    path("offset_plot/*"), emit: offset_plot, optional: true
    path("*.psite.tsv.gz"), emit: psites, optional: true
    path("codon_coverage_rpf.tsv.gz"), emit: codon_coverage_rpf, optional: true
    path("codon_coverage_psite.tsv.gz"), emit: codon_coverage_psite, optional: true
    path("cds_coverage_psite.tsv.gz"), emit: cds_coverage, optional: true
    path("*nt_coverage_psite.tsv.gz"), emit: cds_window_coverage, optional: true
    path("ribowaltz_qc/*.pdf"), emit: ribowaltz_qc, optional: true

    script:

    length_range = params.expected_length
    periodicity_threshold = params.periodicity_threshold
    method = params.psite_method

    // identify_p_sites.R -b $bam_folder -g $gtf -f $fasta -l $length_range --qc --method --periodicity

        """

        INPUT=`echo $bam_list | sed 's/ /,/g'`

        Rscript --vanilla ${workflow.projectDir}/bin/identify_psites.R \$INPUT $gtf $fasta $length_range $periodicity_threshold $method ${params.exclude_start} ${params.exclude_end} $transcript_info

        """

}



process GET_PSITE_TRACKS {

    label 'process_high'

    container 'iraiosub/nf-riboseq:latest'

    publishDir "${params.outdir}/coverage_tracks/psite", pattern: "*.psites.bed.gz", mode: 'copy', overwrite: true
    // publishDir "${params.outdir}/coverage_tracks/psite", pattern: "*.bigWig", mode: 'copy', overwrite: true

    input:
    path(psite_table)
    path(gtf)
    path(fai)

    output:
    path("*.codon_coverage_psite.bed.gz"), emit: psite_bed
    // path("*.bigwig"), emit: psite_bigwig

    script:

        """

        get_psite_tracks.R -p $psite_table -g $gtf -f $fai

        """

}

process RUST_QC {

label 'process_single'

    container 'iraiosub/nf-riboseq-qc:latest'

    publishDir "${params.outdir}/riboseq_qc/rust_analysis/", pattern: "*.rust_analysis.pdf", mode: 'copy', overwrite: true
    
    input:
    path(transcript_fasta)
    path(transcript_info)
    path (psites)


    output:
    path("*.rust_analysis.pdf"), emit: rust_analysis

    script:

        """
        make_rust_ratio_plots.R -f $transcript_fasta -t $transcript_info -p $psites
            
        """

}