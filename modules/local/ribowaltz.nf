#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process IDENTIFY_PSITES {

    tag "${sample_id}"
    label 'process_high'

    // conda "bioconda::ribowaltz=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:1.2.0--r42hdfd78af_1' :
        'quay.io/biocontainers/ribowaltz:1.2.0--r42hdfd78af_1' }"

    publishDir { params.psite_per_sample ? "${params.outdir}/psites/per_sample/${sample_id}" : "${params.outdir}/psites" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam_list)
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
    path("*coverage*.tsv.gz"), emit: coverage_tables, optional: true
    path("ribowaltz_qc/*.pdf"), emit: ribowaltz_qc, optional: true

    script:

    length_range = params.psite_length_range ?: params.expected_length
    periodicity_threshold = params.periodicity_threshold
    method = params.psite_method
    def genomic_coverage = params.genomic_psite_coverage || params.genomic_codon_coverage

    // identify_p_sites.R -b $bam_folder -g $gtf -f $fasta -l $length_range --qc --method --periodicity

        """

        INPUT=`echo $bam_list | sed 's/ /,/g'`

        Rscript --vanilla ${workflow.projectDir}/bin/identify_psites.R \$INPUT $gtf $fasta $length_range $periodicity_threshold $method ${params.exclude_start} ${params.exclude_end} $transcript_info ${params.psite_per_sample} $genomic_coverage

        """

}


process MERGE_PSITE_TABLES {

    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:1.2.0--r42hdfd78af_1' :
        'quay.io/biocontainers/ribowaltz:1.2.0--r42hdfd78af_1' }"

    publishDir "${params.outdir}/psites", mode: 'copy', overwrite: true

    input:
    path(tables, stageAs: 'tables/??/*')

    output:
    path("psite_offset.tsv.gz"), emit: psite_offset
    path("codon_coverage_rpf.tsv.gz"), emit: codon_coverage_rpf
    path("codon_coverage_psite.tsv.gz"), emit: codon_coverage_psite
    path("cds_coverage_psite.tsv.gz"), emit: cds_coverage
    path("*nt_coverage_psite.tsv.gz"), emit: cds_window_coverage

    script:

        """
        Rscript --vanilla ${workflow.projectDir}/bin/merge_psite_tables.R tables
        """

}



process PREPARE_PSITE_TRACKS {

    label 'process_medium'

    container 'iraiosub/nf-riboseq@sha256:0c0c33861d533653fb8fe6091f247e302d68c1e9b55ee5f5faf46979a9750a4b'

    publishDir "${params.outdir}/annotation", pattern: "*.{rds,sizes}", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/pipeline_info/psite_tracks", pattern: "versions.yml", saveAs: { 'annotation.versions.yml' }, mode: 'copy', overwrite: true

    input:
    path(gtf)
    path(fai)
    path(transcript_info)

    output:
    path("psite_exons.rds"), emit: annotation
    path("genome.chrom.sizes"), emit: sizes
    path("versions.yml"), emit: versions

    script:
        """
        Rscript ${workflow.projectDir}/bin/prepare_psite_tracks.R -g $gtf -f $fai -t $transcript_info

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GenomicFeatures: \$(Rscript -e "cat(as.character(packageVersion('GenomicFeatures')))")
        END_VERSIONS
        """
}

process GET_PSITE_TRACKS {

    tag "${meta.id}"
    label 'process_medium'

    container 'iraiosub/nf-riboseq@sha256:0c0c33861d533653fb8fe6091f247e302d68c1e9b55ee5f5faf46979a9750a4b'

    publishDir "${params.outdir}/coverage_tracks/psite", pattern: "*.{bed.gz,mapping_stats.tsv}", mode: 'copy', overwrite: true
    publishDir "${params.outdir}/pipeline_info/psite_tracks", pattern: "versions.yml", saveAs: { 'mapping.versions.yml' }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(psite_table)
    path(annotation)
    path(sizes)

    output:
    tuple val(meta), path("*.codon_coverage_psite.bed.gz"), emit: psite_bed
    tuple val(meta), path("*.psites.mapping_stats.tsv"), emit: mapping_stats
    tuple val(meta), path("*.forward.unsorted.bedgraph"), emit: forward, optional: true
    tuple val(meta), path("*.reverse.unsorted.bedgraph"), emit: reverse, optional: true
    path("versions.yml"), emit: versions

    script:
        """
        Rscript ${workflow.projectDir}/bin/get_psite_tracks.R -p $psite_table -a $annotation -s $sizes -o ${meta.id}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GenomicFeatures: \$(Rscript -e "cat(as.character(packageVersion('GenomicFeatures')))")
            rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        END_VERSIONS
        """
}

process RUST_QC {

label 'process_single'

    container 'iraiosub/nf-riboseq-qc@sha256:719e18799ff01b3071cb2187fcae78efaafe53faa604ca4a0bf224663f0cefba'

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
