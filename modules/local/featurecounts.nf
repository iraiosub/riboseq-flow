#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process GENE_COUNTS_FEATURECOUNTS {

    tag "${sample_id}"
    label 'process_medium'

    // conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'quay.io/biocontainers/subread:2.0.1--hed695b0_0' }"

    publishDir "${params.outdir}/featurecounts", pattern: "*.featureCounts*", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(gtf)

    output:
    tuple val(sample_id), path("*.featureCounts.txt"), emit: counts
    tuple val(sample_id), path("*.featureCounts.txt.summary"), emit: summary
    path("versions.yml"), emit: versions

    script:

    def feature = params.feature ?: 'exon'

    if (params.strandedness == "forward") {

        strandedness = 1
    } else if (params.strandedness == "reverse") {

        strandedness = 2 
    } else if (params.strandedness == "unstranded") {

        strandedness = 0
    } else {

        error "Library strandedness must be a valid argument. Options are 'forward', 'reverse', 'unstranded'. "
    }


        """
        featureCounts -a $gtf -s $strandedness -t $feature -T ${task.cpus} -o ${sample_id}.featureCounts.txt $bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            featureCounts: \$(featureCounts -v 2>&1 | sed -n 's/^featureCounts v//p')
        END_VERSIONS

        """

}


process MERGE_FEATURECOUNTS {

    label 'process_single'

    container 'iraiosub/nf-riboseq-qc@sha256:719e18799ff01b3071cb2187fcae78efaafe53faa604ca4a0bf224663f0cefba'

    publishDir "${params.outdir}/featurecounts", pattern: "*.featureCounts.tsv.gz", mode: 'copy', overwrite: true

    input:
    path(featurecounts_tables)

    output:
    path("merged.featureCounts.tsv.gz"), emit: counts
    path("versions.yml"), emit: versions

    script:

        """
        INPUT=`echo $featurecounts_tables | sed 's/ /,/g'`

        get_featurecounts_tables.R -i \$INPUT
        software_versions.R "${task.process}" tidyverse data.table optparse dplyr purrr
        """


}

