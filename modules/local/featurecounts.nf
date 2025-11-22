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

        """

}


process MERGE_FEATURECOUNTS {

    label 'process_single'

    // conda '/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/riboseq_nf/riboseq/env.yml'
    container 'iraiosub/nf-riboseq-qc:latest'

    publishDir "${params.outdir}/featurecounts", pattern: "*.featureCounts.tsv.gz", mode: 'copy', overwrite: true

    input:
    path(featurecounts_tables)

    output:
    path("merged.featureCounts.tsv.gz"), emit: counts

    script:

        """
        INPUT=`echo $featurecounts_tables | sed 's/ /,/g'`

        get_featurecounts_tables.R -i \$INPUT
        """


}

