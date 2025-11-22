#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process RIBOCUTTER {

    tag "${sample_id}"
    label 'process_single'

    // conda 'bioconda::ribocutter=0.1.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocutter:0.1.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/ribocutter:0.1.1--pyh5e36f6f_0' }"

    publishDir "${params.outdir}/ribocutter"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.csv"), emit: guides
    // tuple val(sample_id), path("*stats.csv"), emit: stats


    script:

    args = "-g " + params.guide_number
    args += " -r " + params.max_reads
    args += " " + params.extra_ribocutter_args


    if (!params.min_length) {

        min_read_length_arg = ""
        suffix = ""

    } else {

        min_read_length_arg = " --min_read_length " + params.min_length
        suffix = ".min" + params.min_length
    }


        """

        ribocutter -i $reads -o ${sample_id}${suffix} $min_read_length_arg $args

        """
}


process GET_PROPORTION_TARGETED {

    label 'process_single'

    container 'iraiosub/nf-riboseq-qc:latest'

    publishDir "${params.outdir}/ribocutter"

    input:
    path(guides)

    output:
    path("*.pdf"), emit: percentage_targeted
    path("ribocutter_mqc.tsv"), emit: ribocutter_mqc


    script:

        """

        INPUT=`echo $guides | sed 's/ /,/g'`

        echo \$INPUT

        analyse_ribocutter_guides.R -g \$INPUT

        """
}
